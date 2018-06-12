%% figure 38

% goal: like figure 38, mirroring Fig 2a of Taheri-Araghi et al. Current Biology (2014)
%       but added width instead of added length or volume


% strategy:
%
%       0. initialize data
%       1. for all conditions in all fluctuating experiments
%       2. gather all added widths and birth widths for tracks...
%               - with birth after 3rd hour
%               - with division before bubble appearance (if any)
%       3. bin added width by birth width
%       4. calculate mean and sem of added width
%       5. plot over birth width bin


% criteria for passable cell cycles:
%       1. all cell cycles born after 3 hrs, dividing before any bubbles
%       2. all cycle cycles with negative added VOLUME removed
%       3. remove cell cycles with V_division greater than 3x V_birth

% last update: jen, 2018 Jun 12

% commit: Taheri-like adder plot with added width

% OK let's go!

%%

clc
clear

% 0. initialize data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')

dataIndex = find(~cellfun(@isempty,storedMetaData));
experimentCount = length(dataIndex);


% 1. initialize environmental conditions for data collection and plotting
palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
environment_order = {'low','30','300','900','3600','ave','high'};
environment_ticks = zeros(length(environment_order),1);

V_birth = cell(1,length(environment_order));
V_added = cell(1,length(environment_order));

W_birth = cell(1,length(environment_order));
W_added = cell(1,length(environment_order));

T_birth = cell(1,length(environment_order));
T_division = cell(1,length(environment_order));

%%
% for all experiments
for e = 1:experimentCount
    
    % 2. collect experiment meta data
    index = dataIndex(e);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    bubbletime = storedMetaData{index}.bubbletime;
    
    % excluding monod experiments and outliers
    if strcmp(date, '2017-10-31') == 1 || strcmp (timescale, 'monod') == 1
        disp(strcat(date,': excluded from analysis'))
        continue
    end
    disp(strcat(date, ': analyze!'))
    
    
    % 3. load measured data
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-window5-width1p4-1p7-jiggle-0p5.mat');
    load(filename,'D5','M','M_va','T');
    
    
    % 4. compile experiment data matrix
    xy_start = min(min(storedMetaData{index}.xys));
    xy_end = max(max(storedMetaData{index}.xys));
    exptData = buildDM(D5, M, M_va, T, xy_start, xy_end,e);
    
    
    for condition = 1:length(bubbletime)
        
        % 5. isolate condition specific data
        conditionData = exptData(exptData(:,23) == condition,:);  % col 23 = cond vals
        
        
        % 6. trim data to full cell cycles ONLY
        ccFraction = conditionData(:,9);            % col 9 = ccFraction
        conditionData_fullOnly = conditionData(~isnan(ccFraction),:);
        
        
        % 7. isolate corrected time, added mass, volume, and birth event data (drop)
        curveFinder = conditionData_fullOnly(:,6);      % col 6   = curve Finder
        widths = conditionData_fullOnly(:,11);          % col 11  = width (um)
        volumes = conditionData_fullOnly(:,12);         % col 12  = volume (Va)
        isBirth = conditionData_fullOnly(:,5);          % col 5   = isDrop (1 = birth event, 0 = not)
        timestamps = conditionData_fullOnly(:,2)/3600;  % col 2   = raw timestamps
        
        
        % 8. prepare to assign condition data into a cell, where:
        % column = environmental condition
        % row = biological replicate
        
        % i. determine column no. of environmental condition
        if condition == 2
            eColumn = find(strcmp(environment_order,'low'));
        elseif condition == 3
            eColumn = find(strcmp(environment_order,'ave'));
        elseif condition == 4
            eColumn = find(strcmp(environment_order,'high'));
        else
            eColumn = find(strcmp(environment_order,num2str(timescale)));
        end
        environment_ticks(eColumn) = environment_ticks(eColumn) + 1;
        
        % ii. determine replicate no. of current condition data
        eRow = environment_ticks(eColumn);
        
        
        
        % 9. for each unique cell cycle, collect birth size, added size, ad
        %    time of birth and division
        unique_cycles = unique(curveFinder);
        curveCounter = 0;
        
        condition_birthVolumes = [];
        condition_addedVolumes = [];
        
        condition_birthWidths = [];
        condition_addedWidths = [];
        
        condition_birthTimes = [];
        condition_divTimes = [];
        
        
        for cc = 1:length(unique_cycles)
            
            currentTimes = timestamps(curveFinder == unique_cycles(cc));
            
            % discard all cell cycles shorter than 10 min
            if length(currentTimes) < 5
                disp(strcat('short cycle: ',num2str(length(currentTimes))))
                continue
            end
            
            % discard all cell cycles born before 3 hr or dividing after bubble
            if currentTimes(1) < 3
                disp(strcat(num2str(currentTimes(1)),': toss, before 3 hrs'))
                continue
            elseif bubbletime(condition) ~= 0 && currentTimes(end) > bubbletime(condition)
                disp(strcat(num2str(currentTimes(end)),': toss, divides after bubble'))
                continue
            end
            
            % for all remaining curves, count and collect data
            curveCounter = curveCounter + 1;
            
            currentVolumes = volumes(curveFinder == unique_cycles(cc));
            condition_birthVolumes(curveCounter,1) = currentVolumes(1);
            condition_addedVolumes(curveCounter,1) = currentVolumes(end)-currentVolumes(1);
            
            currentWidths = widths(curveFinder == unique_cycles(cc));
            condition_birthWidths(curveCounter,1) = currentWidths(1);
            condition_addedWidths(curveCounter,1) = currentWidths(end)-currentWidths(1);
            
            condition_birthTimes(curveCounter,1) = currentTimes(1);
            condition_divTimes(curveCounter,1) = currentTimes(end);
            
        end
        
        V_birth{eRow,eColumn} = condition_birthVolumes;
        V_added{eRow,eColumn} = condition_addedVolumes;
        
        W_birth{eRow,eColumn} = condition_birthWidths;
        W_added{eRow,eColumn} = condition_addedWidths;
        
        T_birth{eRow,eColumn} = condition_birthTimes;
        T_division{eRow,eColumn} = condition_divTimes;
        
    end
end
%%
save('addedSizeData_width.mat','V_birth','V_added','W_birth','W_added','T_birth','T_division')

%%

% the above is solid!
% next steps: concatenate all the column data!

clc
clear

% 0. initialize data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
load('addedSizeData_width.mat')

palette = {'Indigo','DarkTurquoise','SteelBlue','DeepSkyBlue','DodgerBlue','GoldenRod','FireBrick'};
environment_order = {'low','30','300','900','3600','ave','high'};

% 1. for each environment... 
for eCol = 1:length(environment_order)
    
    % 2. isolate column and concatenate data
    withData = cellfun(@isempty,V_birth(:,eCol));
    
    % volume
    currCol_birthvols = V_birth(withData==0,eCol);
    currCol_addedvols = V_added(withData==0,eCol);
    concat_Vbirth = [];
    concat_Vadded = [];
    
    % length
    currCol_birthwidths = W_birth(withData==0,eCol);
    currCol_addedwidths = W_added(withData==0,eCol);
    concat_Wbirth = [];
    concat_Wadded = [];
    
    % times
    currCol_birthtimes = T_birth(withData==0,eCol);
    currCol_divtimes = T_division(withData==0,eCol);
    concat_Tbirth = [];
    concat_Tdiv = [];
    
    for replicate = 1:length(currCol_birthvols)
        
        % volume
        rep_Vbirths = currCol_birthvols{replicate};
        rep_Vadded = currCol_addedvols{replicate};
        concat_Vbirth = [concat_Vbirth; rep_Vbirths];
        concat_Vadded = [concat_Vadded; rep_Vadded];
        
        % length
        rep_Wbirths = currCol_birthwidths{replicate};
        rep_Wadded = currCol_addedwidths{replicate};
        concat_Wbirth = [concat_Wbirth; rep_Wbirths];
        concat_Wadded = [concat_Wadded; rep_Wadded];
        
        % times
        rep_Tbirths = currCol_birthtimes{replicate};
        rep_Tdivs = currCol_divtimes{replicate};
        concat_Tbirth = [concat_Tbirth; rep_Tbirths];
        concat_Tdiv = [concat_Tdiv; rep_Tdivs];
        
    end
    
   
    
    % 4. remove unphysical events, where added size is negative
    inconceivables_vol = concat_Vadded(concat_Vadded < 0);
    inconceivables_width = concat_Wadded(concat_Vadded < 0);
    
    conceivables_Vbirth = concat_Vbirth(concat_Vadded > 0);
    conceivables_Vadded = concat_Vadded(concat_Vadded > 0);
    
    conceivables_Wbirth = concat_Wbirth(concat_Vadded > 0);
    conceivables_Wadded = concat_Wadded(concat_Vadded > 0);
    
    
    % 5. remove non-single cell cycles, where division size > 3* birth size
    division_volume = conceivables_Vbirth + conceivables_Vadded;
    division_width = conceivables_Wbirth + conceivables_Wadded;
    
    vBirth_3x = conceivables_Vbirth * 3;
    
    nonSingles_volume = length(division_volume(division_volume > vBirth_3x))
    
    singles_Vbirth = conceivables_Vbirth(division_volume <= vBirth_3x);
    singles_Vadded = conceivables_Vadded(division_volume <= vBirth_3x);
    
    singles_Wbirth = conceivables_Wbirth(division_volume <= vBirth_3x);
    singles_Wadded = conceivables_Wadded(division_volume <= vBirth_3x);
    
    
    
    % 5. normalize all values by mean of birth size
    mean_Vbirth = mean(singles_Vbirth);
    mean_Wbirth = mean(singles_Wbirth);
    
    normalized_Vbirth = singles_Vbirth./mean_Vbirth;
    normalized_Vadded = singles_Vadded./mean_Vbirth;
    
    normalized_Wbirth = singles_Wbirth./mean_Wbirth;
    normalized_Wadded = singles_Wadded./mean_Wbirth;
    
    
    % 6. bin birth volumes by every 0.1 cubic um 
    birthBins_vol = floor(conceivables_Vbirth*10);
    bins = 1:max(birthBins_vol);
    xbins_vol = bins'/ 10;
    
    birthBins_vol_normalized = floor(normalized_Vbirth*10);
    bins_normalized = 1:max(birthBins_vol_normalized);
    xbins_vol_normalized = bins_normalized'/ 10;
    
    birthBins_width = floor(conceivables_Wbirth*10);
    bins = 1:max(birthBins_width);
    xbins_width = bins'/10;
    
    birthBins_width_normalized = floor(normalized_Wbirth*10);
    bins_normalized = 1:max(birthBins_width_normalized);
    xbins_width_normalized = bins_normalized'/10;
    
    
    % 7. accumulate added volume by birth bin
    binned_addedVols = accumarray(birthBins_vol,conceivables_Vadded,[],@(x) {x});
    binned_v_means = cellfun(@mean,binned_addedVols);
    binned_v_stds = cellfun(@std,binned_addedVols);
    binned_v_counts = cellfun(@length,binned_addedVols);
    binned_v_sems = binned_v_stds./sqrt(binned_v_counts);
    
    binned_addedVols_normalized = accumarray(birthBins_vol_normalized,normalized_Vadded,[],@(x) {x});
    binned_v_norm_means = cellfun(@mean,binned_addedVols_normalized);
    binned_v_norm_stds = cellfun(@std,binned_addedVols_normalized);
    binned_v_norm_counts = cellfun(@length,binned_addedVols_normalized);
    binned_v_norm_sems = binned_v_norm_stds./sqrt(binned_v_norm_counts);
    
    binned_addedWidths = accumarray(birthBins_width,conceivables_Wadded,[],@(x) {x});
    binned_w_means = cellfun(@mean,binned_addedWidths);
    binned_w_stds = cellfun(@std,binned_addedWidths);
    binned_w_counts = cellfun(@length,binned_addedWidths);
    binned_w_sems = binned_w_stds./sqrt(binned_w_counts);
    
    binned_addedWidths_normalized = accumarray(birthBins_width_normalized,normalized_Wadded,[],@(x) {x});
    binned_w_norm_means = cellfun(@mean,binned_addedWidths_normalized);
    binned_w_norm_stds = cellfun(@std,binned_addedWidths_normalized);
    binned_w_norm_counts = cellfun(@length,binned_addedWidths_normalized);
    binned_w_norm_sems = binned_w_norm_stds./sqrt(binned_w_norm_counts);
    
    
    % 8. plot subplots for all conditions
    color = rgb(palette{eCol});
    condition = environment_order{eCol};
    
    
    % volume, not normalized
    figure(1)
    subplot(1,length(environment_order),eCol)
    plot(xbins_vol,binned_v_means,'o','Color',color)
    hold on
    errorbar(xbins_vol,binned_v_means,binned_v_sems,'.','Color',color)
    %axis([0 10 -2 10])
    axis([0 5 -1 6])
    title(condition)
    if eCol == 1
        ylabel('added volume,mean + sem')
    end
    if eCol == 4
        xlabel('birth volume, mean + sem')
    end
    
    % volume, normalized
    figure(2)
    subplot(1,length(environment_order),eCol)
    plot(xbins_vol_normalized,binned_v_norm_means,'o','Color',color)
    hold on
    errorbar(xbins_vol_normalized,binned_v_norm_means,binned_v_norm_sems,'.','Color',color)
    %axis([-1 5 -1 4])
    axis([0 3 -.1 3])
    title(condition)
    if eCol == 1
        ylabel('added volume / mean birth volume')
    end
    if eCol == 4
        xlabel('birth volume / mean birth volume')
    end
    
    
    % width, not normalized
    figure(3)
    subplot(1,length(environment_order),eCol)
    plot(xbins_width,binned_w_means,'o','Color',color)
    hold on
    errorbar(xbins_width,binned_w_means,binned_w_sems,'.','Color',color)
    %axis([0 10 -2 10])
    axis([0 2 -.5 .5])
    title(condition)
    if eCol == 1
        ylabel('added width, mean + sem')
    end
    if eCol == 4
        xlabel('birth width, mean + sem')
    end
    
    
    % length, normalized
    figure(4)
    subplot(1,length(environment_order),eCol)
    plot(xbins_width_normalized,binned_w_norm_means,'o','Color',color)
    hold on
    errorbar(xbins_width_normalized,binned_w_norm_means,binned_w_norm_sems,'.','Color',color)
    %axis([-1 5 -1 4])
    axis([0 2 -.5 .5])
    title(condition)
    if eCol == 1
        ylabel('added width / mean birth width')
    end
    if eCol == 4
        xlabel('birth width / mean birth width')
    end
    
end


%%



