%% figure 28

% goal: plot added mass over binned birth size, mirroring Fig 2a of
%       Taheri-Araghi et al. Current Biology (2014)
%

% strategy:
%
%       0. initialize data
%       1. for all conditions in all fluctuating experiments
%       2. gather all added volumes and birth volumes for tracks...
%               - with birth after 3rd hour
%               - with division before bubble appearance (if any)
%       3. bin added volumes by birth volume
%       4. calculate mean and sem of added volumes
%       5. plot over birth volume bin


% versions of this plot:
%       v1. all cell cycles born after 3 hrs, dividing before any bubbles
%       v2. same as v1, but with x axis cropped to 0-10 cubic um
%       v3. all cycle cycles with negative added volumes removed
%       v4. all cell cycles born after FOUR hrs
%       v5. trim births before 3 hrs, normalize all values by condition mean
%           birth size

% last update: jen, 2018 May 22
% commit: attempts to re-produce added size vs birth size plot in Taheri,
%         Curr bio


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

L_birth = cell(1,length(environment_order));
L_added = cell(1,length(environment_order));

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
        lengthVals = conditionData_fullOnly(:,3);       % col 3  = length (um)
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
        
        condition_birthLengths = [];
        condition_addedLengths = [];
        
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
            
            currentLengths = lengthVals(curveFinder == unique_cycles(cc));
            condition_birthLengths(curveCounter,1) = currentLengths(1);
            condition_addedLengths(curveCounter,1) = currentLengths(end)-currentLengths(1);
            
            condition_birthTimes(curveCounter,1) = currentTimes(1);
            condition_divTimes(curveCounter,1) = currentTimes(end);
            
        end
        
        V_birth{eRow,eColumn} = condition_birthVolumes;
        V_added{eRow,eColumn} = condition_addedVolumes;
        
        L_birth{eRow,eColumn} = condition_birthLengths;
        L_added{eRow,eColumn} = condition_addedLengths;
        
        T_birth{eRow,eColumn} = condition_birthTimes;
        T_division{eRow,eColumn} = condition_divTimes;
        
    end
end
%%
save('addedSizeData.mat','V_birth','V_added','L_birth','L_added','T_birth','T_division')

%%

% the above is solid!
% next steps: concatenate all the column data!

clc
clear

% 0. initialize data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
load('addedSizeData.mat')

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
    currCol_birthlengths = L_birth(withData==0,eCol);
    currCol_addedlengths = L_added(withData==0,eCol);
    concat_Lbirth = [];
    concat_Ladded = [];
    
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
        rep_Lbirths = currCol_birthlengths{replicate};
        rep_Ladded = currCol_addedlengths{replicate};
        concat_Lbirth = [concat_Lbirth; rep_Lbirths];
        concat_Ladded = [concat_Ladded; rep_Ladded];
        
        % times
        rep_Tbirths = currCol_birthtimes{replicate};
        rep_Tdivs = currCol_divtimes{replicate};
        concat_Tbirth = [concat_Tbirth; rep_Tbirths];
        concat_Tdiv = [concat_Tdiv; rep_Tdivs];
        
    end
    
    
    % 3. remove all births before 4 hours
%     concat_Vadded = concat_Vadded(concat_Tbirth > 4);
%     concat_Vbirth = concat_Vbirth(concat_Tbirth > 4);
%     
%     concat_Lbirth = concat_Lbirth(concat_Tbirth > 4);
%     concat_Ladded = concat_Ladded(concat_Tbirth > 4);
%     
    
    % 4. remove unphysical events, where added size is negative
    inconceivables_vol = concat_Vadded(concat_Vadded < 0);
    inconceivables_length = concat_Ladded(concat_Ladded < 0);
    length(inconceivables_length)
    
    conceivables_Vbirth = concat_Vbirth(concat_Vadded > 0);
    conceivables_Vadded = concat_Vadded(concat_Vadded > 0);
    
    conceivables_Lbirth = concat_Lbirth(concat_Ladded > 0);
    conceivables_Ladded = concat_Ladded(concat_Ladded > 0);
    
    
    
    % 5. normalize all values by mean of birth size
    mean_Vbirth = mean(conceivables_Vbirth);
    mean_Lbirth = mean(conceivables_Lbirth);
    
    normalized_Vbirth = conceivables_Vbirth./mean_Vbirth;
    normalized_Vadded = conceivables_Vadded./mean_Vbirth;
    
    normalized_Lbirth = conceivables_Lbirth./mean_Lbirth;
    normalized_Ladded = conceivables_Ladded./mean_Lbirth;
    
    
    % 6. bin birth volumes by every 0.1 cubic um 
    %birthBins_vol = floor(conceivables_Vbirth*10);
    birthBins_vol = floor(normalized_Vbirth*10);
    bins = 1:max(birthBins_vol);
    xbins_vol = bins'/ 10;
    
    %birthBins_length = floor(conceivables_Lbirth*10);
    birthBins_length = floor(normalized_Lbirth*10);
    bins = 1:max(birthBins_length);
    xbins_length = bins'/10;
    
    
    % 7. accumulate added volume by birth bin
    %binned_addedVols = accumarray(birthBins_vol,conceivables_Vadded,[],@(x) {x});
    binned_addedVols = accumarray(birthBins_vol,normalized_Vadded,[],@(x) {x});
    
    binned_v_means = cellfun(@mean,binned_addedVols);
    binned_v_stds = cellfun(@std,binned_addedVols);
    binned_v_counts = cellfun(@length,binned_addedVols);
    binned_v_sems = binned_v_stds./sqrt(binned_v_counts);
    
    %binned_addedLengths = accumarray(birthBins_length,conceivables_Ladded,[],@(x) {x});
    binned_addedLengths = accumarray(birthBins_length,normalized_Ladded,[],@(x) {x});
    binned_l_means = cellfun(@mean,binned_addedLengths);
    binned_l_stds = cellfun(@std,binned_addedLengths);
    binned_l_counts = cellfun(@length,binned_addedLengths);
    binned_l_sems = binned_l_stds./sqrt(binned_l_counts);
    
    
    % 8. plot subplots for all conditions
    color = rgb(palette{eCol});
    condition = environment_order{eCol};
    
    % volume
    figure(1)
    subplot(1,length(environment_order),eCol)
    plot(xbins_vol,binned_v_means,'o','Color',color)
    hold on
    errorbar(xbins_vol,binned_v_means,binned_v_sems,'.','Color',color)
    axis([-1 5 -1 4])
    title(condition)
    
    if eCol == 1
        %ylabel('added volume,mean + sem')
        ylabel('added volume / mean birth volume')
    end
    
    if eCol == 4
        %xlabel('birth volume, mean + sem')
        xlabel('birth volume / mean birth volume')
    end
    
    
    % length
    figure(2)
    subplot(1,length(environment_order),eCol)
    plot(xbins_length,binned_l_means,'o','Color',color)
    hold on
    errorbar(xbins_length,binned_l_means,binned_l_sems,'.','Color',color)
    axis([-1 5 -1 4])
    title(condition)
    
    if eCol == 1
        %ylabel('added length, mean + sem')
        ylabel('added length / mean birth length')
    end
    
    if eCol == 4
        %xlabel('birth length, mean + sem')
        xlabel('birth length / mean birth length')
    end
    
end


%%



