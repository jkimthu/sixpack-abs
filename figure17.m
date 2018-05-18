%% figure 17


%  Goal: determine whether growth in low, ave, high and fluc environments
%        are characterizable as adder, sizer, or timer

%        plot (1) single cell added volume vs birth volume
%        plot (2) single cell added volume vs time of birth
%        plot (3) single cell added volume vs time of division


%  Strategy: 
%
%       a) initialize experimental data
%       b) identify complete cell cycles within each condition 
%       c) compile birth volumes and added mass, using isDrop
%       d) calculate mean and plot condition data
%



%  Last edit: jen, 2018 May 18
%  Commit: population averaged added volume vs birth vol across
%          experiments, considering only data after 3 hrs instead of total



%  OK let's go!

%% initialize

clc
clear

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')

dataIndex = find(~cellfun(@isempty,storedMetaData));
experimentCount = length(dataIndex);


% 1. for all experiments in dataset
exptCounter = 0;
for e = 1:experimentCount
    
    % 1. collect experiment meta data
    index = dataIndex(e);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    bubbletime = storedMetaData{index}.bubbletime;
    
    % exclude outliers from analysis (2017-10-31 and monod experiments)
    if strcmp(date, '2017-10-31') == 1 || strcmp (timescale, 'monod') == 1
        disp(strcat(date,': excluded from analysis'))
        continue
    end
    disp(strcat(date, ': analyze!'))
    exptCounter = exptCounter + 1;

    
    % 2. load measured data
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-window5-width1p4-1p7-jiggle-0p5.mat');
    load(filename,'D5','M','M_va','T');
    
    
    % 3. compile experiment data matrix
    xy_start = min(min(storedMetaData{index}.xys));
    xy_end = max(max(storedMetaData{index}.xys));
    exptData = buildDM(D5, M, M_va, T, xy_start, xy_end,e);
    clear D5 M M_va T xy_start xy_end
    
    
    % 4. initialize parameters for plotting
    palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
    shapes = {'o','x','square','*'};
    binsPerHour = 2;
    
    
    for condition = 1:length(palette)
        
        % 5. isolate condition specific data
        conditionData = exptData(exptData(:,23) == condition,:);  % col 23 = cond vals
        
        
        % 6. trim data to full cell cycles ONLY
        ccFraction = conditionData(:,9);            % col 9 = ccFraction
        conditionData_fullOnly = conditionData(~isnan(ccFraction),:);
        
        
        % 7. isolate data to stabilized regions of growth
        minTime = 3;  % hr
        maxTime = bubbletime(condition);
        timestamp = conditionData_fullOnly(:,2)/3600; % col 2   = raw time
        
        times_trim1 = timestamp(timestamp >= minTime);
        conditionData_trim1 = conditionData_fullOnly(timestamp >= minTime,:);
        
        if maxTime > 0
            conditionData_trim2 = conditionData_trim1(times_trim1 <= maxTime,:);
        else
            conditionData_trim2 = conditionData_trim1;
        end
        clear times_trim1 timestamp minTime maxTime
        
        
        % 7. isolate corrected time, added mass, volume, and birth event data (drop)
        curveFinder = conditionData_trim2(:,6);           % col 6   = curve Finder
        addedVol = conditionData_trim2(:,15);             % col 15  = added vol (calculated from Va)
        volume = conditionData_trim2(:,12);               % col 12  = volume (Va)
        timestamps = conditionData_trim2(:,2)/3600;        % col 2   = raw time
        
        
        % 8. for each unique cell cycle, collect V_birth and time
        unique_cycles = unique(curveFinder);
        for cc = 1:length(unique_cycles)
            
            currentTimes = timestamps(curveFinder == unique_cycles(cc));
            timestamps_division(cc,1) = currentTimes(end);
            timestamps_birth(cc,1) = currentTimes(1);
            
            currentVolumes = volume(curveFinder == unique_cycles(cc));
            V_birth(cc,1) = currentVolumes(1);
            V_added(cc,1) = currentVolumes(end) - currentVolumes(1);
            
            
        end
        
        
        % 9. bin added volumes into time bins
        bins_div = ceil(timestamps_division*binsPerHour);
        bins_birth = ceil(timestamps_birth*binsPerHour);
        
        addedVol_binnedByDivTime_mean = accumarray(bins_div,V_added,[],@mean);
        addedVol_binnedByDivTime_std = accumarray(bins_div,V_added,[],@std);
        addedVol_binnedByDivTime_counts = accumarray(bins_div,V_added,[],@length);
        addedVol_binnedByDivTime_sem = addedVol_binnedByDivTime_std./sqrt(addedVol_binnedByDivTime_counts);
        
        addedVol_binnedByBirthTime_mean = accumarray(bins_birth,V_added,[],@mean);
        addedVol_binnedByBirthTime_std = accumarray(bins_birth,V_added,[],@std);
        addedVol_binnedByBirthTime_counts = accumarray(bins_birth,V_added,[],@length);
        addedVol_binnedByBirthTime_sem = addedVol_binnedByBirthTime_std./sqrt(addedVol_binnedByBirthTime_counts);
        
        timeVector = 1:(10*binsPerHour);
        if length(addedVol_binnedByDivTime_mean) < length(timeVector)
            
            addedVol_binnedByDivTime_mean = [addedVol_binnedByDivTime_mean; zeros(length(timeVector)-length(addedVol_binnedByDivTime_mean),1)];
            addedVol_binnedByDivTime_std = [addedVol_binnedByDivTime_std; zeros(length(timeVector)-length(addedVol_binnedByDivTime_std),1)];
            addedVol_binnedByDivTime_counts = [addedVol_binnedByDivTime_counts; zeros(length(timeVector)-length(addedVol_binnedByDivTime_counts),1)];
            addedVol_binnedByDivTime_sem = [addedVol_binnedByDivTime_sem; zeros(length(timeVector)-length(addedVol_binnedByDivTime_sem),1)];
            
        end
        
        if length(addedVol_binnedByBirthTime_mean) < length(timeVector)
            
            addedVol_binnedByBirthTime_mean = [addedVol_binnedByBirthTime_mean; zeros(length(timeVector)-length(addedVol_binnedByBirthTime_mean),1)];
            addedVol_binnedByBirthTime_std = [addedVol_binnedByBirthTime_std; zeros(length(timeVector)-length(addedVol_binnedByBirthTime_std),1)];
            addedVol_binnedByBirthTime_counts = [addedVol_binnedByBirthTime_counts; zeros(length(timeVector)-length(addedVol_binnedByBirthTime_counts),1)];
            addedVol_binnedByBirthTime_sem = [addedVol_binnedByBirthTime_sem; zeros(length(timeVector)-length(addedVol_binnedByBirthTime_sem),1)];
            
        end
            
        % 10. plot
        color = rgb(palette(condition));
        
        if condition == 1 && timescale == 300
            xmark = shapes{2};
        elseif condition == 1 && timescale == 900
            xmark = shapes{3};
        elseif condition == 1 && timescale == 3600
            xmark = shapes{4};
        else
            xmark = shapes{1};
        end
        
        % added volume vs. birth size
        figure(17)
        errorbar(mean(V_birth),mean(V_added),std(V_added),'Color',color)
        hold on
        plot(mean(V_birth),mean(V_added),'Marker',xmark,'Color',color)
        hold on
        xlabel('birth volume (cubic um)')
        ylabel('added volume (cubic um)')
        title('population averaged added mass vs birth size')
        axis([0 12 0 12])
        
        
        % added volume vs. time of birth
        figure(18)
        errorbar(timeVector,addedVol_binnedByBirthTime_mean,addedVol_binnedByBirthTime_sem,'Color',color)
        hold on
        plot(timeVector,addedVol_binnedByBirthTime_mean,'Marker',xmark,'Color',color)
        hold on
        xlabel('time of birth (hr)')
        ylabel('added volume (cubic um)')
        title('population-averaged added volume binned by time of birth')
        axis([0 timeVector(end)+1 0 12])
        
        % added volume vs. time of division
        figure(19)
        errorbar(timeVector,addedVol_binnedByDivTime_mean,addedVol_binnedByDivTime_sem,'Color',color)
        hold on
        plot(timeVector,addedVol_binnedByDivTime_mean,'Marker',xmark,'Color',color)
        hold on
        xlabel('time of division (hr)')
        ylabel('added volume (cubic um)')
        title('population-averaged added volume binned by time of division')
        axis([0 timeVector(end)+1 0 12])
        
        clear timestamps_birth timestamps_division currentTimes
        
    end
    
end
