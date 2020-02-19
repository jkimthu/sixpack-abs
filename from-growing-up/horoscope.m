%% horoscope


% Goal: are growth phenotypes predicted by when a cell is born? plot multiple stats after binning birth data by nutrient phase
%
%       Mean and scatter plots:
%          
%           1. cell cycle duration
%           2. volume added per cell cycle
%           3. fraction of cell cycle spent in high nutrient


%  Last edit: Jen Nguyen, 2018 Feb 27

%  Commit: revise nScore plots using updated nutrientScore, version 2.
%  trimming out data from cell cycles less than 10 mins long also makes
%  signals cleaner


%  Strategy:
%
%     0. initialize complete meta data
%     0. for experiments of interest...
%           1. collect experiment date and exclude outliers (2017-10-31 and monod experiments)
%           2. initialize binning parameters
%           3. load measured data for stable condition
%           4. for stable average... 
%                  5. isolate isDrop and timestamp data
%                  6. accumulate and trim data to stabilized region based on timestamp
%                  7. isolate birth events (isDrop == 1), corresponding data and timestamps
%                  8. re-define period to begin at start of low nutrient pulse, by subtracting quarter period from corrected timestamp
%                  9. assign elements of timestamp vector to period fraction bins
%                 10. bin cell cycle durations (greater than 10 min) by period fraction
%                 11. plot mean and s.e.m. of cell cycle duration over nutrient period 
%                 12. plot scatter of cell cycle duration over nutrient period 
%          13. repeat analysis for fluctuating environment, plotting fluc data over stable


% OK! Lez go!


%% ONE. CELL CYCLE DURATION

% 0. initialize data

clc
clear

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));
experimentCount = length(dataIndex);
                            

% for target experiments...
for e = 14
    
    % 1. collect experiment date
    index = dataIndex(e);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    bubbleTime = storedMetaData{index}.bubbletime;
    
    % exclude outliers from analysis (2017-10-31 and monod experiments)
    if strcmp(date, '2017-10-31') == 1 || strcmp (timescale, 'monod') == 1
       disp(strcat(date,': excluded from analysis'))
       continue
    end
    disp(strcat(date, ': analyze!'))

    
    % 2. initialize binning parameters
    binsPerPeriod = 12;
    
    
    % 3. load measured data in stable average condition
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-window5-width1p4-1p7-jiggle-0p5.mat');
    load(filename,'D','D5','M','M_va','T');
    
    % 4. for stable average and then fluctuating environment...
    environment = [3;1];
    
    for i = 1:length(environment)
        
        condition = environment(i);
        xy_start = storedMetaData{index}.xys(condition,1);
        xy_end = storedMetaData{index}.xys(condition,end);
        conditionData = buildDM(D5, M, M_va, T, xy_start, xy_end, e);
        
        % 5. isolate relevant data
        isDrops = conditionData(:,5);           % col 5  =  isDrop; 0 during curve, 1 at birth event
        durations = conditionData(:,8);         % col 8  =  curve durations
        correctedTimes = conditionData(:,30);   % col 30 =  true times after lag correction
            

        % 6. accumulate data for easy trimming
        data = [isDrops correctedTimes durations];
        
        % 6i. trim data to timepoints after 3 hrs 
        data_trim1 = data(correctedTimes/3600 >= 3,:);
         
        
        % 6ii. trim data to timepoints before appearance of bubbles
        if bubbleTime(condition) > 0
            
            data_trim2 = data_trim1(data_trim1(:,2)/3600 <= bubbleTime(condition),:);
            
        else
            
            data_trim2 = data_trim1;
                       
        end
        
        % 6iii. recover data as individual
        isDrops_trimmed = data_trim2(:,1);
        
        
        % 7. isolate birth events (isDrop == 1), corresponding data and timestamps
        birthEvents = isDrops_trimmed(isDrops_trimmed == 1,1);
        birthEvent_timestamps = data_trim2(isDrops_trimmed == 1,2); % uses corrected timestamps
        birthEvent_durations = data_trim2(isDrops_trimmed == 1,3)/60; % min

        
        % 8.  re-define period to begin at start of low nutrient pulse, by
        %     subtracting quarter period from corrected timestamp
        shifted_birthTimestamps = birthEvent_timestamps - (timescale/4);
        
        % 9. assign elements of timestamp vector to period fraction bins
        timeInPeriods = shifted_birthTimestamps/timescale; % unit = sec/sec
        timeInPeriodFraction = timeInPeriods - floor(timeInPeriods);
        assignedBin = ceil(timeInPeriodFraction * binsPerPeriod);
        
        % 10. bin cell cycle durations by period fraction
        % 10. (i) first remove all durations of lesser than 10 min
        birthEvent_durations_trimmedByDuration = birthEvent_durations(birthEvent_durations > 10);
        assignedBin_trimmedByDuration = assignedBin(birthEvent_durations > 10);
        
        % 10. (ii) bin
        durations_binnedByPeriodFraction = accumarray(assignedBin_trimmedByDuration, birthEvent_durations_trimmedByDuration, [], @(x) {x});
        
        durations_binned_mean = cellfun(@mean,durations_binnedByPeriodFraction);
        durations_binned_std = cellfun(@std,durations_binnedByPeriodFraction);
        durations_binned_count = cellfun(@length,durations_binnedByPeriodFraction);
        durations_binned_sem = durations_binned_std./durations_binned_count;
        
        
        % 11. plot mean and s.e.m. of cell cycle duration over nutrient period 
        % 11. (i) convert bin # to absolute time (sec)
        timePerBin = timescale/binsPerPeriod;  % in sec
        binPeriod = linspace(1, binsPerPeriod, binsPerPeriod);
        timePeriod = timePerBin*binPeriod';

        % 11. (ii) repeat quarter period on both sides and plot over period fraction
        quarterOne = linspace(1,binsPerPeriod/4,binsPerPeriod/4);
        quarterFour = linspace(binsPerPeriod*3/4+1,binsPerPeriod,binsPerPeriod/4);
        quarterZero = linspace((binsPerPeriod/4-1),0,binsPerPeriod/4)*-1;
        
        signal_quarterOne = durations_binned_mean(quarterOne,1);
        signal_quarterFour = durations_binned_mean(quarterFour,1);

        time_quarterOne = timePeriod(quarterOne,1) + timescale;
        time_quarterZero = quarterZero*timePerBin;
        
        mean_durationSignal_stitched = [signal_quarterFour; durations_binned_mean; signal_quarterOne];
        time_stitched = [time_quarterZero'; timePeriod; time_quarterOne];
  
        figure(1)
        if condition == 3
            subplot(1,2,1)
            stem(time_stitched, mean_durationSignal_stitched,'Color',[0.25 0.25 0.9])
            hold on
            axis([min(time_stitched) max(time_stitched) 0 90])
            title(strcat('stable: ',date))
        elseif condition == 1
            subplot(1,2,2)
            stem(time_stitched, mean_durationSignal_stitched,'Color',[0 0.7 0.7])
            axis([min(time_stitched) max(time_stitched) 0 90])
            title(strcat('fluc: ',num2str(timescale)))
        end
        xlabel('time (sec)')
        ylabel('mean cell cycle duration (min)')
        grid on
       
        
        % 12. plot scatter of cell cycle duration over nutrient period 
        scatter_quarterOne = durations_binnedByPeriodFraction(quarterOne,1);
        scatter_quarterFour = durations_binnedByPeriodFraction(quarterFour,1);
        scatter_stitched = [scatter_quarterFour; durations_binnedByPeriodFraction; scatter_quarterOne];
        
        sz = 20; % circle size
        figure(2)
        if condition == 3
            subplot(1,2,1)
            for b = 1:length(scatter_stitched)
                y = scatter_stitched{b};
                x = ones(length(y),1)*time_stitched(b);
                scatter(x,y,sz,[0.25 0.25 0.9]) % purple
                hold on
            end
            axis([min(time_stitched) max(time_stitched) 0 90])
            title(strcat('stable: ',date))
        else
            subplot(1,2,2)
            for b = 1:length(scatter_stitched)
                y = scatter_stitched{b};
                x = ones(length(y),1)*time_stitched(b);
                scatter(x,y,sz,[0 0.7 0.7]) % green
                hold on
            end
            axis([min(time_stitched) max(time_stitched) 0 90])
            title(strcat('fluc: ',num2str(timescale)))
        end
        xlabel('time (sec)')
        ylabel('cell cycle duration (min)')
        grid on
        
        % 13. repeat analysis for fluctuating environment, plotting fluc data over stable
    end
    clear i environment 
    
    cd('/Users/jen/Documents/StockerLab/Data_analysis/horoscope')
    
    Fig1 = figure(1);
    saveas(Fig1,strcat('horoscope-duration-mean-',num2str(timescale),'-',date),'epsc')
    close(Fig1)
    
    Fig2 = figure(2);
    saveas(Fig2,strcat('horoscope-duration-scatter-',num2str(timescale),'-',date),'epsc')
    close(Fig2)
    
end


%% TWO. VOLUME ADDED

% 0. initialize data

clc
clear

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));
experimentCount = length(dataIndex);
                            

% for target experiments...
for e = 6:13
    
    % 1. collect experiment date
    index = dataIndex(e);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    bubbleTime = storedMetaData{index}.bubbletime;
    
    % exclude outliers from analysis (2017-10-31 and monod experiments)
    if strcmp(date, '2017-10-31') == 1 || strcmp (timescale, 'monod') == 1
       disp(strcat(date,': excluded from analysis'))
       continue
    end
    disp(strcat(date, ': analyze!'))

    
    % 2. initialize binning parameters
    binsPerPeriod = 12;
    
    
    % 3. load measured data in stable average condition
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-window5-width1p4-1p7-jiggle-0p5.mat');
    load(filename,'D','D5','M','M_va','T');
    
    % 4. for stable average and then fluctuating environment...
    environment = [3;1];
    
    for i = 1:length(environment)
        
        condition = environment(i);
        xy_start = storedMetaData{index}.xys(condition,1);
        xy_end = storedMetaData{index}.xys(condition,end);
        conditionData = buildDM(D5, M, M_va, T, xy_start, xy_end);
        
        % 5. isolate relevant data
        isDrops = conditionData(:,5);           % col 5  =  isDrop; 0 during curve, 1 at birth event
        vaAdded = conditionData(:,20);          % col 20 =  volume (cylinder with caps) added over cell cycle
        correctedTimes = conditionData(:,30);   % col 30 =  timestamps corrected for lag between junc and xy positions
        
                  
        % 6. accumulate data for easy trimming
        data = [isDrops correctedTimes vaAdded];
        
        % 6i. trim data to timepoints after 3 hrs 
        data_trim1 = data(correctedTimes/3600 >= 3,:);
         
        
        % 6ii. trim data to timepoints before appearance of bubbles
        if bubbleTime(condition) > 0
            
            data_trim2 = data_trim1(data_trim1(:,2)/3600 <= bubbleTime(condition),:);
            
        else
            
            data_trim2 = data_trim1;
                       
        end
        
        % 6iii. recover data as individual
        isDrops_trimmed = data_trim2(:,1);
        
        
        % 7. isolate birth events (isDrop == 1), corresponding data and timestamps
        birthEvents = isDrops_trimmed(isDrops_trimmed == 1,1);
        birthEvent_timestamps = data_trim2(isDrops_trimmed == 1,2);
        birthEvent_addedVol = data_trim2(isDrops_trimmed == 1,3);
       
   
        % 8.  re-define period to begin at start of low nutrient pulse, by
        %     subtracting quarter period from corrected timestamp
        shifted_birthTimestamps = birthEvent_timestamps - (timescale/4);
        
        % 9. assign elements of timestamp vector to period fraction bins
        timeInPeriods = shifted_birthTimestamps/timescale; % unit = sec/sec
        timeInPeriodFraction = timeInPeriods - floor(timeInPeriods);
        assignedBin = ceil(timeInPeriodFraction * binsPerPeriod);
        
        % 10. bin addedVolumes by period fraction
        % 10. (i) first remove all zero-valued added volumes
        addedVolumes_nonZeros = birthEvent_addedVol(birthEvent_addedVol > 0);
        assignedBin_nonZeros = assignedBin(birthEvent_addedVol > 0);
        
        % 10. (ii) bin
        addedVolumes_binnedByPeriodFraction = accumarray(assignedBin_nonZeros, addedVolumes_nonZeros, [], @(x) {x});
        
        addedVol_binned_mean = cellfun(@mean,addedVolumes_binnedByPeriodFraction);
        addedVol_binned_std = cellfun(@std,addedVolumes_binnedByPeriodFraction);
        addedVol_binned_count = cellfun(@length,addedVolumes_binnedByPeriodFraction);
        addedVol_binned_sem = addedVol_binned_std./addedVol_binned_count;
        
        
        % 11. plot mean and s.e.m. of cell cycle duration over nutrient period 
        % 11. (i) convert bin # to absolute time (sec)
        timePerBin = timescale/binsPerPeriod;  % in sec
        binPeriod = linspace(1, binsPerPeriod, binsPerPeriod);
        timePeriod = timePerBin*binPeriod';

        % 11. (ii) repeat quarter period on both sides and plot over period fraction
        quarterOne = linspace(1,binsPerPeriod/4,binsPerPeriod/4);
        quarterFour = linspace(binsPerPeriod*3/4+1,binsPerPeriod,binsPerPeriod/4);
        quarterZero = linspace((binsPerPeriod/4-1),0,binsPerPeriod/4)*-1;
        
        signal_quarterOne = addedVol_binned_mean(quarterOne,1);
        signal_quarterFour = addedVol_binned_mean(quarterFour,1);

        time_quarterOne = timePeriod(quarterOne,1) + timescale;
        time_quarterZero = quarterZero*timePerBin;
        
        mean_addedSignal_stitched = [signal_quarterFour; addedVol_binned_mean; signal_quarterOne];
        time_stitched = [time_quarterZero'; timePeriod; time_quarterOne];
  
        figure(1)
        if condition == 3
            subplot(1,2,1)
            stem(time_stitched, mean_addedSignal_stitched,'Color',[0.25 0.25 0.9])
            hold on
            axis([min(time_stitched) max(time_stitched) 0 10])
            title(strcat('stable: ',date))
        elseif condition == 1
            subplot(1,2,2)
            stem(time_stitched, mean_addedSignal_stitched,'Color',[0 0.7 0.7])
            axis([min(time_stitched) max(time_stitched) 0 10])
            title(strcat('fluc: ',num2str(timescale)))
        end
        xlabel('time (sec)')
        ylabel('mean added vol (cubic um)')
        grid on
       
        
        % 12. plot scatter of cell cycle duration over nutrient period 
        scatter_quarterOne = addedVolumes_binnedByPeriodFraction(quarterOne,1);
        scatter_quarterFour = addedVolumes_binnedByPeriodFraction(quarterFour,1);
        scatter_stitched = [scatter_quarterFour; addedVolumes_binnedByPeriodFraction; scatter_quarterOne];
        
        sz = 20; % circle size
        figure(2)
        if condition == 3
            subplot(1,2,1)
            for b = 1:length(scatter_stitched)
                y = scatter_stitched{b};
                x = ones(length(y),1)*time_stitched(b);
                scatter(x,y,sz,[0.25 0.25 0.9]) % purple
                hold on
            end
            axis([min(time_stitched) max(time_stitched) 0 10])
            title(strcat('stable: ',date))
        else
            subplot(1,2,2)
            for b = 1:length(scatter_stitched)
                y = scatter_stitched{b};
                x = ones(length(y),1)*time_stitched(b);
                scatter(x,y,sz,[0 0.7 0.7]) % green
                hold on
            end
            axis([min(time_stitched) max(time_stitched) 0 10])
            title(strcat('fluc: ',num2str(timescale)))
        end
        xlabel('time (sec)')
        ylabel('added volume (cubic um)')
        grid on
        
        % 13. repeat analysis for fluctuating environment, plotting fluc data over stable
    end
    clear i environment 
    
    cd('/Users/jen/Documents/StockerLab/Data_analysis/horoscope')
    
    Fig1 = figure(1);
    saveas(Fig1,strcat('horoscope-addedVol-mean-',num2str(timescale),'-',date),'epsc')
    close(Fig1)
    
    Fig2 = figure(2);
    saveas(Fig2,strcat('horoscope-addedVol-scatter-',num2str(timescale),'-',date),'epsc')
    close(Fig2)
    
end

%% THREE. nSCORE


% 0. initialize data

clc
clear

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));
experimentCount = length(dataIndex);
                            

% for target experiments...
for e = 6:14
    
    % 1. collect experiment date
    index = dataIndex(e);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    bubbleTime = storedMetaData{index}.bubbletime;
    
    % exclude outliers from analysis (2017-10-31 and monod experiments)
    if strcmp(date, '2017-10-31') == 1 || strcmp (timescale, 'monod') == 1
       disp(strcat(date,': excluded from analysis'))
       continue
    end
    disp(strcat(date, ': analyze!'))

    
    % 2. initialize binning parameters
    binsPerPeriod = 12;
    
    
    % 3. load measured data in stable average condition
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-window5-width1p4-1p7-jiggle-0p5.mat');
    load(filename,'D','D5','M','M_va','T');
    
    % 4. for stable average and then fluctuating environment...
    environment = [3;1];
    
    for i = 1:length(environment)
        
        condition = environment(i);
        xy_start = storedMetaData{index}.xys(condition,1);
        xy_end = storedMetaData{index}.xys(condition,end);
        conditionData = buildDM(D5, M, M_va, T, xy_start, xy_end, e);
        
        % 5. isolate relevant data
        isDrops = conditionData(:,5);           % col 5  =  isDrop; 0 during curve, 1 at birth event
        correctedTimes = conditionData(:,30);   % col 30 =  true times after lag correction
        durations = conditionData(:,8)/60;      % col 8  =  curve durations
        
        % 6. calculate nScores for condition
        [binaryNutrientSignal, nScore] = nutrientScore(timescale,conditionData);

        % 6. accumulate data for easy trimming
        data = [isDrops correctedTimes nScore durations];
        
        % 6i. trim data to timepoints after 3 hrs 
        data_trim1 = data(correctedTimes/3600 >= 3,:);
         
        
        % 6ii. trim data to timepoints before appearance of bubbles
        if bubbleTime(condition) > 0
            
            data_trim2 = data_trim1(data_trim1(:,2)/3600 <= bubbleTime(condition),:);
            
        else
            
            data_trim2 = data_trim1;
                       
        end
        
        % 6iii. recover data as individual
        isDrops_trimmed = data_trim2(:,1);
        correctedTimes_trimmed = data_trim2(:,2);
        nScores_trimmed = data_trim2(:,3);
        durations_trimmed = data_trim2(:,4);     % min
        
        
        % 7. isolate birth events (isDrop == 1), corresponding data and timestamps
        birthEvents = isDrops_trimmed(isDrops_trimmed == 1);
        birthEvents_timestamps = correctedTimes_trimmed(isDrops_trimmed == 1);
        birthEvents_nScores = nScores_trimmed(isDrops_trimmed == 1);
        birthEvents_durations = durations_trimmed(isDrops_trimmed == 1);
        
        
        % 8.  re-define period to begin at start of low nutrient pulse, by
        %     subtracting quarter period from corrected timestamp
        shifted_birthTimestamps = birthEvents_timestamps - (timescale/4);
        
        
        % 9. assign elements of timestamp vector to period fraction bins
        timeInPeriods = shifted_birthTimestamps/timescale; % unit = sec/sec
        timeInPeriodFraction = timeInPeriods - floor(timeInPeriods);
        assignedBin = ceil(timeInPeriodFraction * binsPerPeriod);
        
        
        % 10. bin nScores by period fraction
        % 10. (i) first remove all nScores with no value (NaN)
        %valued_nScores = birthEvents_nScores(birthEvents_nScores >= 0);
        %assignedBin_4valued_nScores = assignedBin(birthEvents_nScores >= 0);
        
        valued_nScores = birthEvents_nScores(birthEvents_durations >= 10);
        assignedBin_4valued_nScores = assignedBin(birthEvents_durations >= 10);
        
        % 10. (ii) bin
        nScores_binnedByPeriodFraction = accumarray(assignedBin_4valued_nScores, valued_nScores, [], @(x) {x});
        
        nScores_binned_mean = cellfun(@mean,nScores_binnedByPeriodFraction);
        nScores_binned_std = cellfun(@std,nScores_binnedByPeriodFraction);
        nScores_binned_count = cellfun(@length,nScores_binnedByPeriodFraction);
        nScores_binned_sem = nScores_binned_std./nScores_binned_count;
        
        
        % 11. plot mean and s.e.m. of cell cycle duration over nutrient period 
        % 11. (i) convert bin # to absolute time (sec)
        timePerBin = timescale/binsPerPeriod;  % in sec
        binPeriod = linspace(1, binsPerPeriod, binsPerPeriod);
        timePeriod = timePerBin*binPeriod';

        % 11. (ii) repeat quarter period on both sides and plot over period fraction
        quarterOne = linspace(1,binsPerPeriod/4,binsPerPeriod/4);
        quarterFour = linspace(binsPerPeriod*3/4+1,binsPerPeriod,binsPerPeriod/4);
        quarterZero = linspace((binsPerPeriod/4-1),0,binsPerPeriod/4)*-1;
        
        signal_quarterOne = nScores_binned_mean(quarterOne,1);
        signal_quarterFour = nScores_binned_mean(quarterFour,1);

        time_quarterOne = timePeriod(quarterOne,1) + timescale;
        time_quarterZero = quarterZero*timePerBin;
        
        mean_nScores_stitched = [signal_quarterFour; nScores_binned_mean; signal_quarterOne];
        time_stitched = [time_quarterZero'; timePeriod; time_quarterOne];
  
        figure(1)
        if condition == 3
            subplot(1,2,1)
            stem(time_stitched, mean_nScores_stitched,'Color',[0.25 0.25 0.9])
            hold on
            axis([min(time_stitched) max(time_stitched) 0 1])
            title(strcat('stable: ',date))
        elseif condition == 1
            subplot(1,2,2)
            stem(time_stitched, mean_nScores_stitched,'Color',[0 0.7 0.7])
            axis([min(time_stitched) max(time_stitched) 0 1])
            title(strcat('fluc: ',num2str(timescale)))
        end
        xlabel('time (sec)')
        ylabel('mean nScore (min)')
        grid on
       
        
        % 12. plot scatter of cell cycle duration over nutrient period 
        scatter_quarterOne = nScores_binnedByPeriodFraction(quarterOne,1);
        scatter_quarterFour = nScores_binnedByPeriodFraction(quarterFour,1);
        scatter_stitched = [scatter_quarterFour; nScores_binnedByPeriodFraction; scatter_quarterOne];
        
        sz = 20; % circle size
        figure(2)
        if condition == 3
            subplot(1,2,1)
            for b = 1:length(scatter_stitched)
                y = scatter_stitched{b};
                x = ones(length(y),1)*time_stitched(b);
                scatter(x,y,sz,[0.25 0.25 0.9]) % purple
                hold on
            end
            axis([min(time_stitched) max(time_stitched) 0 1])
            title(strcat('stable: ',date))
        else
            subplot(1,2,2)
            for b = 1:length(scatter_stitched)
                y = scatter_stitched{b};
                x = ones(length(y),1)*time_stitched(b);
                scatter(x,y,sz,[0 0.7 0.7]) % green
                hold on
            end
            axis([min(time_stitched) max(time_stitched) 0 1])
            title(strcat('fluc: ',num2str(timescale)))
        end
        xlabel('time (sec)')
        ylabel('nScore')
        grid on
        
        % 13. repeat analysis for fluctuating environment, plotting fluc data over stable
    end
    clear i environment 
    
    cd('/Users/jen/Documents/StockerLab/Data_analysis/horoscope')
    
    Fig1 = figure(1);
    saveas(Fig1,strcat('horoscope-nScore-mean-',num2str(timescale),'-',date),'epsc')
    close(Fig1)
    
    Fig2 = figure(2);
    saveas(Fig2,strcat('horoscope-nScore-scatter-',num2str(timescale),'-',date),'epsc')
    close(Fig2)
    
end


