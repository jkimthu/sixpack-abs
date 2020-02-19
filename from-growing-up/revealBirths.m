%% revealBirths


% Goal: are births synchronized? plot multiple looks at birth over time, normalized by current cell count.
%
%           1. absolute track count over time
%           2. births normalized by current track count over time - accounts for changing population size
%           3. births normalized by cell count over period fraction - looks for environmentally induced trends
%   


%  Last edit: Jen Nguyen, 2018 Feb 23



%  Strategy:
%
%     0. initialize complete meta data
%     0. for experiments of interest...
%           1. collect experiment date and exclude outliers (2017-10-31 and monod experiments)
%           2. initialize binning parameters
%           3. load measured data for stable condition
%           4. for stable average... 
%                  5. isolate isDrop and timestamp data
%                  6. correct time based on calculated lag in signal between junc and xy position
%                  7. isolate birth events (isDrop == 1) and corresponding timestamps
%                  8. PLOT ONE: number of tracked cells over time
%                  9. PLOT TWO: birth events over time, normalized by tracks per bin
%                 10. PLOT THREE: births normalized by cell count over period fraction
%          11. repeat analysis for fluctuating environment, plotting fluc data over stable



%           5. find average growth rate of stable average condition
%                       i. isolate data
%                      ii. remove data not in stabilized region
%                     iii. remove zeros from mu data (always bounding start and end of tracks)
%                      iv. calculate mean value for of mu and bvpr in stable
%           6. for fluctuating condition, load measured data
%                       i. isolate data of interest
%                      ii. normalize mu and bvpr data by mean of stable average condition
%                     iii. remove data not in stabilized region
%                      iv. remove zeros from mu data (always bounding start and end of tracks)
%           7. accumulate data by shifted time bin (period fraction)
%                       i. from original timestamp, subtract shift = period/4 + offset
%                      ii. re-define period to begin at start of high nutrient pulse
%                     iii. bin data by period fraction
%           8.  convert bin # to absolute time (in seconds)
%           9.  calculate average and s.e.m. per timebin
%          10.  plot, with repeated high nutrient half period at end
%    11. repeat for all fluctuating experiments


% OK! Lez go!


%%
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
    binsPerPeriod_resolved = 12;
    
    if timescale == 30
        
        binsPerPeriod_smooth = 0.05; % for smooth curves
        
    elseif timescale == 300
        
        binsPerPeriod_smooth = 0.5;
        
    elseif timescale == 900
        
        binsPerPeriod_smooth = 1.5;
        
    elseif timescale == 3600
        
        binsPerPeriod_smooth = 6;
        
    end
    
    
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
        

        % 5. isolate birth event, track number, and timestamp data
        isDrops = conditionData(:,5);           % col 5  =  isDrop; 0 during curve, 1 at birth event
        Time = conditionData(:,2);              % col 2  =  timestamps in sec
        stagePositions = conditionData(:,24);   % col 24 =  stage number, xy position from which data originates
        trackNum = conditionData(:,27);         % col 27 =  track count (not ID from particle tracking)
        
        
        % 6. correct time based on calculated lag in signal between junc and xy position
        if condition == 1 % IN FLUCTUATING ENVIRONMENT

            % 6i. initialize xy positions, lag times, and empty vector for compiling corrected timestamps
            XYs = unique(stagePositions);
            [lagTimes,~] = calculateLag(e);
            
            % 6ii. for each xy position
            correctedTimes = [];
            for position = 1:length(XYs)
                
                % iii. identify position and corresponding lag time
                currentXY = XYs(position);
                currentLag = lagTimes(position);
                
                % iv. subtract lag time from timestamp, to re-align cell experience (xy) with generated signal (junc)
                edits = Time(stagePositions==currentXY) - currentLag;
                
                % v. re-assign
                correctedTimes = [correctedTimes; edits];
                
            end
            clear currentXY currentLag edits XYs lagTimes position
                
        else % in stable environemnt, use original timestamps
            
            correctedTimes = Time;
            
        end
        
        
        % 7. isolate birth events (isDrop == 1) and corresponding timestamps
        birthEvents = isDrops(isDrops == 1);
        birthEvent_timestamps = correctedTimes(isDrops == 1);
        
        
        % 8. PLOT ONE: total tracks over time
        
        % 8i. accumulate track numbers by time bin
        periodsPerHour = 3600/timescale;
        binsPerHour = binsPerPeriod_smooth * periodsPerHour;
        timeBins_all = ceil(correctedTimes/3600*binsPerHour);
        binnedTracks = accumarray(timeBins_all,trackNum,[],@(x) {x});
        
        % 8ii. find unique track numbers in each time bin
        binnedTracks_unique = cellfun(@unique,binnedTracks,'UniformOutput',0);
        trackCountsPerBin = cellfun(@length,binnedTracks_unique);
        
        % 8iii. convert bin # to absolute time
        binVector = linspace(1, 10*binsPerHour, 10*binsPerHour);
        timeVector = binVector'./binsPerHour;
        
        % 8iv. extend births vector to full experiment time and plot
        padding = zeros(length(timeVector)-max(timeBins_all),1);

        figure(1)
        if condition == 3
            stem(timeVector,[trackCountsPerBin; padding],'Color',[0.25 0.25 0.9]) % dark purple
            hold on
        elseif condition == 1
            stem(timeVector,[trackCountsPerBin; padding],'Color',[0 0.7 0.7]) % teal green
        end
        axis([0 10.2 0 ceil(mean(trackCountsPerBin))*1.7])
        title(strcat('tracked cells over time : ',num2str(date),', timescale:',num2str(timescale)))
        xlabel('time (hr)')
        ylabel('unique tracks')
        legend('stable','fluc')
        
        
        % 9. PLOT TWO: births normalized by current track count over time

        % 9i. accumulate birth events by time bin
        timeBins_births = ceil(birthEvent_timestamps/3600*binsPerHour);
        binnedDrops = accumarray(timeBins_births,birthEvents,[],@(x) {x});
        
        % 9ii. count birth events per timebin
        birthsPerBin = cellfun(@sum,binnedDrops);
        
        % 9iii. normalize birth event count by track count, for each time bin
        birthsPerBin_normalized = birthsPerBin./trackCountsPerBin;
        
        % 9iv. extend births vector to full experiment time and plot
        padding = zeros(length(timeVector)-max(timeBins_births),1);
        
        figure(2)
        if condition == 3
            stem(timeVector,[birthsPerBin_normalized; padding],'Color',[0.25 0.25 0.9]) % dark purple
            hold on
        elseif condition == 1
            stem(timeVector,[birthsPerBin_normalized; padding],'Color',[0 0.7 0.7]) % teal green
        end
        axis([0 10.2 0 0.5])
        title(strcat('normalized birth events over time : ',num2str(date),', timescale:',num2str(timescale)))
        xlabel('time (hr)')
        ylabel('reproducing population fraction')
        legend('stable','fluc')
        

        % 10. PLOT THREE: births normalized by cell count over period fraction
        
        % 10i. re-define period to begin at start of low nutrient pulse, by
        %      subtracting quarter period from corrected timestamp
        shifted_birthTimestamps = birthEvent_timestamps - (timescale/4);
        
        % 10ii. trim data to timepoints between 3 hrs and bubbles
        shifted_birthTimestamps_trim1 = shifted_birthTimestamps(shifted_birthTimestamps/3600 >= 3);
        birthEvents_trim1 = birthEvents(shifted_birthTimestamps/3600 >= 3);
        
        if bubbleTime(condition) > 0
            
            shifted_birthTimestamps_trim2 = shifted_birthTimestamps_trim1(shifted_birthTimestamps_trim1/3600 <= bubbleTime(condition));
            birthEvents_trim2 = birthEvents_trim1(shifted_birthTimestamps_trim1/3600 <= bubbleTime(condition));
            
        else
            
            shifted_birthTimestamps_trim2 = shifted_birthTimestamps_trim1;
            birthEvents_trim2 = birthEvents_trim1;
            
        end
        
        % 10iii. bin birth data by period fraction
        timeInPeriods_births = shifted_birthTimestamps_trim2/timescale; % unit = sec/sec
        timeInPeriodFraction_births = timeInPeriods_births - floor(timeInPeriods_births);
        assignedBin_birthTimestamps = ceil(timeInPeriodFraction_births * binsPerPeriod_resolved);
        
        births_binnedByPeriodFraction = accumarray(assignedBin_birthTimestamps, birthEvents_trim2, [], @(x) {x});
        birthsPerPeriodFraction = cellfun(@sum,births_binnedByPeriodFraction);
        
        % 10iv.  convert bin # to absolute time (sec)
        timePerBin = timescale/binsPerPeriod_resolved;  % in sec
        binPeriod = linspace(1, binsPerPeriod_resolved, binsPerPeriod_resolved);
        timePeriod = timePerBin*binPeriod';
         
        % 10v. normalize binned births by total births
        birthsPerPeriodFraction_normalized = birthsPerPeriodFraction./sum(birthEvents_trim2);
        
        
        % 10 vi. repeat quarter period on both sides and plot over period fraction
        quarterOne = linspace(1,binsPerPeriod_resolved/4,binsPerPeriod_resolved/4);
        quarterFour = linspace(binsPerPeriod_resolved*3/4+1,binsPerPeriod_resolved,binsPerPeriod_resolved/4);
        quarterZero = linspace((binsPerPeriod_resolved/4-1),0,binsPerPeriod_resolved/4)*-1;
        
        signal_quarterOne = birthsPerPeriodFraction_normalized(quarterOne,1);
        signal_quarterFour = birthsPerPeriodFraction_normalized(quarterFour,1);
        
        time_quarterOne = timePeriod(quarterOne,1) + timescale;
        time_quarterZero = quarterZero*timePerBin;
        
        birthSignal_stitched = [signal_quarterFour; birthsPerPeriodFraction_normalized; signal_quarterOne];
        time_stitched = [time_quarterZero'; timePeriod; time_quarterOne];
  
        figure(3)
        if condition == 3
            subplot(1,2,1)
            stem(time_stitched, birthSignal_stitched,'Color',[0.25 0.25 0.9])
            hold on
            axis([min(time_stitched) max(time_stitched) 0 0.2])
            title(strcat('stable: ',date))
        elseif condition == 1
            subplot(1,2,2)
            stem(time_stitched, birthSignal_stitched,'Color',[0 0.7 0.7])
            axis([min(time_stitched) max(time_stitched) 0 0.2])
            title(strcat('fluc: ',num2str(timescale)))
        end
        xlabel('time (sec)')
        ylabel('birth events per period fraction, normalized by total')
        grid on
        
        %clear birthEvents birthEvent_timestamps binVector timeVector padding
        %clear binnedDrops isDrops stagePositions Time condition xy_start xy_end 
        
    end
    clear i environment 
    
    cd('/Users/jen/Documents/StockerLab/Data_analysis/birthEvents_v_time')
    
    Fig1 = figure(1);
    saveas(Fig1,strcat('trackCountvTime-',num2str(timescale),'-',date),'epsc')
    close(Fig1)
    
    Fig2 = figure(2);
    saveas(Fig2,strcat('normalizedBirthEventsvTime-',num2str(timescale),'-',date),'epsc')
    close(Fig2)
    
    Fig3 = figure(3);
    saveas(Fig3,strcat('normalizedBirthEventsvPeriodFrac-',num2str(timescale),'-',date),'epsc')
    close(Fig3)

end





