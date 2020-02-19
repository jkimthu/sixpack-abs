%% distribute

%  Goal: plot distributions of... normalized by population average
%           1. cell cycle duration
%           2. cell volume at birth
%           3. added mass per cell cycle

%  Strategy:
%
%       0. initialize data & specify target concentration
%       1. create a directory of experiments with target concentration
%       2. for all experiments in target directory... accumulate cell size and curve duration data
%               3. move to experiment folder and build data matrix
%               4. for each condition with target concentration...
%                       5. build data maxtrix from data for current condition
%                       6. isolate volume(Va), curve duration, added volume per cell cycle, drop and time data
%                       7. isolate only data during which drop == 1 (birth event)
%                               - birth size = volume at birth event
%                               - one value for duration and added volume is gathered per cell cycle
%                                 (durations and added volume are assembled
%                                 in data matrix as final values, repeated
%                                 for all timepoints in a curve)
%                       8. remove data from cell cycles with added volume > 0
%                       9. remove data from cell cycles with durations shorter than 10 min (zero values are incomplete)
%                      10. trim data to stabilized / non-bubble timestamps
%                      11. trim outliers (those 3 std dev away from media) from final dataset
%                      12. calculate count number of data points per bin
%                      13. bin data and normalize bin counts by total counts
%                      14. plot pdf and histograms per experiment
%              15. repeat for each condition, storing fluc and stable data in their own variables
%      16.  after assembling data from all experiments, plot violin plots
%           comparing distrubtions between fluc and stable for each experiment
%


%  Last edit: jen, 2018 Feb 13
%  commit: remove data that is at least three standard deviations away from
%          median of distribution. also expanded longest div time allowed
%          to 90 mins (instead of 60)




% OK! Lez go!


%%
% 0. initialize data & specify target concentration

clear
clc
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')

dataIndex = find(~cellfun(@isempty,storedMetaData));
birthSizeData = cell(size(storedMetaData));
ccDurationsData = cell(size(storedMetaData));

% initialize summary vectors for calculated data
experimentCount = length(dataIndex);

% determine target concentration
targetConcentration = 0.0105; % average


% 1. create a directory of conditions with target concentration
targetConditions = cell(experimentCount,1);

for e = 1:experimentCount
    
    % identify conditions with target concentration
    index = dataIndex(e);
    concentrations = storedMetaData{index}.concentrations;
    
    % each cell represents an experiment, each value a condition of target concentration
    targetConditions{e} = find(concentrations == targetConcentration);
    
end
clear e index

%%
% 2. for all experiments in target directory... accumulate cell size and curve duration data
compiled_birthVol_fluc = [];
compiled_birthVol_stable = [];

compiled_durations_fluc = [];
compiled_durations_stable = [];

compiled_addedVol_fluc = [];
compiled_addedVol_stable = [];

for e = 1:experimentCount
    
    % exclude all experiments without specified nutrient concentration of interest
    if isempty(targetConditions{e})
        continue
    end
    
    % 3. move to experiment folder and load data
    
    % identify experiment by date
    index = dataIndex(e);
    date = storedMetaData{index}.date;
    
    % exclude outlier from analysis
    if strcmp(date, '2017-10-31') == 1 %|| strcmp (timescale, 'monod') == 1
        disp(strcat(date,': excluded from analysis'))
        continue
    end
    disp(strcat(date, ': analyze!'))
    
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    
    % load data
    timescale = storedMetaData{index}.timescale;
    if ischar(timescale) == 0
        filename = strcat('lb-fluc-',date,'-window5-width1p4-1p7-jiggle-0p5.mat');
    elseif strcmp(date,'2017-09-26') == 1
        filename = 'lb-monod-2017-09-26-window5-va-jiggle-c12-0p1-c3456-0p5-bigger1p8.mat';
    elseif strcmp(date, '2017-11-09') == 1
        filename = 'lb-control-2017-11-09-window5-width1p4-jiggle-0p5.mat';
    end
    load(filename,'D5','M','M_va','T')
    
    
    % 4. for each condition with target concentration...
    for i = 1:length(targetConditions{e})
        c = targetConditions{e}(i);
        
        % 5. build experiment data matrix
        display(strcat('Condition (',num2str(i),') of (',num2str(length(targetConditions{e})),'); experiment (', num2str(e),') of (', num2str(length(dataIndex)),')'))
        xy_start = storedMetaData{index}.xys(c,1);
        xy_end = storedMetaData{index}.xys(c,end);
        conditionData = buildDM(D5,M,M_va,T,xy_start,xy_end);
         
        % 6. isolate volume, duration, addedVolume, drop and time data
        durations = conditionData(:,8)/60;       % col 8 = calculated curve durations (sec) converted to min
        volumes = conditionData(:,14);        % col 14 = calculated va_vals (cubic um)
        addedVolume = conditionData(:,20);    % col 20 = added volume per cell cycle
        drops = conditionData(:,5);           % col 5 = 1 at birth, zero otherwise
        timestamps = conditionData(:,2)/3600; % time in seconds converted to hours
        clear conditionData

        % 7. isolate only data during which drop == 1 (birth event)
        birthVolumes = volumes(drops == 1);
        uniqueDurations = durations(drops == 1);
        birthTimes_all = timestamps(drops == 1);
        uniqueAddedVolumes = addedVolume(drops == 1);
        
        % 8. remove data from cell cycles with added volume > 0
        positiveVolumes = birthVolumes(uniqueAddedVolumes > 0);
        positiveDurations = uniqueDurations(uniqueAddedVolumes > 0);
        positiveAddedVolumes = uniqueAddedVolumes(uniqueAddedVolumes > 0);
        positive_birthTimes_fullCyclesOnly = birthTimes_all(uniqueAddedVolumes > 0);
        
        % 9. remove data from cell cycles with durations shorter than 10 min (zero values are incomplete)
        finalBirthVolumes = positiveVolumes(positiveDurations > 10);%
        finalDurations = positiveDurations(positiveDurations > 10);
        finalAddedVolumes = positiveAddedVolumes(positiveDurations > 10);
        birthTimes_fullCyclesOnly = positive_birthTimes_fullCyclesOnly(positiveDurations > 10);
        
        % 10. trim data to stabilized / non-bubble timestamps
        minTime = 3;  % hr
        maxTime = storedMetaData{index}.bubbletime(c);
        
        % trim data after minimum time
        % note: birth size data is trimmed using birth times
        %       duration(divTime) and addedVol data is trimmed using div_times
        %     - this is because any birth size within time window can count in dataset,
        %       whereas only full cycle cycles (bounded by drops) count for "per cell cycle" data
        birthTimes_fullCycle_trim1 = birthTimes_fullCyclesOnly(birthTimes_fullCyclesOnly >= minTime);
        finalBirthVolumes_trim1 = finalBirthVolumes(birthTimes_fullCyclesOnly >= minTime);%
        finalDurations_trim1 = finalDurations(birthTimes_fullCyclesOnly >= minTime);
        finalAddedVolumes_trim1 = finalAddedVolumes(birthTimes_fullCyclesOnly >= minTime);
        
        clear birthTimes birthVolumes divisionTimestamps finalDurations finalAddedVolumes
        
        % trim data after bubble appearance, if applicable (non-zero)
        if maxTime > 0
            birthVolumes_trim2 = finalBirthVolumes_trim1(birthTimes_fullCycle_trim1 <= maxTime);%
            durations_trim2 = finalDurations_trim1(birthTimes_fullCycle_trim1 <= maxTime);
            addedVolumes_trim2 = finalAddedVolumes_trim1(birthTimes_fullCycle_trim1 <= maxTime);
        else
            birthVolumes_trim2 = finalBirthVolumes_trim1;
            durations_trim2 = finalDurations_trim1;
            addedVolumes_trim2 = finalAddedVolumes_trim1;
        end
        clear birth_times_trim1 div_times_trim1 finalDurations_trim1 finalAddedVolumes_trim1
        

        %%
        % 11. trim outliers (those 3 std dev away from median) from final dataset
        medianBirthVol = median(birthVolumes_trim2);
        std_birthVol = std(birthVolumes_trim2);
        birthVolumes_minusUpper3std = birthVolumes_trim2(birthVolumes_trim2 <= (medianBirthVol+std_birthVol*3));
        trueBirthVolumes = birthVolumes_minusUpper3std(birthVolumes_minusUpper3std >= (medianBirthVol-std_birthVol*3));
        
        medianDuration = median(durations_trim2);
        std_duration = std(durations_trim2);
        durations_minusUpper3std = durations_trim2(durations_trim2 <= (medianDuration+std_duration*3));
        trueDurationTimes = durations_minusUpper3std(durations_minusUpper3std >= (medianDuration-std_duration*3));
        
        medianAddedVol = median(addedVolumes_trim2);
        std_addedVol = std(addedVolumes_trim2);
        addedVol_minusUpper3std = addedVolumes_trim2(addedVolumes_trim2 <= (medianAddedVol+std_addedVol*3));
        trueAddedVolumes = addedVol_minusUpper3std(addedVol_minusUpper3std >= (medianAddedVol-std_addedVol*3));
        
        % 11. calculate count number of data points per bin
        count_birthSize_true = length(trueBirthVolumes); % birth size
        count_duration = length(trueDurationTimes);      % curve duration
        count_addedVolume = length(trueAddedVolumes);    % added volume
        
        
        %%
        % 13. bin data and normalize bin counts by total counts
        % 14. plot pdf and histograms per experiment
        
        % birth size, for full cycles only
        binSize = 0.2; % cubic microns
        largestBirthSize = 10; % cubic microns
        
        binVector = linspace(0, largestBirthSize/binSize, largestBirthSize/binSize +1)' *0.2;
        assignedBins = ( ceil(trueBirthVolumes * 5) );
        
        binnedBirthVolumes_true = accumarray(assignedBins, trueBirthVolumes, [], @(x) {x});
        binCounts = cellfun(@length,binnedBirthVolumes_true);
        pdf_birthVolume = binCounts/count_birthSize_true;
        
        if length(pdf_birthVolume) < length(binVector)
            binVector = binVector(1:length(pdf_birthVolume));
        end
        
        fig_birthSize_true = figure(1);
        if c == 1
            birthVol_fluc{e} = trueBirthVolumes;
            compiled_birthVol_fluc = [compiled_birthVol_fluc; trueBirthVolumes];
            bar(binVector,pdf_birthVolume(1:length(binVector)),'FaceColor',[0 0.7 0.7])% green
            hold on
        else
            birthVol_stable{e} = trueBirthVolumes;
            compiled_birthVol_stable = [compiled_birthVol_stable; trueBirthVolumes];
            bar(binVector,pdf_birthVolume(1:length(binVector)),'FaceColor',[0.25 0.25 0.9])% purple
        end
        title(strcat(date,': n=',num2str(count_birthSize_true)))
        legend('fluc','stable')
        xlabel('cell volume at birth, full cycles only (cubic um)')
        ylabel('pdf')
        axis([0 largestBirthSize 0 .3])
        

        % cell cycle duration
        binSize = 2; % min
        longestDivTime = 90; % min
        
        binVector = linspace(0, longestDivTime/binSize, longestDivTime/binSize + 1)'*2;
        assignedBins = ( ceil(trueDurationTimes/2) );
        
        binnedDurations = accumarray(assignedBins, trueDurationTimes, [], @(x) {x});
        binCounts = cellfun(@length,binnedDurations);
        pdf_duration = binCounts/count_duration;
        
        if length(pdf_duration) < length(binVector)
            binVector = binVector(1:length(pdf_duration));
        end
        
        fig_duration = figure(2);
        if c == 1
            duration_fluc{e} = trueDurationTimes;
            compiled_durations_fluc = [compiled_durations_fluc; trueDurationTimes];
            bar(binVector,pdf_duration(1:length(binVector)),'FaceColor',[0 0.7 0.7])% green
            hold on
        else
            duration_stable{e} = trueDurationTimes;
            compiled_durations_stable = [compiled_durations_stable; trueDurationTimes];
            bar(binVector,pdf_duration(1:length(binVector)),'FaceColor',[0.25 0.25 0.9])% purple
        end
        title(strcat(date,': n=',num2str(count_duration)))
        legend('fluc','stable')
        xlabel('cell cycle duration (min)')
        ylabel('pdf')
        axis([0 longestDivTime 0 .3])
        
        
        % added volume per cell cycle
        binSize = 0.2; % cubic microns
        greatestAdd = 10;
        
        binVector = linspace(0, greatestAdd/binSize, greatestAdd/binSize +1)' *0.2;
        assignedBins = ( ceil(trueAddedVolumes * 5) );
        
        binnedAddedVol = accumarray(assignedBins, trueAddedVolumes, [], @(x) {x});
        binCounts = cellfun(@length,binnedAddedVol);
        pdf_addedVol = binCounts/count_addedVolume;
        if length(pdf_addedVol) < length(binVector)
            binVector = binVector(1:length(pdf_addedVol));
        end
        
        fig_addedVol = figure(3);
        if c == 1
            addedVol_fluc{e} = trueAddedVolumes;
            compiled_addedVol_fluc = [compiled_addedVol_fluc; trueAddedVolumes];
            bar(binVector,pdf_addedVol(1:length(binVector)),'FaceColor',[0 0.7 0.7])% green
            hold on
        else
            addedVol_stable{e} = trueAddedVolumes;
            compiled_addedVol_stable = [compiled_addedVol_stable; trueAddedVolumes];
            bar(binVector,pdf_addedVol(1:length(binVector)),'FaceColor',[0.25 0.25 0.9])% purple
        end
        title(strcat(date,': n=',num2str(count_addedVolume)))
        legend('fluc','stable')
        xlabel('added volume per cell cycle (cubic um)')
        ylabel('pdf')
        axis([0 greatestAdd 0 .3])
        
        %15. repeat for each condition, storing fluc and stable data in their own variables
    end
    
    saveas(fig_birthSize_true,strcat('pdf-birthSize-',date),'epsc')
    close(fig_birthSize_true)
    
    saveas(fig_duration,strcat('pdf-duration-',date),'epsc')
    close(fig_duration)
    
    saveas(fig_addedVol,strcat('pdf-addedVolume-',date),'epsc')
    close(fig_addedVol)
    
    clear fig_birthSize fig_duration fig_addedVol
    clear D5 M M_va T filename experimentFolder
    clear mean_birthSize std_birthSize count_birthSize sem_birthSize
    clear mean_duration std_duration count_duration sem_duration
    
end

cd('/Users/jen/Documents/StockerLab/Data_analysis/Monod compiled/LB')

% 16.  after assembling data from all experiments,
%      plot violin plots comparing distrubtions between fluc and stable for each experiment

% birth size
violin_birthVol = figure(20);
distributionPlot(birthVol_fluc,'widthDiv',[2 1],'histOri','left','color',[0 0.7 0.7],'showMM',2) % green
distributionPlot(gca,birthVol_stable,'widthDiv',[2 2],'histOri','right','color',[0.25 0.25 0.9],'showMM',2) % purple
ylim([0 largestBirthSize])
ylabel('volume at birth (cubic microns)')
saveas(violin_birthVol,strcat('violin-birthVol-upTil-',date),'epsc')

% cell cycle duraiton
violin_duration = figure(21);
distributionPlot(duration_fluc,'widthDiv',[2 1],'histOri','left','color',[0 0.7 0.7],'showMM',2) % green
distributionPlot(gca,duration_stable,'widthDiv',[2 2],'histOri','right','color',[0.25 0.25 0.9],'showMM',2) % purple
ylim([0 longestDivTime])
ylabel('cell cycle duration (min)')
saveas(violin_duration,strcat('violin-duration-upTil-',date),'epsc')

% added volume
violin_addedVol = figure(22);
distributionPlot(addedVol_fluc,'widthDiv',[2 1],'histOri','left','color',[0 0.7 0.7],'showMM',2) % green
distributionPlot(gca,addedVol_stable,'widthDiv',[2 2],'histOri','right','color',[0.25 0.25 0.9],'showMM',2) % purple
ylim([0 greatestAdd])
ylabel('volume added per cell cycle (cubic microns)')
saveas(violin_addedVol,strcat('violin-addedVolume-upTil-',date),'epsc')



