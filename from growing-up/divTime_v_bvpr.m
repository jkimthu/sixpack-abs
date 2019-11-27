% divisionTime_v_bvpr

% goal: plot division time (cell cycle duration) against biovolume production rate,
%       for a specific concentration

% strategy:
%
%       0. initialize data
%       0. determine target concentration
%       1. create a directory of experiments with target concentration
%       2. for all experiments in target directory... accumulate cell size stats
%               3. move to experiment folder and build data matrix
%               4. for each condition with target concentration...
%                       5. isolate data from current condition
%                       6. isolate curve durations, drop and time data
%                       7. isolate only data during which drop == 1 (birth event)
%                       8. trim data to stabilized / non-bubble timestamps
%                       9. calculate mean and s.e.m. of size
%                      10. accumulate data for storage and plotting
%              11. store data from all conditions into measured data structure
%      12.  plot average cell cycle duration against corresponding biovol production rate
%


% last updated: 2018 January 19

% OK let's go!!

%% 0. initialize data

clear
clc

cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')

dataIndex = find(~cellfun(@isempty,storedMetaData));
ccDurationsData = cell(size(storedMetaData));

% initialize summary vectors for calculated data
experimentCount = length(dataIndex);

% determine target concentration
targetConcentration = 0.0105; % average

%% 1. create a directory of experiments with target concentration
targetConditions = cell(1,experimentCount);

for e = 1:experimentCount
    
    % identify conditions with target concentration
    index = dataIndex(e);
    concentrations = storedMetaData{index}.concentrations;
    
    targetConditions{e} = find(concentrations == targetConcentration);
    
end
clear e index

%% 2. for all experiments in target directory... calculate cell size stats

for e = 1:experimentCount
    
    if isempty(targetConditions{e})
        continue
    end

    % 3. move to experiment folder and build data matrix
    
    % identify experiment by date
    index = dataIndex(e);
    date = storedMetaData{index}.date;
    
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
    
    % build experiment data matrix
    display(strcat('Experiment (', num2str(e),') of (', num2str(length(dataIndex)),')'))
    exptData = buildDM(D5,M,M_va,T);
    
    clear D5 M M_va T filename experimentFolder
    
    
    % 4. for each condition with target concentration...
    for i = 1:length(targetConditions{e})
        c = targetConditions{e}(i);
        
        % 5. isolate data from current condition
        conditionData = exptData(exptData(:,35) == c,:);
        
        % 6. isolate volume, drop and time data
        durations = conditionData(:,8);       % col 8 = calculated curve durations (sec)
        drops = conditionData(:,5);           % col 5 = 1 at birth, zero otherwise
        timestamps = conditionData(:,2)/3600; % time in seconds converted to hours
        clear conditionData
        
        % 7. isolate only data during which drop == 1 (birth event)
        uniqueDurations = durations(drops == 1);
        birthTimes = timestamps(drops == 1);
        
        % 8. remove zeros from duration data (these are for births w/o associated division)
        finalDurations = uniqueDurations(uniqueDurations > 0);
        divisionTimestamps = birthTimes(uniqueDurations > 0);
        
        % 8. trim data to stabilized / non-bubble timestamps
        minTime = 3;  % hr
        maxTime = storedMetaData{index}.bubbletime(c);
        
        times_trim1 = divisionTimestamps(divisionTimestamps >= minTime);
        finalDurations_trim1 = finalDurations(divisionTimestamps >= minTime);
        clear divisionTimestamps finalDurations
        
        if maxTime > 0
            trueDurations_trim2 = finalDurations_trim1(times_trim1 <= maxTime);
        else
            trueDurations_trim2 = finalDurations_trim1;
        end
        clear times_trim1
        
        % 9. calculate mean and s.e.m. of size
        mean_duration = mean(trueDurations_trim2);
        count_duration = length(trueDurations_trim2);
        std_duration = std(trueDurations_trim2);
        sem_duration = std_duration./sqrt(count_duration);
        
        % 10. accumulate data for storage / plotting
        compiledDuration{c}.date = date;
        compiledDuration{c}.timescale = timescale;
        compiledDuration{c}.condition = c;
        compiledDuration{c}.mean = mean_duration;
        compiledDuration{c}.std = std_duration;
        compiledDuration{c}.count = count_duration;
        compiledDuration{c}.sem = sem_duration;
        
    end
    
    clear mean_duration std_duration count_duration sem_duration
    
    % 11. store data from all conditions into measured data structure
    ccDurationsData{index} = compiledDuration;
    clear compiledDuration

end
%% 12.  plot average cell size at birth against corresponding biovol production rate

cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('bioProdRateData.mat')
load('measuredData.mat')

%%
% initialize data vector for easy plotting
counter = 0;
summaryDurations = zeros(1,(experimentCount-1)*2);
summaryStds_durations = zeros(1,(experimentCount-1)*2);
summarySems_durations = zeros(1,(experimentCount-1)*2);
summaryBVPRs_durations = zeros(1,(experimentCount-1)*2);
summaryMus_durations = zeros(1,(experimentCount-1)*2);
summaryTimescales = cell(1,(experimentCount-1)*2);


for e = 1:experimentCount
    
    % identify conditions with target concentration
    index = dataIndex(e);
    
    for i = 1:length(targetConditions{e})
        c = targetConditions{e}(i);
        
        counter = counter + 1;
        summaryDurations(counter) = ccDurationsData{index}{c}.mean/60;
        summaryStds_durations(counter) = ccDurationsData{index}{c}.std;
        summarySems_durations(counter) = ccDurationsData{index}{c}.sem;
        summaryBVPRs_durations(counter) = bioProdRateData{index}{c}.mean;
        summaryMus_durations(counter) = measuredData{index}.individuals{c}.muMean;
        
        summaryTimescales{counter} = ccDurationsData{index}{c}.timescale;
    end
end



figure(1)
for p = 1:counter
    
    if mod(p,2) == 0
        h(p) = errorbar(summaryBVPRs_durations(p),summaryDurations(p),summarySems_durations(p),'o','Color',[0.25 0.25 0.9],'MarkerSize',10);
        hold on
   
    elseif summaryTimescales{p} == 30
        h(p) = errorbar(summaryBVPRs_durations(p),summaryDurations(p),summarySems_durations(p),'o','Color',[1 0 0],'MarkerSize',10); % red
        hold on

    elseif summaryTimescales{p} == 300
        h(p) = errorbar(summaryBVPRs_durations(p),summaryDurations(p),summarySems_durations(p),'o','Color',[1 0.85 0.01],'MarkerSize',10); % sunflower yellow
        hold on
        
    elseif summaryTimescales{p} == 900
        h(p) = errorbar(summaryBVPRs_durations(p),summaryDurations(p),summarySems_durations(p),'o','Color',[0 0.7 0.7],'MarkerSize',10); % green
        hold on
        
    end
    
end
xlabel('biovol prod rate (cubic um/hr)')
ylabel('cell cycle duration (min)')
legend('30 sec','5 min','15 min','stable')
axis([2 10 -20 70])

figure(2)
for p = 1:counter
    
    if mod(p,2) == 0
        h(p) = errorbar(summaryMus_durations(p),summaryDurations(p),summarySems_durations(p),'o','Color',[0.25 0.25 0.9],'MarkerSize',10);
        axis([0.5 3 -20 70])
        hold on
   
    elseif summaryTimescales{p} == 30
        h(p) = errorbar(summaryMus_durations(p),summaryDurations(p),summarySems_durations(p),'o','Color',[1 0 0],'MarkerSize',10); % red
        hold on

    elseif summaryTimescales{p} == 300
        h(p) = errorbar(summaryMus_durations(p),summaryDurations(p),summarySems_durations(p),'o','Color',[1 0.85 0.01],'MarkerSize',10); % sunflower yellow
        hold on
        
    elseif summaryTimescales{p} == 900
        h(p) = errorbar(summaryMus_durations(p),summaryDurations(p),summarySems_durations(p),'o','Color',[0 0.7 0.7],'MarkerSize',10); % green
        hold on
        
    end
    
end
xlabel('doubling rate of volume (1/hr)')
ylabel('cell cycle duration (min)')
legend('30 sec','5 min','15 min','stable')

