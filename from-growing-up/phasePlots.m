% phasePlots

% goal: using data from a specified concentration (average), generate four phase plots
%           1. cell volume at birth vs. biovolume production rate
%           2. cell volume at birth vs. mu
%           3. cell cycle duration vs. biovolume production rate
%           4. cell cycle duration vs. mu
%           5. mu vs bvpr
%           6. bvpr vs. mu


% strategy:
%
%       0. initialize data & specify target concentration
%       1. create a directory of experiments with target concentration
%       2. for all experiments in target directory... accumulate cell size and curve duration data
%               3. move to experiment folder and build data matrix
%               4. for each condition with target concentration...
%                       5. isolate data from current condition
%                       6. isolate size, curve duration, drop and time data
%                       7. isolate only data during which drop == 1 (gives size at birth event, and one curve duration per curve)
%                               note: separate size and duration time vectors, as lengths may differ
%                       8. trim data to stabilized / non-bubble timestamps
%                       9. calculate mean and s.e.m. of size, curve duration
%                      10. accumulate data for storage and plotting
%              11. store data from all conditions into measured data structure
%      12.  plot average and s.e.m. against corresponding biovol production rate
%


% last updated: 2018 Mar 12
% commit: add figures (5) and (6) which plot mu and bvpr against each other

% OK let's go!!

%% 0. initialize data & specify target concentration

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
targetConditions = cell(1,experimentCount);

for e = 1:experimentCount
    
    % identify conditions with target concentration
    index = dataIndex(e);
    concentrations = storedMetaData{index}.concentrations;
    
    % each cell represents an experiment, each value a condition of target concentration
    targetConditions{e} = find(concentrations == targetConcentration);
    
end
clear e index

%% 2. accumulate cell size and cell cycle duration data

for e = 1:experimentCount
    
    % exclude all experiments without specified nutrient concentration of interest
    if isempty(targetConditions{e})
        continue
    end

    % 3. move to experiment folder and build data matrix
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
        
        % build experiment data matrix
        display(strcat('Experiment (', num2str(e),') of (', num2str(length(dataIndex)),'), condition: ', num2str(c)))
        xy_start = storedMetaData{index}.xys(c,1);
        xy_end = storedMetaData{index}.xys(c,end);
        exptData = buildDM(D5, M, M_va, T, xy_start, xy_end, e);
        
        % 5. isolate data from current condition
        conditionData = exptData(exptData(:,28) == c,:);
        
        % 6. isolate volume, duration, drop and time data
        durations = conditionData(:,8);       % col 8 = calculated curve durations (sec)
        volumes = conditionData(:,14);        % col 14 = calculated va_vals (cubic um)
        drops = conditionData(:,5);           % col 5 = 1 at birth, zero otherwise
        timestamps = conditionData(:,2)/3600; % time in seconds converted to hours
        clear conditionData
        
        % 7. isolate only data during which drop == 1 (birth event)
        birthVolumes = volumes(drops == 1);
        uniqueDurations = durations(drops == 1);
        birthTimes = timestamps(drops == 1);
        
        % 8. remove zeros from duration data (these are for births w/o associated division)
        finalDurations = uniqueDurations(uniqueDurations > 0);
        divisionTimestamps = birthTimes(uniqueDurations > 0);
        
        % 9. trim data to stabilized / non-bubble timestamps
        minTime = 3;  % hr
        maxTime = storedMetaData{index}.bubbletime(c);
        
        birth_times_trim1 = birthTimes(birthTimes >= minTime);
        div_times_trim1 = divisionTimestamps(divisionTimestamps >= minTime);
        
        volumes_trim1 = birthVolumes(birthTimes >= minTime);
        finalDurations_trim1 = finalDurations(divisionTimestamps >= minTime);
        
        clear birthTimes birthVolumes divisionTimestamps finalDurations
        
        if maxTime > 0
            volumes_trim2 = volumes_trim1(birth_times_trim1 <= maxTime);
            durations_trim2 = finalDurations_trim1(div_times_trim1 <= maxTime);
        else
            volumes_trim2 = volumes_trim1;
            durations_trim2 = finalDurations_trim1;
        end
        clear birth_times_trim1 div_times_trim1
        
        
        % 10. trim outliers (those 3 std dev away from sampled median) from final dataset
        sigma = 3;
        
        medianBirthVol = median(volumes_trim2);
        std_birthVol = std(volumes_trim2);
        birthVolumes_minusUpper3std = volumes_trim2(volumes_trim2 <= (medianBirthVol+std_birthVol*sigma));
        trueBirthVolumes = birthVolumes_minusUpper3std(birthVolumes_minusUpper3std >= (medianBirthVol-std_birthVol*sigma));
        
        medianDuration = median(durations_trim2);
        std_duration = std(durations_trim2);
        durations_minusUpper3std = durations_trim2(durations_trim2 <= (medianDuration+std_duration*sigma));
        trueDurationTimes = durations_minusUpper3std(durations_minusUpper3std >= (medianDuration-std_duration*sigma));
        
        
        % 10. calculate mean and s.e.m. of size
        mean_birthSize = mean(trueBirthVolumes);
        count_birthSize = length(trueBirthVolumes);
        std_birthSize = std(trueBirthVolumes);
        sem_birthSize = std_birthSize./sqrt(count_birthSize);
        
        % calculate mean and s.e.m. of curve duration
        mean_duration = mean(trueDurationTimes);
        count_duration = length(trueDurationTimes);
        std_duration = std(trueDurationTimes);
        sem_duration = std_duration./sqrt(count_duration);
        
       
        % 10. accumulate data for storage / plotting
        compiledBirthSize{c}.date = date;
        compiledBirthSize{c}.timescale = timescale;
        compiledBirthSize{c}.condition = c;
        compiledBirthSize{c}.mean = mean_birthSize;
        compiledBirthSize{c}.std = std_birthSize;
        compiledBirthSize{c}.count = count_birthSize;
        compiledBirthSize{c}.sem = sem_birthSize;
        
        compiledDuration{c}.date = date;
        compiledDuration{c}.timescale = timescale;
        compiledDuration{c}.condition = c;
        compiledDuration{c}.mean = mean_duration;
        compiledDuration{c}.std = std_duration;
        compiledDuration{c}.count = count_duration;
        compiledDuration{c}.sem = sem_duration;
        
    end
    clear D5 M M_va T filename experimentFolder
    clear mean_birthSize std_birthSize count_birthSize sem_birthSize
    clear mean_duration std_duration count_duration sem_duration
    
    % 11. store data from all conditions into measured data structure
    birthSizeData{index} = compiledBirthSize;
    ccDurationsData{index} = compiledDuration;
    clear compiledBirthSize compiledDuration

end
%% 12.  plot !!

cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('bioProdRateData.mat')

% measuredData is a variable created in plotMonod.m
load('measuredData.mat')

% initialize data vector for easy plotting
counter = 0;
summaryDates = cell(1,(experimentCount-1)*2);
summaryBVPRs = zeros(1,(experimentCount-1)*2);
summaryMus = zeros(1,(experimentCount-1)*2);
summaryTimescales = cell(1,(experimentCount-1)*2);

summarySizes = zeros(1,(experimentCount-1)*2);
summaryStds_size = zeros(1,(experimentCount-1)*2);
summarySems_size = zeros(1,(experimentCount-1)*2);

summaryDurations = zeros(1,(experimentCount-1)*2);
summaryStds_durations = zeros(1,(experimentCount-1)*2);
summarySems_durations = zeros(1,(experimentCount-1)*2);


for e = 1:experimentCount
    
    % identify conditions with target concentration
    index = dataIndex(e);
    date = storedMetaData{index}.date;
    
    % exclude outlier from analysis
    if strcmp(date, '2017-10-31') == 1 %|| strcmp (timescale, 'monod') == 1
        disp(strcat(date,': excluded from analysis'))
        continue
    end
    disp(strcat(date, ': analyze!'))
    
    for i = 1:length(targetConditions{e})
        c = targetConditions{e}(i);
        counter = counter + 1;
        
        % growth rate and meta data
        summaryDates{counter} = date;
        summaryBVPRs(counter) = bioProdRateData{index}{c}.mean;
        summaryMus(counter) = measuredData{index}.individuals{c}.muMean;
        summaryTimescales{counter} = birthSizeData{index}{c}.timescale;
        
        % cell volume
        summarySizes(counter) = birthSizeData{index}{c}.mean;
        summaryStds_size(counter) = birthSizeData{index}{c}.std;
        summarySems_size(counter) = birthSizeData{index}{c}.sem;  
        
        % cell cycle duration
        summaryDurations(counter) = ccDurationsData{index}{c}.mean/60;
        summaryStds_durations(counter) = ccDurationsData{index}{c}.std/60;
        summarySems_durations(counter) = ccDurationsData{index}{c}.sem/60;
        
    end
end


shift = 0.01;

% figure 1. cell volume at birth vs. biovolume production rate
figure(1) 
for p = 1:counter
    
    if mod(p,2) == 0
        h(p) = errorbar(summaryBVPRs(p),summarySizes(p),summarySems_size(p),'o','Color',[0.25 0.25 0.9],'MarkerSize',10);
        text(summaryBVPRs(p)+shift,summarySizes(p)+shift, summaryDates(p),'FontSize',10);
        hold on
   
    elseif summaryTimescales{p} == 30
        h(p) = errorbar(summaryBVPRs(p),summarySizes(p),summarySems_size(p),'o','Color',[1 0 0],'MarkerSize',10); % red
        text(summaryBVPRs(p)+shift,summarySizes(p)+shift, summaryDates(p),'FontSize',10);
        hold on

    elseif summaryTimescales{p} == 300
        h(p) = errorbar(summaryBVPRs(p),summarySizes(p),summarySems_size(p),'o','Color',[1 0.85 0.01],'MarkerSize',10); % sunflower yellow
        text(summaryBVPRs(p)+shift,summarySizes(p)+shift, summaryDates(p),'FontSize',10);
        hold on
        
    elseif summaryTimescales{p} == 900
        h(p) = errorbar(summaryBVPRs(p),summarySizes(p),summarySems_size(p),'o','Color',[0 0.7 0.7],'MarkerSize',10); % green
        text(summaryBVPRs(p)+shift,summarySizes(p)+shift, summaryDates(p),'FontSize',10);
        hold on
        
    elseif summaryTimescales{p} == 3600
        h(p) = errorbar(summaryBVPRs(p),summarySizes(p),summarySems_size(p),'o','Color',[1 0.5 0.5],'MarkerSize',10); % peach
        text(summaryBVPRs(p)+shift,summarySizes(p)+shift, summaryDates(p),'FontSize',10);
        hold on
        
    end
    
end
xlabel('biovol prod rate (cubic um/hr)')
ylabel('cell volume at birth (cubic um)')
legend('30 sec','5 min','15 min','60min','stable')
axis([2 10 1 5])


% figure 2. cell volume at birth vs. mu
figure(2)
for p = 1:counter
    
    if mod(p,2) == 0
        h(p) = errorbar(summaryMus(p),summarySizes(p),summarySems_size(p),'o','Color',[0.25 0.25 0.9],'MarkerSize',10);
        text(summaryMus(p)+shift,summarySizes(p)+shift, summaryDates(p),'FontSize',10);
        hold on
   
    elseif summaryTimescales{p} == 30
        h(p) = errorbar(summaryMus(p),summarySizes(p),summarySems_size(p),'o','Color',[1 0 0],'MarkerSize',10); % red
        text(summaryMus(p)+shift,summarySizes(p)+shift, summaryDates(p),'FontSize',10);
        hold on

    elseif summaryTimescales{p} == 300
        h(p) = errorbar(summaryMus(p),summarySizes(p),summarySems_size(p),'o','Color',[1 0.85 0.01],'MarkerSize',10); % sunflower yellow
        text(summaryMus(p)+shift,summarySizes(p)+shift, summaryDates(p),'FontSize',10);
        hold on
        
    elseif summaryTimescales{p} == 900
        h(p) = errorbar(summaryMus(p),summarySizes(p),summarySems_size(p),'o','Color',[0 0.7 0.7],'MarkerSize',10); % green
        text(summaryMus(p)+shift,summarySizes(p)+shift, summaryDates(p),'FontSize',10);
        hold on
        
    elseif summaryTimescales{p} == 3600
        h(p) = errorbar(summaryMus(p),summarySizes(p),summarySems_size(p),'o','Color',[1 0.5 0.5],'MarkerSize',10); % peach
        text(summaryMus(p)+shift,summarySizes(p)+shift, summaryDates(p),'FontSize',10);
        hold on
        
    end
    
end
xlabel('doubling rate of volume (1/hr)')
ylabel('cell volume at birth (cubic um)')
legend('30 sec','5 min','15 min','60 min','stable')
axis([0.5 3 1 5])


% figure 3. cell cycle duration vs. biovolume production rate
figure(3) 
for p = 1:counter
    
    if mod(p,2) == 0
        h(p) = errorbar(summaryBVPRs(p),summaryDurations(p),summarySems_durations(p),'o','Color',[0.25 0.25 0.9],'MarkerSize',10);
        text(summaryBVPRs(p)+shift,summaryDurations(p)+shift, summaryDates(p),'FontSize',10);
        hold on
   
    elseif summaryTimescales{p} == 30
        h(p) = errorbar(summaryBVPRs(p),summaryDurations(p),summarySems_durations(p),'o','Color',[1 0 0],'MarkerSize',10); % red
        text(summaryBVPRs(p)+shift,summaryDurations(p)+shift, summaryDates(p),'FontSize',10);
        hold on

    elseif summaryTimescales{p} == 300
        h(p) = errorbar(summaryBVPRs(p),summaryDurations(p),summarySems_durations(p),'o','Color',[1 0.85 0.01],'MarkerSize',10); % sunflower yellow
        text(summaryBVPRs(p)+shift,summaryDurations(p)+shift, summaryDates(p),'FontSize',10);
        hold on
        
    elseif summaryTimescales{p} == 900
        h(p) = errorbar(summaryBVPRs(p),summaryDurations(p),summarySems_durations(p),'o','Color',[0 0.7 0.7],'MarkerSize',10); % green
        text(summaryBVPRs(p)+shift,summaryDurations(p)+shift, summaryDates(p),'FontSize',10);
        hold on
        
    elseif summaryTimescales{p} == 3600
        h(p) = errorbar(summaryBVPRs(p),summaryDurations(p),summarySems_durations(p),'o','Color',[1 0.5 0.5],'MarkerSize',10); % peach
        text(summaryBVPRs(p)+shift,summaryDurations(p)+shift, summaryDates(p),'FontSize',10);
        hold on
        
    end
    
end
xlabel('biovol prod rate (cubic um/hr)')
ylabel('cell cycle duration (min)')
legend('30 sec','5 min','15 min','60 min','stable')
axis([2 10 0 60])


% figure 4. cell cycle duration vs. mu
figure(5)
for p = 1:counter
    
    if mod(p,2) == 0
        h(p) = errorbar(summaryMus(p),summaryDurations(p),summarySems_durations(p),'o','Color',[0.25 0.25 0.9],'MarkerSize',10);
        text(summaryMus(p)+shift,summaryDurations(p)+shift, summaryDates(p),'FontSize',10);
        axis([0 4 1 5])
        hold on
   
    elseif summaryTimescales{p} == 30
        h(p) = errorbar(summaryMus(p),summaryDurations(p),summarySems_durations(p),'o','Color',[1 0 0],'MarkerSize',10); % red
        text(summaryMus(p)+shift,summaryDurations(p)+shift, summaryDates(p),'FontSize',10);
        hold on

    elseif summaryTimescales{p} == 300
        h(p) = errorbar(summaryMus(p),summaryDurations(p),summarySems_durations(p),'o','Color',[1 0.85 0.01],'MarkerSize',10); % sunflower yellow
        text(summaryMus(p)+shift,summaryDurations(p)+shift, summaryDates(p),'FontSize',10);
        hold on
        
    elseif summaryTimescales{p} == 900
        h(p) = errorbar(summaryMus(p),summaryDurations(p),summarySems_durations(p),'o','Color',[0 0.7 0.7],'MarkerSize',10); % green
        text(summaryMus(p)+shift,summaryDurations(p)+shift, summaryDates(p),'FontSize',10);
        hold on

    elseif summaryTimescales{p} == 3600
        h(p) = errorbar(summaryMus(p),summaryDurations(p),summarySems_durations(p),'o','Color',[1 0.5 0.5],'MarkerSize',10); % peach
        text(summaryMus(p)+shift,summaryDurations(p)+shift, summaryDates(p),'FontSize',10);
        hold on
    
    end
    
end
xlabel('doubling rate of volume (1/hr)')
ylabel('cell cycle duration (min)')
legend('30 sec','5 min','15 min','60 min','stable')
axis([0.5 3 0 60])

% figure 5.  mu vs bvpr
figure(5)
for p = 1:counter
    
    if mod(p,2) == 0
        h(p) = plot(summaryBVPRs(p),summaryMus(p),'o','Color',[0.25 0.25 0.9],'MarkerSize',10);
        text(summaryBVPRs(p)+shift,summaryMus(p)+shift, summaryDates(p),'FontSize',10);
        axis([0 4 1 5])
        hold on
   
    elseif summaryTimescales{p} == 30
        h(p) = plot(summaryBVPRs(p),summaryMus(p),'o','Color',[1 0 0],'MarkerSize',10); % red
        text(summaryBVPRs(p)+shift,summaryMus(p)+shift, summaryDates(p),'FontSize',10);
        hold on

    elseif summaryTimescales{p} == 300
        h(p) = plot(summaryBVPRs(p),summaryMus(p),'o','Color',[1 0.85 0.01],'MarkerSize',10); % sunflower yellow
        text(summaryBVPRs(p)+shift,summaryMus(p)+shift, summaryDates(p),'FontSize',10);
        hold on
        
    elseif summaryTimescales{p} == 900
        h(p) = plot(summaryBVPRs(p),summaryMus(p),'o','Color',[0 0.7 0.7],'MarkerSize',10); % green
        text(summaryBVPRs(p)+shift,summaryMus(p)+shift, summaryDates(p),'FontSize',10);
        hold on

    elseif summaryTimescales{p} == 3600
        h(p) = plot(summaryBVPRs(p),summaryMus(p),'o','Color',[1 0.5 0.5],'MarkerSize',10); % peach
        text(summaryBVPRs(p)+shift,summaryDurations(p)+shift, summaryDates(p),'FontSize',10);
        hold on
    
    end
    
end
xlabel('biovolume production rate (cubic cm/hr)')
ylabel('doubling rate of volume (1/hr)')
legend('30 sec','5 min','15 min','60 min','stable')
axis([2 10 0.5 3])

% figure 6.  bvpr vs mu
figure(6)
for p = 1:counter
    
    if mod(p,2) == 0
        h(p) = plot(summaryMus(p),summaryBVPRs(p),'o','Color',[0.25 0.25 0.9],'MarkerSize',10);
        text(summaryMus(p)+shift,summaryBVPRs(p)+shift, summaryDates(p),'FontSize',10);
        axis([0 4 1 5])
        hold on
   
    elseif summaryTimescales{p} == 30
        h(p) = plot(summaryMus(p),summaryBVPRs(p),'o','Color',[1 0 0],'MarkerSize',10); % red
        text(summaryMus(p)+shift,summaryBVPRs(p)+shift, summaryDates(p),'FontSize',10);
        hold on

    elseif summaryTimescales{p} == 300
        h(p) = plot(summaryMus(p),summaryBVPRs(p),'o','Color',[1 0.85 0.01],'MarkerSize',10); % sunflower yellow
        text(summaryMus(p)+shift,summaryBVPRs(p)+shift, summaryDates(p),'FontSize',10);
        hold on
        
    elseif summaryTimescales{p} == 900
        h(p) = plot(summaryMus(p),summaryBVPRs(p),'o','Color',[0 0.7 0.7],'MarkerSize',10); % green
        text(summaryMus(p)+shift,summaryBVPRs(p)+shift, summaryDates(p),'FontSize',10);
        hold on

    elseif summaryTimescales{p} == 3600
        h(p) = plot(summaryMus(p),summaryBVPRs(p),'o','Color',[1 0.5 0.5],'MarkerSize',10); % peach
        text(summaryMus(p)+shift,summaryBVPRs(p)+shift, summaryDates(p),'FontSize',10);
        hold on
    
    end
    
end
xlabel('doubling rate of volume (1/hr)')
ylabel('biovolume production rate (cubic cm/hr)')
legend('30 sec','5 min','15 min','60 min','stable')
axis([0.5 3 2 10])


