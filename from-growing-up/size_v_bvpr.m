% size_v_bvpr

% goal: plot cell size at birth against biovolume production rate,
%       for a specific concentration (i.e. average)

% strategy:
%
%       0. initialize data
%       0. determine target concentration
%       1. create a directory of experiments with target concentration
%       2. for all experiments in target directory... accumulate cell size stats
%               3. move to experiment folder and build data matrix
%               4. for each condition with target concentration...
%                       5. isolate data from current condition
%                       6. isolate size, drop and time data
%                       7. isolate only data during which drop == 1 (birth event)
%                       8. trim data to stabilized / non-bubble timestamps
%                       9. calculate mean and s.e.m. of size
%                      10. accumulate data for storage and plotting
%              11. store data from all conditions into measured data structure
%      12.  plot average cell size at birth against corresponding biovol production rate
%


% last updated: 2018 Jan 19
% commit: update to include new experiments, Jan 4, 11, 12, 16 and 17
%         add experimental date to point plotted  


% OK let's go!!

%% 0. initialize data

clear
clc

cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')

dataIndex = find(~cellfun(@isempty,storedMetaData));
birthSizeData = cell(size(storedMetaData));

% initialize summary vectors for calculated data
experimentCount = length(dataIndex);

% determine target concentration
targetConcentration = 0.0105; % average

%% 1. create a directory of conditions with target concentration
targetConditions = cell(1,experimentCount);

for e = 1:experimentCount
    
    % identify conditions with target concentration
    index = dataIndex(e);
    concentrations = storedMetaData{index}.concentrations;
    
    % each cell represents an experiment, each value a condition of target concentration
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
    
    % exclude outlier from analysis
%     if strcmp(date, '2017-10-31') == 1 %|| strcmp (timescale, 'monod') == 1
%         disp(strcat(date,': excluded from analysis'))
%         continue
%     end
%     disp(strcat(date, ': analyze!'))
    
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
        volumes = conditionData(:,15);        % col 15 = calculated va_vals (cubic um)
        drops = conditionData(:,5);           % col 5 = 1 at birth, zero otherwise
        timestamps = conditionData(:,2)/3600; % time in seconds converted to hours
        clear conditionData
        
        % 7. isolate only data during which drop == 1 (birth event)
        birthVolumes = volumes(drops == 1);
        birthTimes = timestamps(drops == 1);
        
        % 8. trim data to stabilized / non-bubble timestamps
        minTime = 3;  % hr
        maxTime = storedMetaData{index}.bubbletime(c);
        
        times_trim1 = birthTimes(birthTimes >= minTime);
        volumes_trim1 = birthVolumes(birthTimes >= minTime);
        clear birthTimes birthVolumes
        
        if maxTime > 0
            volumes_trim2 = volumes_trim1(times_trim1 <= maxTime);
        else
            volumes_trim2 = volumes_trim1;
        end
        clear times_trim1
        
        % 9. calculate mean and s.e.m. of size
        mean_birthSize = mean(volumes_trim2);
        count_birthSize = length(volumes_trim2);
        std_birthSize = std(volumes_trim2);
        sem_birthSize = std_birthSize./sqrt(count_birthSize);
        
        % 10. accumulate data for storage / plotting
        compiledBirthSize{c}.date = date;
        compiledBirthSize{c}.timescale = timescale;
        compiledBirthSize{c}.condition = c;
        compiledBirthSize{c}.mean = mean_birthSize;
        compiledBirthSize{c}.std = std_birthSize;
        compiledBirthSize{c}.count = count_birthSize;
        compiledBirthSize{c}.sem = sem_birthSize;
        
    end
    
    clear mean_birthSize std_birthSize count_birthSize sem_birthSize
    
    % 11. store data from all conditions into measured data structure
    birthSizeData{index} = compiledBirthSize;
    clear compiledBirthSize

end
%% 12.  plot average cell size at birth against corresponding biovol production rate

cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('bioProdRateData.mat')

% measuredData is a variable created in plotMonod.m
load('measuredData.mat')

%%
% initialize data vector for easy plotting
counter = 0;
summaryDates = cell(1,(experimentCount-1)*2);
summarySizes = zeros(1,(experimentCount-1)*2);
summaryStds = zeros(1,(experimentCount-1)*2);
summarySems = zeros(1,(experimentCount-1)*2);
summaryBVPRs = zeros(1,(experimentCount-1)*2);
summaryMus = zeros(1,(experimentCount-1)*2);
summaryTimescales = cell(1,(experimentCount-1)*2);


for e = 1:experimentCount
    
    % identify conditions with target concentration
    index = dataIndex(e);
    date = storedMetaData{index}.date;
    
    % exclude outlier from analysis
%     if strcmp(date, '2017-10-31') == 1 %|| strcmp (timescale, 'monod') == 1
%         disp(strcat(date,': excluded from analysis'))
%         continue
%     end
%     disp(strcat(date, ': analyze!'))
    
    for i = 1:length(targetConditions{e})
        c = targetConditions{e}(i);
        
        counter = counter + 1;
        summaryDates{counter} = date;
        summarySizes(counter) = birthSizeData{index}{c}.mean;
        summaryStds(counter) = birthSizeData{index}{c}.std;
        summarySems(counter) = birthSizeData{index}{c}.sem;
        summaryBVPRs(counter) = bioProdRateData{index}{c}.mean;
        summaryMus(counter) = measuredData{index}.individuals{c}.muMean;
        
        
        summaryTimescales{counter} = birthSizeData{index}{c}.timescale;
    end
end

%%
shift = 0.1;

figure(1)
for p = 1:counter
    
    if mod(p,2) == 0
        h(p) = errorbar(summaryBVPRs(p),summarySizes(p),summarySems(p),'o','Color',[0.25 0.25 0.9],'MarkerSize',10);
        text(summaryBVPRs(p)+shift,summarySizes(p)+shift, summaryDates(p),'FontSize',10);
        hold on
   
    elseif summaryTimescales{p} == 30
        h(p) = errorbar(summaryBVPRs(p),summarySizes(p),summarySems(p),'o','Color',[1 0 0],'MarkerSize',10); % red
        text(summaryBVPRs(p)+shift,summarySizes(p)+shift, summaryDates(p),'FontSize',10);
        hold on

    elseif summaryTimescales{p} == 300
        h(p) = errorbar(summaryBVPRs(p),summarySizes(p),summarySems(p),'o','Color',[1 0.85 0.01],'MarkerSize',10); % sunflower yellow
        text(summaryBVPRs(p)+shift,summarySizes(p)+shift, summaryDates(p),'FontSize',10);
        hold on
        
    elseif summaryTimescales{p} == 900
        h(p) = errorbar(summaryBVPRs(p),summarySizes(p),summarySems(p),'o','Color',[0 0.7 0.7],'MarkerSize',10); % green
        text(summaryBVPRs(p)+shift,summarySizes(p)+shift, summaryDates(p),'FontSize',10);
        hold on
        
    end
    
end
xlabel('biovol prod rate (cubic um/hr)')
ylabel('cell volume at birth (cubic um)')
legend('30 sec','5 min','15 min','stable')
axis([2 10 1 5])

figure(2)
for p = 1:counter
    
    if mod(p,2) == 0
        h(p) = errorbar(summaryMus(p),summarySizes(p),summarySems(p),'o','Color',[0.25 0.25 0.9],'MarkerSize',10);
        text(summaryBVPRs(p)+shift,summarySizes(p)+shift, summaryDates(p),'FontSize',10);
        axis([0 4 1 5])
        hold on
   
    elseif summaryTimescales{p} == 30
        h(p) = errorbar(summaryMus(p),summarySizes(p),summarySems(p),'o','Color',[1 0 0],'MarkerSize',10); % red
        text(summaryBVPRs(p)+shift,summarySizes(p)+shift, summaryDates(p),'FontSize',10);
        hold on

    elseif summaryTimescales{p} == 300
        h(p) = errorbar(summaryMus(p),summarySizes(p),summarySems(p),'o','Color',[1 0.85 0.01],'MarkerSize',10); % sunflower yellow
        text(summaryBVPRs(p)+shift,summarySizes(p)+shift, summaryDates(p),'FontSize',10);
        hold on
        
    elseif summaryTimescales{p} == 900
        h(p) = errorbar(summaryMus(p),summarySizes(p),summarySems(p),'o','Color',[0 0.7 0.7],'MarkerSize',10); % green
        text(summaryBVPRs(p)+.2,summarySizes(p)+.2, summaryDates(p),'FontSize',10);
        hold on
        
    end
    
end
xlabel('doubling rate of volume (1/hr)')
ylabel('cell volume at birth (cubic um)')
legend('30 sec','5 min','15 min','stable')

