% figure 23

%  Goal: plot growth rate vs nutrient phase,
%        binning rates by period fraction (20th of a period)




%  Strategy:
%
%  Part A:
%     0. initialize analysis parameters
%     0. initialize complete meta data

%  Part B:
%     1. for all experiments in dataset:
%           2. initialize experiment meta data
%           3. load measured data
%           4. gather specified condition data
%           5. isolate parameters for growth rate calculations
%           6. calculate growth rate
%           7. isolate data to stabilized regions of growth
%           8. isolate selected specific growth rate
%           9. bin growth rates by 20th of period
%          10. calculate average growth rate per timebin
%          11. plot
%    12. repeat for all experiments



%  Last edit: jen, 2018 November 17
%  commit: repeat data plotted at period fraction = 1 at 0 as well


% Okie, go go let's go!

%% (A) initialize analysis
clc
clear

% 0. initialize analysis parameters
condition = 1; % 1 = fluctuating; 3 = ave nutrient condition
binsPerPeriod = 20;


% 0. define growth rate of interest
prompt = 'Enter specific growth rate definition as string (raw / norm / log2 / lognorm): ';
specificGrowthRate = input(prompt);


% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')

exptArray = [5,6,7,9,10,11,12,13,14,15]; % list experiments by index


%% (B) bin volumes into 20th of period timescale

% 1. for all experiments in dataset
exptCounter = 0;
for e = 1:length(exptArray)
       
    % 2. initialize experiment meta data
    index = exptArray(e);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    expType = storedMetaData{index}.experimentType;
    bubbletime = storedMetaData{index}.bubbletime;
    xys = storedMetaData{index}.xys;
    
    disp(strcat(date, ': analyze!'))
    
    exptCounter = exptCounter + 1;
    datesForLegend{exptCounter} = date;

    
    % 3. load measured data
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-width1p7-jiggle-0p5.mat');
    load(filename,'D','D5','T');
    
    
    % 4. gather specified condition data
    condition = 1; % 1 = fluctuating; 3 = ave nutrient condition
    xy_start = storedMetaData{index}.xys(condition,1);
    xy_end = storedMetaData{index}.xys(condition,end);
    conditionData = buildDM(D5, T, xy_start, xy_end, index, expType);
    clear D D5 T xy_start xy_end xys
     
    
    % 5. isolate parameters for growth rate calculations
    volumes = conditionData(:,11);        % col 11 = calculated va_vals (cubic um)
    timestamps_sec = conditionData(:,2);  % col 2  = timestamp in seconds
    isDrop = conditionData(:,4);          % col 4  = isDrop, 1 marks a birth event
    curveFinder = conditionData(:,5);     % col 5  = curve finder (ID of curve in condition)
    trackNum = conditionData(:,20);       % col 20 = track number (not ID from particle tracking)
    
    
    % 6. calculate growth rate
    growthRates = calculateGrowthRate(volumes,timestamps_sec,isDrop,curveFinder,trackNum);
    clear trackNum curveFinder isDrop volumes
    
    
    % 7. isolate data to stabilized regions of growth
    %    NOTE: errors (excessive negative growth rates) occur at trimming
    %          point if growth rate calculation occurs AFTER time trim.
    
    minTime = 3;  % hr
    maxTime = bubbletime(condition); % limit analysis to whole integer # of periods
    timestamps_hr = timestamps_sec/3600; % time in seconds converted to hours
    
    % trim to minumum
    times_trim1 = timestamps_hr(timestamps_hr >= minTime);
    conditionData_trim1 = conditionData(timestamps_hr >= minTime,:);
    growthRates_trim1 = growthRates(timestamps_hr >= minTime,:);
    
    % trim to maximum
    if maxTime > 0
        conditionData_trim2 = conditionData_trim1(times_trim1 <= maxTime,:);
        growthRates_trim2 = growthRates_trim1(times_trim1 <= maxTime,:);
    else
        conditionData_trim2 = conditionData_trim1;
        growthRates_trim2 = growthRates_trim1;
    end
    clear times_trim1 timestamps_hr minTime maxTime timestamps_sec
    clear growthRates conditionData
   

    
    % 8. isolate selected specific growth rate
    if strcmp(specificGrowthRate,'raw') == 1
        specificColumn = 1;         % for selecting appropriate column in growthRates
    elseif strcmp(specificGrowthRate,'norm') == 1
        specificColumn = 2;
    elseif strcmp(specificGrowthRate,'log2') == 1
        specificColumn = 3;
    elseif strcmp(specificGrowthRate,'lognorm') == 1
        specificColumn = 4;
    end
    
    growthRt = growthRates_trim2(:,specificColumn);
    clear specificColumn
    
    
    
    
    % 9. bin growth rates by 20th of period, first assigning growth rate to
    %    time value of the middle of two timepoints (not end value)!
    timestep_sec = 60+57;
    timeInSeconds = conditionData_trim2(:,22);  % col 22 = signal corrected time
    if strcmp(date,'2017-10-10') == 1
        timeInSeconds = conditionData_trim2(:,2);
    end
    
    timeInSeconds_middle = timeInSeconds - (timestep_sec/2);
    timeInPeriods = timeInSeconds_middle/timescale;    % units = sec/sec
    timeInPeriods_floors = floor(timeInPeriods);
    timeInPeriodFraction = timeInPeriods - timeInPeriods_floors;
    assignedBin = timeInPeriodFraction * binsPerPeriod;
    assignedBin = ceil(assignedBin);

%     timeInSeconds = conditionData_trim2(:,22);  % col 22 = signal corrected time
%     if strcmp(date,'2017-10-10') == 1
%         timeInSeconds = conditionData_trim2(:,2);
%     end
%     timeInPeriods = timeInSeconds/timescale;    % units = sec/sec
%     timeInPeriods_floors = floor(timeInPeriods);
%     timeInPeriodFraction = timeInPeriods - timeInPeriods_floors;
%     assignedBin = timeInPeriodFraction * binsPerPeriod;
%     assignedBin = ceil(assignedBin);
    
    growthRt_binned = accumarray(assignedBin,growthRt,[],@(x) {x});
    clear timeInSeconds timeInPeriods timeInPeriodFraction assignedBin
    
    
    % 10. calculate average growth rate per timebin
    %growthRt_means(exptCounter,:) = cellfun(@nanmean,growthRt_binned);
    growthRt_means = cellfun(@nanmean,growthRt_binned);
    growthRt_means_zeroed(exptCounter,:) = [growthRt_means(end); growthRt_means];
    
    % 11. plot
    if timescale == 30
        color = rgb('FireBrick');
    elseif timescale == 300
        color = rgb('Gold');
    elseif timescale == 900
        color = rgb('MediumSeaGreen');
    elseif timescale == 3600
        color = rgb('MediumSlateBlue');
    end
    
    shape = '.';

    figure(1)
    plot(growthRt_means_zeroed(exptCounter,:),'Color',color,'Marker',shape)
    hold on
    grid on
    title('growth rate: log2 mean')
    xlabel('period bin (1/20)')
    ylabel('growth rate (1/hr)')
    axis([1,21,-1,4])
    legend(datesForLegend)
    

% 12. repeat for all experiments
end
