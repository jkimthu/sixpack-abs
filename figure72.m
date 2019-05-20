%  Figure 72 - successive period overal

%  Goal: plot growth rate vs nutrient phase,
%        binning rates by period fraction (20th of a period)
%        plotting each period in succession, colored by time




%  Strategy:
%
%  Part 1:
%     0. initialize analysis parameters
%     0. initialize complete meta data

%  Part 2:
%     1. for all experiments in dataset:
%           2. initialize experiment meta data
%           3. load measured data
%           4. gather specified condition data
%           5. isolate parameters for growth rate calculations
%           6. calculate growth rate
%           7. generate a period vector by which to isolate data
%           8. isolate growth rate of interest
%           9. for each period, bin growth rates into 20ths of period timescale
%                   10. isolate data from current period,
%                   11. bin growth rates by 20th of period,
%                       assigning growth rate to time value of the middle of two timepoints (not end value)!
%                   12. calculate average growth rate per timebin
%                   13. plot
%    14. repeat for all experiments



%  Last edit: jen, 2019 May 20
%  commit: first commit, looking for stabilization of growth response

% Okie, go go let's go!

%% A. initialize analysis
clc
clear

% 0. initialize analysis parameters
condition = 1; % 1 = fluctuating; 3 = ave nutrient condition
binsPerPeriod = 20;


% 0. define growth rate of interest
specificGrowthRate = 'log2';
specificColumn = 3;


% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
exptArray = [13,14,15]; % list experiments by index


% 0. initialize colors per successive period
colorSpectrum = {'Indigo','MediumSlateBlue','DodgerBlue','DeepSkyBlue','Teal','DarkGreen','MediumSeaGreen','GoldenRod','DarkOrange','Red'};
shape = '.';


%% B. loop through each period and isolate growth rates

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

    
    
    % 3. load measured data
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
    load(filename,'D5','T');
    
        
    
    % 4. gather specified condition data
    condition = 1; % 1 = fluctuating; 3 = ave nutrient condition
    xy_start = storedMetaData{index}.xys(condition,1);
    xy_end = storedMetaData{index}.xys(condition,end);
    conditionData = buildDM(D5, T, xy_start, xy_end, index, expType);
    clear D5 T xy_start xy_end xys
    

    
    % 5. isolate parameters for growth rate calculations
    volumes = getGrowthParameter(conditionData,'volume');             % volume = calculated va_vals (cubic um)
    timestamps_sec = getGrowthParameter(conditionData,'timestamp');   % ND2 file timestamp in seconds
    isDrop = getGrowthParameter(conditionData,'isDrop');              % isDrop == 1 marks a birth event
    curveFinder = getGrowthParameter(conditionData,'curveFinder');    % col 5  = curve finder (ID of curve in condition)
    trackNum = getGrowthParameter(conditionData,'trackNum');          % track number, not ID from particle tracking
    
    
        
    % 6. calculate growth rate
    growthRates = calculateGrowthRate(volumes,timestamps_sec,isDrop,curveFinder,trackNum);
    clear trackNum curveFinder isDrop volumes
    
    
    
    % 7. generate a period vector by which to isolate data
    timestamps_corrected = getGrowthParameter(conditionData,'correctedTime'); % signal corrected time
    if strcmp(date,'2017-10-10') == 1
        timestamps_corrected = timestamps_sec;
    end
        
    timestamps_hr = timestamps_sec/3600;
    maxTime = bubbletime(condition);
    timestamps_trimmed = timestamps_corrected(timestamps_hr <= maxTime);
    growthRates_trimmed = growthRates(timestamps_hr <= maxTime,:);
    period = ceil(timestamps_trimmed/3600); % current binning only works with 60 min period data
    clear timestamps_hr timestamps_sec maxTime growthRates
    
    
  
    % 8. isolate growth rate of interest
    growthRt = growthRates_trimmed(:,specificColumn); % log2
    clear growthRates_trimmed

    
    
    % 9. for each period, bin growth rates into 20ths of period timescale
    for pp = 1:max(period)
        
        
        % 10. isolate data from current period
        %     NOTE: errors (excessive negative growth rates) occur at trimming
        %           point if growth rate calculation occurs AFTER time trim.

        periodMus = growthRt(period == pp);
        periodTimes = timestamps_trimmed(period == pp);
        
        
        % 11. bin growth rates by 20th of period,
        %     assigning growth rate to time value of the middle of two timepoints (not end value)!
        timestep_sec = 60+57;
        timeInSeconds_middle = periodTimes - (timestep_sec/2);
        timeInPeriods = timeInSeconds_middle/timescale;    % units = sec/sec
        timeInPeriods_floors = floor(timeInPeriods);
        timeInPeriodFraction = timeInPeriods - timeInPeriods_floors;
        assignedBin = timeInPeriodFraction * binsPerPeriod;
        assignedBin = ceil(assignedBin);
        
        
        periodMus_binned = accumarray(assignedBin,periodMus,[],@(x) {x});
        clear timeInPeriods timeInPeriods_floors timeInPeriodFraction assignedBin timeInSeconds_middle
        

   
        % 12. calculate average growth rate per timebin
        periodMus_means = cellfun(@nanmean,periodMus_binned);
        periodMus_means_zeroed(pp,:) = [periodMus_means(end); periodMus_means];
        
        
        
        % 13. plot
        color = rgb(colorSpectrum{pp});
        figure(e)
        plot(periodMus_means_zeroed(pp,:),'Color',color,'LineWidth',1)
        hold on
        grid on
        title('growth rate: log2 mean')
        xlabel('period bin (1/20)')
        ylabel('growth rate (1/hr)')
        axis([1,21,-1,3])
        title(date)
    
        
    end
    
    periodMus_compiled{exptCounter} = periodMus_means_zeroed;
    clear periodMus_means_zeroed
    
% 14. repeat for all experiments
end




