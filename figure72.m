%%  Figure 72 - successive period overlay

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



%  Last edit: jen, 2019 July 3
%  commit: stabilization of growth response from steady to fluc experiments

% Okie, go go let's go!

%% A. initialize analysis
clc
clear

% 0. initialize analysis parameters
condition = 1; % 1 = fluctuating; 3 = ave nutrient condition
binsPerPeriod = 30;


% 0. define growth rate of interest
specificGrowthRate = 'log2';
specificColumn = 3;


% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
exptArray = [37,38]; % list experiments by index


% 0. initialize colors per successive period
%colorSpectrum = {'Indigo','MediumSlateBlue','DodgerBlue','DeepSkyBlue','Teal','DarkGreen','MediumSeaGreen','GoldenRod','DarkOrange','Red'};
colorSpectrum = {'Indigo','DodgerBlue','DeepSkyBlue','Teal','MediumSeaGreen','GoldenRod','Red'};
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
    shiftTime = storedMetaData{index}.shiftTime;
    
    disp(strcat(date, ': analyze!'))
    exptCounter = exptCounter + 1;

    
    
    % 3. load measured data
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
    load(filename,'D5','T');
    
        
    
    % 4. gather specified condition data
    condition = 1; % 1 = fluctuating; 3 = ave nutrient condition
    if strcmp(expType,'steady2fluc') == 1
            xy_start = xys{condition}(1);
            xy_end = xys{condition}(end);
        else
            xy_start = xys(condition,1);
            xy_end = xys(condition,end);
    end
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
    
    
    
    % 7. trim timestamps and growth rate by bubble time
    %     NOTE: errors (excessive negative growth rates) occur at trimming
    %           point if growth rate calculation occurs AFTER time trim.
    timestamps_corrected = getGrowthParameter(conditionData,'correctedTime'); % signal corrected time
    if strcmp(date,'2017-10-10') == 1
        timestamps_corrected = timestamps_sec;
    end
        
    timestamps_hr = timestamps_sec/3600;
    maxTime = bubbletime(condition);
    
    if maxTime ~= 0
        timestamps_trimmed = timestamps_corrected(timestamps_hr <= maxTime);
        growthRates_trimmed = growthRates(timestamps_hr <= maxTime,:);
    else
        timestamps_trimmed = timestamps_corrected;
        growthRates_trimmed = growthRates;
    end
    
    
    
    % 8. generate a bin vector for entire time series
    %    assigning growth rate to time value of the middle of two timepoints (not end value)!
    timestep_sec = 60+57;
    timeInSeconds_middle = timestamps_trimmed - (timestep_sec/2);
    timeInPeriods = timeInSeconds_middle/timescale;
    binVector = ceil(timeInPeriods*binsPerPeriod);
    clear timestamps_hr timestamps_sec maxTime growthRates timestamps_corrected
    clear timeInPeriods timeInSeconds_middle timeInPeriods timestamps_trimmed timestep_sec
    
    
  
    % 9. isolate growth rate of interest and bin by timestamp
    growthRt = growthRates_trimmed(:,specificColumn); % log2
    growthRt(binVector == 0) = [];
    binVector(binVector == 0) = [];
    growthRt_binned = accumarray(binVector,growthRt,[],@(x) {x});
    clear growthRates_trimmed

    
    
    % 10. initialize first period, based on shift time (time of fluctuation onset)
    shiftBin = (shiftTime/timescale)*binsPerPeriod;
    period_shift = [shiftBin-1, shiftBin:shiftBin+binsPerPeriod]; % includes two bins before upshift
    period_first = period_shift - binsPerPeriod;
    clear period_shift shiftBin
    
    
    % 11. for each period, bin growth rates into period fraction
    %numPeriods = max(ceil(binVector/binsPerPeriod));
    for pp = 1:7 %numPeriods
        
        period_current = period_first + (30*(pp-1));
        
        
        % 12. isolate data from current period
        isEnd = period_current > length(growthRt_binned);
        if sum(isEnd) > 0
            period_current = period_current(isEnd == 0);
        end
        
        period_mus = growthRt_binned(period_current);

        
   
        % 13. calculate average growth rate per timebin
        signal_current = cellfun(@nanmean,period_mus);
        signals_compiled{pp,:} = signal_current;

        
        
        % 13. plot
        color = rgb(colorSpectrum{pp});
        figure(e)
        plot(signal_current,'Color',color,'LineWidth',1)
        hold on
        grid on
        title('growth rate: log2 mean')
        xlabel(strcat('period bin (',num2str(binsPerPeriod),' binsPerPeriod)'))
        ylabel('growth rate (1/hr)')
        axis([0,binsPerPeriod+3,-1,3])
        title(date)
        
        clear signal_current
    
        
    end
    
    signals_all{exptCounter} = signals_compiled;
    clear signals_compiled
    
    
% 14. repeat for all experiments
end


%% C. compile replicate data with errorbars
