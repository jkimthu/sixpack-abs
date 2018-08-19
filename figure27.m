% figure 27

%  Goal: plot timescale-specific responses to repeated nutrient shifts,
%        gathering all data (after first 3 hours) and binning by time after
%        shift.

%        compares response between timescales. differences suggest
%        distinct physiologies.

%        outputs three figures:
%               1. two subplots: one each for upshift and downshift
%               2. upshift: each timescale with its own subplot
%               3. downshift: each timescale with its own subplot
%
%        for comparison to single upshift, see Figure 41


%  Strategy:
%       0. initialize complete meta data
%       0. define growth rate and time bin of interest
%       1. create array of experiments of interest, then loop through each:
%               2. initialize experiment meta data
%               3. load measured experiment data
%               4. specify condition of interest (fluc) and build data matrix
%               5. isolate condition data to those from full cell cycles
%               6. isolate data to stabilized regions of growth
%               7. isolate volume (Va), timestamp, mu, drop and curveID data
%               8. calculate growth rate
%               9. isolate selected specific growth rate
%              10. isolate corrected timestamp
%              11. remove nans from data analysis
%              12. compute nutrient signal, where 1 = high and 0 = low
%                     (i) translate timestamps into quarters of nutrient signal
%                    (ii) from nutrient signal quarters, generate a binary nutrient signal where, 1 = high and 0 = low
%              13. assign corrected timestamps to bins, by which to accumulate growth rate data
%              14. list bins belonging to high or low phase chronogically
%              15. choose how many bins of pre-shift data to plot
%              16. accumulate growth data into bins and calculate stats
%              17. plot!
%      18. repeat for all indeces in experiments-of-interest array



%  Last edit: jen, 2018 Aug 19


%  commit: repaired error due to calculating growth after trimming by time 


% OK let's go!

%% (A) initialize analysis
clc
clear

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')

%dataIndex = find(~cellfun(@isempty,storedMetaData));


% 0. define growth rate and time bin of interest
prompt = 'Enter specific growth rate definition as string (raw / norm / log / lognorm / mu): ';
specificGrowthRate = input(prompt);

prompt = 'Enter time per bin in seconds as double (i.e. 25): ';
timePerBin = input(prompt);



% 1. create array of experiments of interest, then loop through each:
exptArray = [5,6,7,10:15]; % use corresponding dataIndex values


%%
for e = 1:length(exptArray)
    
    % 2. initialize experiment meta data
    index = exptArray(e); %dataIndex(e);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    bubbletime = storedMetaData{index}.bubbletime;
    expType = storedMetaData{index}.experimentType;
    
    disp(strcat(date, ': analyze!'))


    
    % 3. load measured experiment data
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-window5-width1p4-1p7-jiggle-0p5.mat');
    load(filename,'D5','M','M_va','T');
    
    
    
    % 4. specify condition of interest (fluc) and build data matrix
    condition = 1;                % 1 = fluctuating
    xy_start = storedMetaData{index}.xys(condition,1);
    xy_end = storedMetaData{index}.xys(condition,end);
    conditionData = buildDM(D5, M, M_va, T, xy_start, xy_end,index,expType);
    clear xy_start xy_end
    
    
    
    % 5. isolate condition data to those from full cell cycles
    curveIDs = conditionData(:,6);           % col 6 = curve ID
    conditionData_fullOnly = conditionData(curveIDs > 0,:);
    clear curveFinder
    
    
    
    % 6. isolate volume (Va), timestamp, mu, drop and curveID data
    volumes = conditionData_fullOnly(:,12);            % col 12 = calculated va_vals (cubic um)
    timestamps_sec = conditionData_fullOnly(:,2);      % col 2  = timestamp in seconds
    isDrop = conditionData_fullOnly(:,5);              % col 5  = isDrop, 1 marks a birth event
    curveFinder = conditionData_fullOnly(:,6);         % col 6  = curve finder (ID of curve in condition)
    mus = conditionData_fullOnly(:,14);           % col 14 = mu, calculated from volume tracks
    
    
    
    % 7. calculate growth rate
    growthRates = calculateGrowthRate(volumes,timestamps_sec,isDrop,curveFinder,mus);
    
    
    
    
    % 8. isolate data to stabilized regions of growth
    %    NOTE: errors (excessive negative growth rates) occur at trimming
    %          point if growth rate calculation occurs AFTER time trim.
    %
    minTime = 3;  % hr
    maxTime = bubbletime(condition);
    timestamps_sec = conditionData_fullOnly(:,2); % time in seconds converted to hours
    timestamps_hr = timestamps_sec / 3600;
    
    % trim to minumum
    times_trim1 = timestamps_hr(timestamps_hr >= minTime);
    conditionData_trim1 = conditionData_fullOnly(timestamps_hr >= minTime,:);
    growthRates_trim1 = growthRates(timestamps_hr >= minTime,:);
    
    % trim to maximum
    if maxTime > 0
        conditionData_trim2 = conditionData_trim1(times_trim1 <= maxTime,:);
        growthRates_trim2 = growthRates_trim1(times_trim1 <= maxTime,:);
    else
        conditionData_trim2 = conditionData_trim1;
        growthRates_trim2 = growthRates_trim1;
    end
    clear growthRates conditionData_fullOnly

    
     
    % 9. isolate selected specific growth rate
    if strcmp(specificGrowthRate,'raw') == 1
        specificColumn = 1;         % for selecting appropriate column in growthRates
    elseif strcmp(specificGrowthRate,'norm') == 1
        specificColumn = 2;
    elseif strcmp(specificGrowthRate,'log') == 1
        specificColumn = 3;
    elseif strcmp(specificGrowthRate,'lognorm') == 1
        specificColumn = 4;
    elseif strcmp(specificGrowthRate,'mu') == 1;
        specificColumn = 5;
    end
    
    growthRt = growthRates_trim2(:,specificColumn);
    
    
    

    % 10. isolate corrected timestamp
    if strcmp(date, '2017-10-10') == 1
        correctedTime = conditionData_trim2(:,2);
    else
        correctedTime = conditionData_trim2(:,25); % col 25 = timestamps corrected for signal lag
    end
    clear D5 M M_va T isDrop timestamps_sec   
    
    
    
    % 11. remove nans from data analysis
    growthRt_noNaNs = growthRt(~isnan(growthRt),:);
    correctedTime_noNans = correctedTime(~isnan(growthRt),:);

   
        
    % 12. compute nutrient signal, where 1 = high and 0 = low
    %       (i) translate timestamps into quarters of nutrient signal
    timeInPeriods = correctedTime_noNans/timescale;        % unit = sec/sec
    timeInPeriodFraction = timeInPeriods - floor(timeInPeriods);
    timeInQuarters = ceil(timeInPeriodFraction * 4);
    
    %      (ii) from nutrient signal quarters, generate a binary nutrient signal where, 1 = high and 0 = low
    binaryNutrientSignal = zeros(length(timeInQuarters),1);
    binaryNutrientSignal(timeInQuarters == 1) = 1;
    binaryNutrientSignal(timeInQuarters == 4) = 1;
    
    
    
    % 13. assign corrected timestamps to bins, by which to accumulate growth rate data
    timeInPeriodFraction_inSeconds = timeInPeriodFraction * timescale;
    timeInPeriodFraction_inBins = ceil(timeInPeriodFraction_inSeconds/timePerBin);
    
    
    
    % 14. list bins belonging to high or low phase chronogically
    firstBinAfter_downshift = (timescale/4)/timePerBin + 1;
    lastBin_downshift = (timescale*3/4)/timePerBin;
    firstBinAfter_upshift = (timescale*3/4)/timePerBin + 1;
    lastBin_ofPeriod = timescale/timePerBin;
    lastBin_Q1 = (timescale/timePerBin)/4;
    
    downshiftBins = firstBinAfter_downshift:lastBin_downshift;
    upshiftBins = [firstBinAfter_upshift:lastBin_ofPeriod, 1:lastBin_Q1];
    
    
    
    % 15. choose how many bins of pre-shift data to plot
    if length(upshiftBins) >= 5
        preShift_bins = 4;
    else
        preShift_bins = 2;
    end
    
    % shorter timescales (less bins) require pulling from Q4 growth data,
    % in order to have 5 pre-shift points
    if lastBin_Q1 - preShift_bins <= 0
        first_preshift = lastBin_Q1 - preShift_bins + lastBin_ofPeriod;
        pre_downshiftBins = [first_preshift,lastBin_ofPeriod,1:lastBin_Q1];
    else
        % otherwise no need to tap into Q4 data
        pre_downshiftBins = lastBin_Q1 - preShift_bins : lastBin_Q1;
    end
    
    pre_upshiftBins = lastBin_downshift - preShift_bins : lastBin_downshift;
    
    
    
    
    % 16. accumulate growth data into bins and calculate stats
    binned_growthRates = accumarray(timeInPeriodFraction_inBins,growthRt_noNaNs,[],@(x) {x});
    binned_mean = accumarray(timeInPeriodFraction_inBins,growthRt_noNaNs,[],@mean);
    binned_std = accumarray(timeInPeriodFraction_inBins,growthRt_noNaNs,[],@std);
    binned_counts = accumarray(timeInPeriodFraction_inBins,growthRt_noNaNs,[],@length);
    binned_sems = binned_std./sqrt(binned_counts);
    
    
    
    % 17. plot
    if timescale == 300
        sp = 1;
        color_high = rgb('DarkSlateBlue');
        color_low = rgb('DarkMagenta');
    elseif timescale == 900
        sp = 2;
        color_high = rgb('Aquamarine');
        color_low = rgb('Teal');
    else
        sp = 3;
        color_high = rgb('Chocolate');
        color_low = rgb('DodgerBlue');
    end
    
    
    % overlay of all experiments
    figure(1)
    subplot(2,1,1) % upshift
    plot((preShift_bins*-1:0)*timePerBin,binned_mean(pre_upshiftBins),'Color',color_low,'LineWidth',1)
    hold on
    plot((1:length(binned_mean(upshiftBins)))*timePerBin,binned_mean(upshiftBins),'Color',color_high,'LineWidth',1)
    grid on
    hold on
    title(strcat('response to upshift: binned every (',num2str(timePerBin),') sec'))
    xlabel('time (sec)')
    ylabel(strcat('growth rate: (',specificGrowthRate,')'))
    axis([preShift_bins*-1*timePerBin,1800,-2,6])
    
    subplot(2,1,2) % downshift
    plot((preShift_bins*-1:0)*timePerBin,binned_mean(pre_downshiftBins),'Color',color_high,'LineWidth',1)
    hold on
    plot((1:length(binned_mean(downshiftBins)))*timePerBin,binned_mean(downshiftBins),'Color',color_low,'LineWidth',1)
    grid on
    hold on
    title(strcat('response to downshift: binned every (',num2str(timePerBin),') sec'))
    xlabel('time (sec)')
    ylabel(strcat('growth rate: (',specificGrowthRate,')'))
    axis([preShift_bins*-1*timePerBin,1800,-2,6])
    
    
    % upshift subplots separating timescale
    figure(2)
    subplot(3,1,sp) % upshift
    plot((preShift_bins*-1:0)*timePerBin,binned_mean(pre_upshiftBins),'Color',color_low,'LineWidth',1)
    hold on
    plot((1:length(binned_mean(upshiftBins)))*timePerBin,binned_mean(upshiftBins),'Color',color_high,'LineWidth',1)
    grid on
    hold on
    title(strcat(num2str(timescale),': response to upshift'))
    xlabel('time (sec)')
    ylabel(strcat('growth rate: (',specificGrowthRate,')'))
    axis([preShift_bins*-1*timePerBin,1800,-2,6])
    
    
    % downshift subplots separating timescale, dV/dt and sem
    figure(3)
    subplot(3,1,sp) % downshift
    plot((preShift_bins*-1:0)*timePerBin,binned_mean(pre_downshiftBins),'Color',color_high,'LineWidth',1)
    hold on
    plot((1:length(binned_mean(downshiftBins)))*timePerBin,binned_mean(downshiftBins),'Color',color_low,'LineWidth',1)
    grid on
    hold on
    title(strcat(num2str(timescale),': response to downshift'))
    xlabel('time (sec)')
    ylabel(strcat('growth rate: (',specificGrowthRate',')'))
    axis([preShift_bins*-1*timePerBin,1800,-2,6])
    
    
    % 18. repeat for all indeces in experiments-of-interest array
end



