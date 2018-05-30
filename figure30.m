% figure 30

% goal: are negative dV/dt values real? what are they associated with?

%       1. are there nutrient period bins that have a higher frequency of
%          negatives than others?
%
%       2. what is the average value of these negatives?


% last update: jen, 2018 May 30
% commit: edit comments to reflect true units

% OK let's go!

%%
clc
clear

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
load('dVdtData_fullOnly_newdVdt.mat')

dataIndex = find(~cellfun(@isempty,storedMetaData));


% 1. for all experiments in dataset
timePerBin = 25; % sec
ec = 0; % experiment counter
exptsToInclude = [6,7,10:14];

%%
for i = 1:length(exptsToInclude)
    
    % 2. collect experiment date
    e = exptsToInclude(i);
    index = dataIndex(e);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    

    disp(strcat(date, ': analyze!'))
    ec = ec + 1;
    
    % reset experiment counter at start of new timescale
    if e == 10 || e == 12
        ec = 1;
    end
    
    
    % 3. initialize experiment meta data
    xys = storedMetaData{index}.xys;
    bubbletime = storedMetaData{index}.bubbletime;
    
    
    % 4. load measured data
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-window5-width1p4-1p7-jiggle-0p5.mat');
    load(filename,'D5','M','M_va','T');
    
    
    % 5. build data matrix from specified condition
    condition = 1; % 1 = fluctuating; 3 = ave nutrient condition
    xy_start = storedMetaData{index}.xys(condition,1);
    xy_end = storedMetaData{index}.xys(condition,end);
    conditionData = buildDM(D5, M, M_va, T, xy_start, xy_end,e);
    
    
    % 6. isolate condition data to those with full cell cycles
    curveIDs = conditionData(:,6);           % col 6 = curve ID
    conditionData_fullOnly = conditionData(curveIDs > 0,:);
    clear curveFinder
    
    
    % 7. isolate volume (Va) and timestamp data
    volumes = conditionData_fullOnly(:,12);        % col 12 = calculated va_vals (cubic um)
    timestamps = conditionData_fullOnly(:,2);      % col 2  = timestamp in seconds
    isDrop = conditionData_fullOnly(:,5);          % col 5  = isDrop, 1 marks a birth event
    curveFinder = conditionData_fullOnly(:,6);     % col 6  = curve finder (ID of curve in condition)
    
 
    % 8. calculate mean timestep and dVdt    
    curveIDs = unique(curveFinder);
    firstFullCurve = curveIDs(2);
    if length(firstFullCurve) > 1
        firstFullCurve_timestamps = timestamps(curveFinder == firstFullCurve);
    else
        firstFullCurve = curveIDs(3);
        firstFullCurve_timestamps = timestamps(curveFinder == firstFullCurve);
    end
    dt = mean(diff(firstFullCurve_timestamps)); % timestep in seconds
    
    dV_raw_noNan = diff(volumes);
    dV_norm = [NaN; dV_raw_noNan./volumes(1:end-1)];
    dVdt_overV = dV_norm/dt * 3600;                    % final units = 1/hr
    
    dVdt_overV(isDrop == 1) = NaN;

    
    
    % 9. isolate data to stabilized regions of growth
    minTime = 3;  % hr
    maxTime = bubbletime(condition);
    
    times_trim1 = timestamps(timestamps >= minTime);
    conditionData_trim1 = conditionData_fullOnly(timestamps >= minTime,:);
    dVdt_overV_trim1 = dVdt_overV(timestamps >= minTime,:);
    
    if maxTime > 0
        conditionData_trim2 = conditionData_trim1(times_trim1 <= maxTime,:);
        dVdt_overV_trim2 = dVdt_overV_trim1(times_trim1 <= maxTime,:);
    else
        conditionData_trim2 = conditionData_trim1;
        dVdt_overV_trim2 = dVdt_overV_trim1;
    end
    clear times_trim1 timestamps minTime maxTime bubbletime
    
    
     % 10. isolate corrected timestamp
    if strcmp(date, '2017-10-10') == 1
        correctedTime = conditionData_trim2(:,2);
    else
        correctedTime = conditionData_trim2(:,25); % col 25 = timestamps corrected for signal lag
    end
    clear D5 M M_va T xy_start xy_end xys
    clear isDrop timestamps dV_raw firstFullCurve firstFullCurve_timestamps
    
    
    
    % 11. compute nutrient signal, where 1 = high and 0 = low
    %       (i) translate timestamps into quarters of nutrient signal
    timeInPeriods = correctedTime/timescale; % unit = sec/sec
    timeInPeriodFraction = timeInPeriods - floor(timeInPeriods);
    timeInQuarters = ceil(timeInPeriodFraction * 4);
    
    %      (ii) from nutrient signal quarters, generate a binary nutrient signal where, 1 = high and 0 = low
    binaryNutrientSignal = zeros(length(timeInQuarters),1);
    binaryNutrientSignal(timeInQuarters == 1) = 1;
    binaryNutrientSignal(timeInQuarters == 4) = 1;
    
    
    % 11. assign corrected timestamps to bins, by which to accumulate volume and dV/dt data
    timeInPeriodFraction_inSeconds = timeInPeriodFraction * timescale;
    timeInPeriodFraction_inBins = ceil(timeInPeriodFraction_inSeconds/timePerBin);
    
    
    firstBinAfter_downshift = (timescale/4)/timePerBin + 1;
    lastBin_downshift = (timescale*3/4)/timePerBin;
    firstBinAfter_upshift = (timescale*3/4)/timePerBin + 1;
    lastBin_ofPeriod = timescale/timePerBin;
    lastBin_Q1 = (timescale/timePerBin)/4;
    
    downshiftBins = firstBinAfter_downshift:lastBin_downshift;
    upshiftBins = [firstBinAfter_upshift:lastBin_ofPeriod, 1:lastBin_Q1];
    
    % determine how many bins of pre-shift data to plot
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
    
    
    
    % 12. remove data associated with NaN (these are in dVdt as birth events)
    growthData = [curveFinder timeInPeriodFraction_inBins volumes dVdt_overV_trim2];
    
    dVdt_negs = dVdt_overV_trim2;
    dVdt_negs(dVdt_overV_trim2 < 0) = NaN;
    growthData_nans = growthData(isnan(dVdt_negs),:);
    growthData_none = growthData(~isnan(dVdt_negs),:);

    nanReporter = size(growthData_nans)
    valuesReporter = size(growthData_none)
    
    
    % 13. collect dV/dt data into bins and calculate stats
    binned_neg_means = accumarray(growthData_nans(:,2),growthData_nans(:,4),[],@nanmean);
    binned_neg_stds = accumarray(growthData_nans(:,2),growthData_nans(:,4),[],@nanstd);
    
    binned_neg_counts = accumarray(growthData_nans(:,2),growthData_nans(:,4),[],@length);
    binned_none_counts = accumarray(growthData_none(:,2),growthData_none(:,4),[],@length);
    binned_total_counts = binned_neg_counts + binned_none_counts;
    binned_fraction_neg = binned_neg_counts./binned_total_counts;
   
    
    % 14. plot
    shapes = {'o','*','square'};
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
    
    
    % counts across bin
    figure(1)
    subplot(2,1,1) % upshift
    plot((preShift_bins*-1:0)*timePerBin,binned_fraction_neg(pre_upshiftBins),'Color',color_low,'LineWidth',1)
    hold on
    plot((1:length(binned_fraction_neg(upshiftBins)))*timePerBin,binned_fraction_neg(upshiftBins),'Color',color_high,'LineWidth',1)
    grid on
    hold on
    title(strcat('upshift: fraction of negative values, binned every (',num2str(timePerBin),') sec'))
    xlabel('time (sec)')
    ylabel('% negative dV/dt')
    axis([preShift_bins*-1*timePerBin,1800,0,1])
    
    subplot(2,1,2) % downshift
    plot((preShift_bins*-1:0)*timePerBin,binned_fraction_neg(pre_downshiftBins),'Color',color_high,'LineWidth',1)
    hold on
    plot((1:length(binned_fraction_neg(downshiftBins)))*timePerBin,binned_fraction_neg(downshiftBins),'Color',color_low,'LineWidth',1)
    grid on
    hold on
    title(strcat('downshift: fraction of negative values, binned every (',num2str(timePerBin),') sec'))
    xlabel('time (sec)')
    ylabel('% negative dV/dt')
    axis([preShift_bins*-1*timePerBin,1800,0,1])
    
    
    % mean negative value across bin
    figure(2)
    subplot(2,1,1) % upshift
    errorbar((preShift_bins*-1:0)*timePerBin,binned_neg_means(pre_upshiftBins),binned_neg_stds(pre_upshiftBins),'Color',color_low,'LineWidth',1,'Marker','.')
    hold on
    errorbar((1:length(binned_neg_means(upshiftBins)))*timePerBin,binned_neg_means(upshiftBins),binned_neg_stds(upshiftBins),'Color',color_high,'LineWidth',1,'Marker','.')
    grid on
    hold on
    title(strcat('upshift: mean negative value + std dev, binned every (',num2str(timePerBin),') sec'))
    xlabel('time (sec)')
    ylabel('mean dV/dt')
    axis([preShift_bins*-1*timePerBin,1800,-2,1])
    
    subplot(2,1,2) % downshift
    errorbar((preShift_bins*-1:0)*timePerBin,binned_neg_means(pre_downshiftBins),binned_neg_stds(pre_downshiftBins),'Color',color_high,'LineWidth',1,'Marker','.')
    hold on
    errorbar((1:length(binned_neg_means(downshiftBins)))*timePerBin,binned_neg_means(downshiftBins),binned_neg_stds(downshiftBins),'Color',color_low,'LineWidth',1,'Marker','.')
    grid on
    hold on
    title(strcat('downshift: mean negative value + std dev, binned every (',num2str(timePerBin),') sec'))
    xlabel('time (sec)')
    ylabel('mean dV/dt')
    axis([preShift_bins*-1*timePerBin,1800,-2,1])
    
 
    clearvars -except dVdtData_fullOnly_newdVdt storedMetaData ec timePerBin datesForLegend dataIndex exptsToInclude
    
end