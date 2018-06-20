%% figure 41
%
%  Goals: overlay growth rate over time from timescale upshifts and control
%         upshift experiment

%
%         Part 1: dV/dt vs time
%         Part 2: (dV/dt)/V vs time
% 

%  General strategy:
%
%         0. initialize folder with stored meta data
%         1. run through upshift code from figures 23 (raw) or 27 (norm)
%         2. run through upshift code from figure 37 (control experiment)



%  last updated: jen, 2018 June 20

%  commit: overlay growth rate over time from timescale upshifts and control
%         upshift experiment, plot both raw and normalized dV/dt


% OK let's go!

%% Part One. raw dV/dt

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
    timestamps_sec = conditionData_fullOnly(:,2);      % col 2  = timestamp in seconds
    isDrop = conditionData_fullOnly(:,5);          % col 5  = isDrop, 1 marks a birth event
    curveFinder = conditionData_fullOnly(:,6);     % col 6  = curve finder (ID of curve in condition)
    
 
    % 8. calculate mean timestep and dVdt    
    curveIDs = unique(curveFinder);
    firstFullCurve = curveIDs(2);
    if length(firstFullCurve) > 1
        firstFullCurve_timestamps = timestamps_sec(curveFinder == firstFullCurve);
    else
        firstFullCurve = curveIDs(3);
        firstFullCurve_timestamps = timestamps_sec(curveFinder == firstFullCurve);
    end
    dt = mean(diff(firstFullCurve_timestamps)); % timestep in seconds
    
    dV_raw = [NaN; diff(volumes)];
    dVdt = dV_raw/dt * 3600;                    % final units = cubic um/sec
    
    dVdt(isDrop == 1) = NaN;

    
    
    % 9. isolate data to stabilized regions of growth
    minTime = 3;  % hr
    maxTime = bubbletime(condition);
    timestamps_hr = timestamps_sec/3600;
    
    times_trim1 = timestamps_hr(timestamps_hr >= minTime);
    conditionData_trim1 = conditionData_fullOnly(timestamps_hr >= minTime,:);
    dVdt_trim1 = dVdt(timestamps_hr >= minTime,:);
    
    if maxTime > 0
        conditionData_trim2 = conditionData_trim1(times_trim1 <= maxTime,:);
        dVdt_trim2 = dVdt_trim1(times_trim1 <= maxTime,:);
    else
        conditionData_trim2 = conditionData_trim1;
        dVdt_trim2 = dVdt_trim1;
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
    bins_byUpshift = timeInPeriodFraction * timescale;
    timeInPeriodFraction_inBins = ceil(bins_byUpshift/timePerBin);
    
    
    %firstBinAfter_downshift = (timescale/4)/timePerBin + 1;
    lastBin_downshift = (timescale*3/4)/timePerBin;
    firstBinAfter_upshift = (timescale*3/4)/timePerBin + 1;
    lastBin_ofPeriod = timescale/timePerBin;
    lastBin_Q1 = (timescale/timePerBin)/4;
    
    %downshiftBins = firstBinAfter_downshift:lastBin_downshift;
    upshiftBins = [firstBinAfter_upshift:lastBin_ofPeriod, 1:lastBin_Q1];
    
    % determine how many bins of pre-shift data to plot
    if length(upshiftBins) >= 5
        preShift_bins = 4;
    else
        preShift_bins = 2;
    end
    
    pre_upshiftBins = lastBin_downshift - preShift_bins : lastBin_downshift;
    

    
    % 12. remove data associated with NaN (these are in dVdt as birth events)    
    growthData = [timeInPeriodFraction_inBins dVdt_trim2];
    
    growthData_nans = growthData(isnan(dVdt_trim2),:);
    growthData_none = growthData(~isnan(dVdt_trim2),:);


    
    % 13. collect dV/dt data into bins and calculate stats
    binned_dVdt = accumarray(growthData_none(:,1),growthData_none(:,2),[],@(x) {x});
    binned_dVdt_mean = accumarray(growthData_none(:,1),growthData_none(:,2),[],@mean);
    binned_dVdt_std = accumarray(growthData_none(:,1),growthData_none(:,2),[],@std);
    binned_dVdt_counts = accumarray(growthData_none(:,1),growthData_none(:,2),[],@length);
    binned_dVdt_sems = binned_dVdt_std./sqrt(binned_dVdt_counts);
    
    
    
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
    
    
    % overlay of all experiments, dV/dt and sem
    figure(2)
    plot((preShift_bins*-1:0)*timePerBin,binned_dVdt_mean(pre_upshiftBins),'Color',color_low,'LineWidth',1)
    hold on
    plot((1:length(binned_dVdt_mean(upshiftBins)))*timePerBin,binned_dVdt_mean(upshiftBins),'Color',color_high,'LineWidth',1)
    grid on
    hold on
    title(strcat('upshift: mean dV/dt, binned every (',num2str(timePerBin),') sec'))
    xlabel('time (sec)')
    ylabel('dV/dt, unsynchronized')
    axis([preShift_bins*-1*timePerBin,1800,-10,25])
     
end



clear
clc

% 15. load control experiment data
load('storedMetaData_controls.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData_controls));
exptsToInclude = 1;

timePerBin = 25; % sec

for i = 1:length(exptsToInclude)
    
    % 16. collect experiment date
    e = exptsToInclude(i);
    index = dataIndex(e);
    date = storedMetaData_controls{index}.date;
    shiftTime = storedMetaData_controls{index}.shiftTime; % in sec
    disp(strcat(date, ': analyze!'))

    
    % 17. initialize experiment meta data
    bubbletime = storedMetaData_controls{index}.bubbletime;
    
    
    % 18. load measured data
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-window5-width1p4-1p7-jiggle-0p5.mat');
    load(filename,'D5','M','M_va','T');
    
    
    % 19. build data matrix from specified condition
    condition = 1; % 1 = fluctuating; 3 = ave nutrient condition
    xy_start = storedMetaData_controls{index}.xys(condition,1);
    xy_end = storedMetaData_controls{index}.xys(condition,end);
    conditionData = buildDM(D5, M, M_va, T, xy_start, xy_end,e);
    
    
    % 20. isolate condition data to those with full cell cycles
    curveIDs = conditionData(:,6);           % col 6 = curve ID
    conditionData_fullOnly = conditionData(curveIDs > 0,:);
    clear curveFinder
    
    
    % 21. isolate volume (Va) and timestamp data
    volumes = conditionData_fullOnly(:,12);        % col 12 = calculated va_vals (cubic um)
    timestamps_sec = conditionData_fullOnly(:,2);  % col 2  = timestamp in seconds
    isDrop = conditionData_fullOnly(:,5);          % col 5  = isDrop, 1 marks a birth event
    curveFinder = conditionData_fullOnly(:,6);     % col 6  = curve finder (ID of curve in condition)
    
 
    % 22. calculate mean timestep and dVdt    
    curveIDs = unique(curveFinder);
    firstFullCurve = curveIDs(2);
    if length(firstFullCurve) > 1
        firstFullCurve_timestamps = timestamps_sec(curveFinder == firstFullCurve);
    else
        firstFullCurve = curveIDs(3);
        firstFullCurve_timestamps = timestamps_sec(curveFinder == firstFullCurve);
    end
    dt = mean(diff(firstFullCurve_timestamps)); % timestep in seconds
    
    dV_raw = [NaN; diff(volumes)];
    dVdt = dV_raw/dt * 3600;                    % final units = cubic um/sec
    
    dVdt(isDrop == 1) = NaN;

    
    
    % 23. isolate data to stabilized regions of growth
    minTime = 3;  % hr
    maxTime = bubbletime(condition);
    timestamps_hr = timestamps_sec/3600;
    
    times_trim1 = timestamps_hr(timestamps_hr >= minTime);
    conditionData_trim1 = conditionData_fullOnly(timestamps_hr >= minTime,:);
    dVdt_trim1 = dVdt(timestamps_hr >= minTime,:);
    
    if maxTime > 0
        conditionData_trim2 = conditionData_trim1(times_trim1 <= maxTime,:);
        dVdt_trim2 = dVdt_trim1(times_trim1 <= maxTime,:);
    else
        conditionData_trim2 = conditionData_trim1;
        dVdt_trim2 = dVdt_trim1;
    end
    clear times_trim1 timestamps_sec timestamps_hr minTime maxTime bubbletime
    
    
     % 24. isolate corrected timestamp
    if strcmp(date, '2018-06-15') == 1
        correctedTime = conditionData_trim2(:,2);
    else
        correctedTime = conditionData_trim2(:,25); % col 25 = timestamps corrected for signal lag
    end
    clear D5 M M_va T xy_start xy_end xys
    clear isDrop timestamps dV_raw firstFullCurve firstFullCurve_timestamps
    
    
    
    % 25. compute nutrient signal, where 1 = high and 0 = low
    binaryNutrientSignal = zeros(length(correctedTime),1);
    binaryNutrientSignal(correctedTime > shiftTime) = 1;
    
    
    % 26. assign corrected timestamps to bins, by which to accumulate dV/dt data
    bins = ceil(correctedTime/timePerBin);      % bin 1 = first 25 sec of experiment
    bins_unique = unique(bins);
    
    upshiftBins = bins(bins*25 > shiftTime);
    upshiftBins_unique = unique(upshiftBins);
    firstBinAfter_upshift = upshiftBins_unique(1);
    
    % determine how many bins of pre-shift data to plot
    if length(upshiftBins_unique) >= 5
        preShift_bins = 4;
    else
        preShift_bins = 2;
    end
    
    index_upshift = find(bins_unique == firstBinAfter_upshift);
    pre_upshiftBins = bins_unique(index_upshift-preShift_bins-1 : index_upshift-1);
    

    
    % 27. remove data associated with NaN (these are in dVdt as birth events)    
    growthData = [bins dVdt_trim2];
    
    growthData_nans = growthData(isnan(dVdt_trim2),:);
    growthData_none = growthData(~isnan(dVdt_trim2),:);

    
    
    % 28. collect dV/dt data into bins and calculate stats    
    binned_dVdt = accumarray(growthData_none(:,1),growthData_none(:,2),[],@(x) {x});
    binned_dVdt_mean = accumarray(growthData_none(:,1),growthData_none(:,2),[],@mean);
    binned_dVdt_std = accumarray(growthData_none(:,1),growthData_none(:,2),[],@std);
    binned_dVdt_counts = accumarray(growthData_none(:,1),growthData_none(:,2),[],@length);
    binned_dVdt_sems = binned_dVdt_std./sqrt(binned_dVdt_counts);
    
    
    
    % 29. plot
    color_high = rgb('MediumVioletRed');
    color_low = rgb('Pink');
 
    % overlay of all experiments, dV/dt and sem
    figure(2)
    plot((preShift_bins*-1:0)*timePerBin,binned_dVdt_mean(pre_upshiftBins),'Color',color_low,'LineWidth',1)
    hold on
    plot((1:length(binned_dVdt_mean(upshiftBins_unique)))*timePerBin,binned_dVdt_mean(upshiftBins_unique),'Color',color_high,'LineWidth',1)
    grid on
    hold on
    title(strcat('upshift: mean dV/dt, binned every (',num2str(timePerBin),') sec'))
    xlabel('time (sec)')
    ylabel('dV/dt (cubic um/sec)')
    axis([preShift_bins*-1*timePerBin,4000,-10,25])
    
  
end


%% Part Two. dV/dt over V

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
    timestamps_sec = conditionData_fullOnly(:,2);      % col 2  = timestamp in seconds
    isDrop = conditionData_fullOnly(:,5);          % col 5  = isDrop, 1 marks a birth event
    curveFinder = conditionData_fullOnly(:,6);     % col 6  = curve finder (ID of curve in condition)
    
 
    % 8. calculate mean timestep and dVdt    
    curveIDs = unique(curveFinder);
    firstFullCurve = curveIDs(2);
    if length(firstFullCurve) > 1
        firstFullCurve_timestamps = timestamps_sec(curveFinder == firstFullCurve);
    else
        firstFullCurve = curveIDs(3);
        firstFullCurve_timestamps = timestamps_sec(curveFinder == firstFullCurve);
    end
    dt = mean(diff(firstFullCurve_timestamps)); % timestep in seconds
    
    dV_raw_noNan = diff(volumes);
    dV_norm = [NaN; dV_raw_noNan./volumes(1:end-1)];
    dVdt_overV = dV_norm/dt * 3600;                    % final units = cubic um/hr
    
    dVdt_overV(isDrop == 1) = NaN;

    
    
    % 9. isolate data to stabilized regions of growth
    minTime = 3;  % hr
    maxTime = bubbletime(condition);
    timestamps_hr = timestamps_sec * 3600;
    
    times_trim1 = timestamps_hr(timestamps_hr >= minTime);
    conditionData_trim1 = conditionData_fullOnly(timestamps_hr >= minTime,:);
    dVdt_overV_trim1 = dVdt_overV(timestamps_hr >= minTime,:);
    
    if maxTime > 0
        conditionData_trim2 = conditionData_trim1(times_trim1 <= maxTime,:);
        dVdt_overV_trim2 = dVdt_overV_trim1(times_trim1 <= maxTime,:);
    else
        conditionData_trim2 = conditionData_trim1;
        dVdt_overV_trim2 = dVdt_overV_trim1;
    end
    clear times_trim1 timestamps_hr timestamps_sec minTime maxTime bubbletime
    
    
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
    
    
    
    % 12. assign corrected timestamps to bins, by which to accumulate volume and dV/dt data
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
    growthData = [timeInPeriodFraction_inBins dVdt_overV_trim2];
    
    growthData_nans = growthData(isnan(dVdt_overV_trim2),:);
    growthData_none = growthData(~isnan(dVdt_overV_trim2),:);
   
    
    % 13. collect volume and dV/dt data into bins and calculate stats
    binned_dVdt = accumarray(growthData_none(:,1),growthData_none(:,2),[],@(x) {x});
    binned_dVdt_mean = accumarray(growthData_none(:,1),growthData_none(:,2),[],@mean);
    binned_dVdt_std = accumarray(growthData_none(:,1),growthData_none(:,2),[],@std);
    binned_dVdt_counts = accumarray(growthData_none(:,1),growthData_none(:,2),[],@length);
    binned_dVdt_sems = binned_dVdt_std./sqrt(binned_dVdt_counts);
    
    
    % 14. plot
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
    
    
    % overlay of all experiments, dV/dt over V
    figure(1)
    plot((preShift_bins*-1:0)*timePerBin,binned_dVdt_mean(pre_upshiftBins),'Color',color_low,'LineWidth',1)
    hold on
    plot((1:length(binned_dVdt_mean(upshiftBins)))*timePerBin,binned_dVdt_mean(upshiftBins),'Color',color_high,'LineWidth',1)
    grid on
    hold on
    title(strcat('upshift: mean (dV/dt)/V, binned every (',num2str(timePerBin),') sec'))
    xlabel('time (sec)')
    ylabel('(dV/dt)/V, unsynchronized')
    axis([preShift_bins*-1*timePerBin,1800,-2,6])
    
    
    clearvars -except dVdtData_fullOnly_newdVdt storedMetaData ec timePerBin datesForLegend dataIndex exptsToInclude
    
end

%%


clear
clc

% 15. load control experiment data
load('storedMetaData_controls.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData_controls));
exptsToInclude = 1;

timePerBin = 25; % sec

for i = 1:length(exptsToInclude)
    
    % 16. collect experiment date
    e = exptsToInclude(i);
    index = dataIndex(e);
    date = storedMetaData_controls{index}.date;
    shiftTime = storedMetaData_controls{index}.shiftTime; % in sec
    disp(strcat(date, ': analyze!'))

    
    % 17. initialize experiment meta data
    bubbletime = storedMetaData_controls{index}.bubbletime;
    
    
    % 18. load measured data
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-window5-width1p4-1p7-jiggle-0p5.mat');
    load(filename,'D5','M','M_va','T');
    
    
    % 19. build data matrix from specified condition
    condition = 1; % 1 = fluctuating; 3 = ave nutrient condition
    xy_start = storedMetaData_controls{index}.xys(condition,1);
    xy_end = storedMetaData_controls{index}.xys(condition,end);
    conditionData = buildDM(D5, M, M_va, T, xy_start, xy_end,e);
    
    
    % 20. isolate condition data to those with full cell cycles
    curveIDs = conditionData(:,6);           % col 6 = curve ID
    conditionData_fullOnly = conditionData(curveIDs > 0,:);
    clear curveFinder
    
    
    % 21. isolate volume (Va) and timestamp data
    volumes = conditionData_fullOnly(:,12);        % col 12 = calculated va_vals (cubic um)
    timestamps_sec = conditionData_fullOnly(:,2);  % col 2  = timestamp in seconds
    isDrop = conditionData_fullOnly(:,5);          % col 5  = isDrop, 1 marks a birth event
    curveFinder = conditionData_fullOnly(:,6);     % col 6  = curve finder (ID of curve in condition)
    
 
    % 22. calculate mean timestep and dVdt    
    curveIDs = unique(curveFinder);
    firstFullCurve = curveIDs(2);
    if length(firstFullCurve) > 1
        firstFullCurve_timestamps = timestamps_sec(curveFinder == firstFullCurve);
    else
        firstFullCurve = curveIDs(3);
        firstFullCurve_timestamps = timestamps_sec(curveFinder == firstFullCurve);
    end
    dt = mean(diff(firstFullCurve_timestamps)); % timestep in seconds
    
    dV_raw_noNan = diff(volumes);
    dV_norm = [NaN; dV_raw_noNan./volumes(1:end-1)];
    
    dVdt_overV = dV_norm/dt * 3600;             % final units = cubic um/hr
    dVdt_overV(isDrop == 1) = NaN;              % replace all growth rates at division events with NaN


    
    % 23. isolate data to stabilized regions of growth
    minTime = 3;  % hr
    maxTime = bubbletime(condition);
    timestamps_hr = timestamps_sec/3600;
    
    times_trim1 = timestamps_hr(timestamps_hr >= minTime);
    conditionData_trim1 = conditionData_fullOnly(timestamps_hr >= minTime,:);
    dVdt_overV_trim1 = dVdt_overV(timestamps_hr >= minTime,:);
    
    if maxTime > 0
        conditionData_trim2 = conditionData_trim1(times_trim1 <= maxTime,:);
        dVdt_overV_trim2 = dVdt_overV_trim1(times_trim1 <= maxTime,:);
    else
        conditionData_trim2 = conditionData_trim1;
        dVdt_overV_trim2 = dVdt_overV_trim1;
    end
    clear times_trim1 timestamps_sec timestamps_hr minTime maxTime bubbletime
    
    
     % 24. isolate corrected timestamp
    if strcmp(date, '2018-06-15') == 1
        correctedTime = conditionData_trim2(:,2);
    else
        correctedTime = conditionData_trim2(:,25); % col 25 = timestamps corrected for signal lag
    end
    clear D5 M M_va T xy_start xy_end xys
    clear isDrop timestamps dV_raw_noNan dV_norm firstFullCurve firstFullCurve_timestamps
    
    
    
    % 25. compute nutrient signal, where 1 = high and 0 = low
    binaryNutrientSignal = zeros(length(correctedTime),1);
    binaryNutrientSignal(correctedTime > shiftTime) = 1;
    
    
    % 26. assign corrected timestamps to bins, by which to accumulate dV/dt data
    bins = ceil(correctedTime/timePerBin);      % bin 1 = first 25 sec of experiment
    bins_unique = unique(bins);
    
    upshiftBins = bins(bins*25 > shiftTime);
    upshiftBins_unique = unique(upshiftBins);
    firstBinAfter_upshift = upshiftBins_unique(1);
    
    % determine how many bins of pre-shift data to plot
    if length(upshiftBins_unique) >= 5
        preShift_bins = 4;
    else
        preShift_bins = 2;
    end
    
    index_upshift = find(bins_unique == firstBinAfter_upshift);
    pre_upshiftBins = bins_unique(index_upshift-preShift_bins-1 : index_upshift-1);
    

    
    % 27. remove data associated with NaN (these are in dVdt as birth events)    
    growthData = [bins dVdt_overV_trim2];
    
    growthData_nans = growthData(isnan(dVdt_overV_trim2),:);
    growthData_none = growthData(~isnan(dVdt_overV_trim2),:);

    
    
    % 28. collect dV/dt data into bins and calculate stats    
    binned_dVdt = accumarray(growthData_none(:,1),growthData_none(:,2),[],@(x) {x});
    binned_dVdt_mean = accumarray(growthData_none(:,1),growthData_none(:,2),[],@mean);
    binned_dVdt_std = accumarray(growthData_none(:,1),growthData_none(:,2),[],@std);
    binned_dVdt_counts = accumarray(growthData_none(:,1),growthData_none(:,2),[],@length);
    binned_dVdt_sems = binned_dVdt_std./sqrt(binned_dVdt_counts);
    
    
    
    % 29. plot
    color_high = rgb('MediumVioletRed');
    color_low = rgb('Pink');
 
    % overlay of all experiments, dV/dt and sem
    figure(1)
    plot((preShift_bins*-1:0)*timePerBin,binned_dVdt_mean(pre_upshiftBins),'Color',color_low,'LineWidth',1)
    hold on
    plot((1:length(binned_dVdt_mean(upshiftBins_unique)))*timePerBin,binned_dVdt_mean(upshiftBins_unique),'Color',color_high,'LineWidth',1)
    grid on
    hold on
    title(strcat('upshift: mean dV/dt normalizeb by initial volume, binned every (',num2str(timePerBin),') sec'))
    xlabel('time (sec)')
    ylabel('(dV/dt)/ V (1/sec)')
    axis([preShift_bins*-1*timePerBin,4000,-10,25])
    
  
end

