% figure 24

%  Goal: histograms of instantaneous dV/dt that go into each bin of figure
%        23, an effort to determine source of consistent dips in 15 and 60
%        min replicates

%  Strategy:


%  Last edit: jen, 2018 May 18

%  commit: violin plots with dV/dt calculation performed BEFORE 3hr trim




% OK let's go!

%% (A) initialize analysis
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
exptsToInclude = [6,7,10:13];
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
    clear curveIDs
    
    
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
    
    dV_raw = [NaN; diff(volumes)];
    dVdt = dV_raw/dt * 3600;                    % final units = cubic um/sec
    dVdt(isDrop == 1) = NaN;
    
    clear curveFinder isDrop volumes

    
    % 9. isolate data to stabilized regions of growth
    minTime = 3;  % hr
    maxTime = bubbletime(condition);

    times_trim1 = timestamps(timestamps >= minTime);
    conditionData_trim1 = conditionData_fullOnly(timestamps >= minTime,:);
    dVdt_trim1 = dVdt(timestamps >= minTime,:);
    
    if maxTime > 0
        conditionData_fullOnly = conditionData_trim1(times_trim1 <= maxTime,:);
        dVdt_trim2 = dVdt_trim1(times_trim1 <= maxTime,:);
    else
        conditionData_fullOnly = conditionData_trim1;
        dVdt_trim2 = dVdt_trim1;
    end
    clear times_trim1 dVdt_trim1 timestamps minTime maxTime 
    
    
    % 10. isolate corrected timestamp
    if strcmp(date, '2017-10-10') == 1
        correctedTime = conditionData_fullOnly(:,2);
    else
        correctedTime = conditionData_fullOnly(:,25); % col 25 = timestamps corrected for signal lag
    end
    clear D5 M M_va T xy_start xy_end xys
    clear isDrop timestamps dV_raw firstFullCurve firstFullCurve_timestamps
    
    
    % 10. compute nutrient signal, where 1 = high and 0 = low
    timeInPeriods = correctedTime/timescale; % unit = sec/sec
    timeInPeriodFraction = timeInPeriods - floor(timeInPeriods);

    % 11. assign corrected timestamps to bins, by which to accumulate volume and dV/dt data
    timeInPeriodFraction_inSeconds = timeInPeriodFraction * timescale;
    timeInPeriodFraction_inBins = ceil(timeInPeriodFraction_inSeconds/timePerBin);

    
    % 12. remove data associated with NaN (these are in dVdt as birth events)
    growthData = [timeInPeriodFraction_inBins dVdt_trim2];
    growthData_nans = growthData(isnan(dVdt_trim2),:);
    growthData_none = growthData(~isnan(dVdt_trim2),:);
    
    % 13. collectdV/dt data into bins
    binned_dVdt = accumarray(growthData_none(:,1),growthData_none(:,2),[],@(x) {x});
    
    
    % 14. plot
%     shapes = {'o','*','square'};
    if timescale == 300
        sp = 1;
%         color_high = rgb('DarkSlateBlue');
%         color_low = rgb('DarkMagenta');
    elseif timescale == 900
        sp = 2;
%         color_high = rgb('Aquamarine');
%         color_low = rgb('Teal');
    else
        sp = 3;
%         color_high = rgb('Chocolate');
%         color_low = rgb('DodgerBlue');
    end
%     
    
    % plot violin plots across transition
    figure(sp)
    if ec == 1
        distributionPlot(binned_dVdt,'widthDiv',[2 1],'histOri','left','color',[0 0.7 0.7],'showMM',2) % green
    else
        distributionPlot(gca,binned_dVdt,'widthDiv',[2 2],'histOri','right','color',[0.25 0.25 0.9],'showMM',2) % purple
    end
    xlabel('bins (25 sec each)')
    ylabel('dV/dt (cubic um/hr)')
    title('histograms of growth rates per period bin: does not include division events')
    if timescale == 300
        legend('left: 2017-11-15','right: 2018-01-11')
    elseif timescale == 900
        legend('left: 2018-01-16','right: 2018-01-17')
    elseif timescale == 3600
        legend('left: 2018-01-29','right: 2018-01-31')
    end
    
    
    clearvars -except dVdtData_fullOnly_newdVdt storedMetaData ec timePerBin datesForLegend dataIndex exptsToInclude
    
end



