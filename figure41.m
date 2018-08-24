%% figure 41
%
%  Goals: compare growth rate response from cells experiencing a single nutrient shift
%         to responses from cells experiencing repeated nutrient shifts. differences suggest
%         distinct physiologies.

%         plots growth rate over time, with t = 0 as shift time

%         outputs three figures:
%               1. two subplots: one each for upshift and downshift
%               2. upshift: each timescale with its own subplot
%               3. downshift: each timescale with its own subplot



%  General strategy:
%
%         Part A. initialize folder with stored meta data
%         Part B. upshift or downshift




%  last updated: jen, 2018 August 24

%  commit: edit to use new buildDM (exptType) and calculateGrowthRate
%          functions, and to plot both upshift and downshift cases


% OK let's go!

%% Part A. initialize

clc
clear

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')



% 0. define shift type, growth rate and time bin of interest
prompt = 'Enter shift type as string (upshift / downshift): ';
shiftType = input(prompt);

prompt = 'Enter specific growth rate definition as string (raw / norm / log / lognorm / mu): ';
specificGrowthRate = input(prompt);

prompt = 'Enter time per bin in seconds as double (i.e. 25): ';
timePerBin = input(prompt);


%% Part B. overlay up or downshift data


% 1. create array of experiments of interest, then loop through each

if strcmp(shiftType,'upshift');
    exptArray = [5,6,7,10:15,21,22]; % use corresponding dataIndex values
else
    exptArray = [5,6,7,10:15,26,27];
end


for e = 1:length(exptArray)
    
    % 2. initialize experiment meta data
    index = exptArray(e);                               % previous, dataIndex(e);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    bubbletime = storedMetaData{index}.bubbletime;
    expType = storedMetaData{index}.experimentType;
    shiftTime = storedMetaData{index}.shiftTime;        % sec

    disp(strcat(date, ': analyze!'))
    
    
    
    % 3. load measured data
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    if isnan(shiftTime) == 1
        filename = strcat('lb-fluc-',date,'-window5-width1p4-1p7-jiggle-0p5.mat');
    else
        filename = strcat('lb-fluc-',date,'-window5-width1p7-jiggle-0p5.mat');
    end
    load(filename,'D5','M','M_va','T');
    
    
    
    
    % 4. specify condition of interest (fluc) and build data matrix
    condition = 1;                                      % 1 = fluctuating
    xy_start = storedMetaData{index}.xys(condition,1);
    xy_end = storedMetaData{index}.xys(condition,end);
    conditionData = buildDM(D5, M, M_va, T, xy_start, xy_end,index,expType);
    
    
    
    
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
    
    if isnan(shiftTime) == 1
        minTime = 3; % hr
    else
        minTime = 2.5;
    end
    
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
        xmin = -5;                  % lower limit for plotting x axis
        xmax = 25;                  % upper limit for plotting x axis
    elseif strcmp(specificGrowthRate,'norm') == 1
        specificColumn = 2;
        xmin = -1;
        xmax = 3;
    elseif strcmp(specificGrowthRate,'log') == 1
        specificColumn = 3;
        xmin = -2;
        xmax = 3;
    elseif strcmp(specificGrowthRate,'lognorm') == 1
        specificColumn = 4;
        xmin = -0.5;
        xmax = 1;
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
    clear growthRt correctedTime
    
    
    
    
    % 12. compute nutrient signal, where 1 = high and 0 = low
    %       (i) translate timestamps into quarters of nutrient signal, such
    %           that each timepoint has a value between 1-4, which
    %           conveys the quartile of the nutrient signal in which that
    %           timepoint resides
    
    %           in original fluctuating experiments:
    %           Q1 & Q4 are part of the high nutrient phase, whereas
    %           Q2 & Q3 are part of the low
    
    timeInPeriods = correctedTime_noNans/timescale; % unit = sec/sec
    timeInPeriodFraction = timeInPeriods - floor(timeInPeriods);

    
    
    
    
    % 13. assign corrected timestamps to bins, by which to accumulate growth data
    
    if isnan(shiftTime) == 1 % if original fluctuating experiment
        
        timeInPeriodFraction_inSeconds = timeInPeriodFraction * timescale;
        bins = ceil(timeInPeriodFraction_inSeconds/timePerBin);
        
        % 14. find which bins are boundaries signal phases 
        lastBin_Q1 = (timescale/timePerBin)/4;                      % last bin before downshift
        firstBin_downshift = (timescale/4)/timePerBin + 1;
        lastBin_downshift = (timescale*3/4)/timePerBin;             % last bin before upshift
        
        firstBin_single_shift = (timescale*3/4)/timePerBin + 1;          % first bin of upshift
        lastBin_ofPeriod = timescale/timePerBin;                    % total bins in signal period
        
        
        % 15. list bins chronologically to combine broken up high nutrient phase
        %       i.e. start of upshift is Q4, concatenate Q1
        downshiftBins = firstBin_downshift:lastBin_downshift;
        upshiftBins = [firstBin_single_shift:lastBin_ofPeriod, 1:lastBin_Q1];
        clear timeInPeriodFraction timeInPeriodFraction_inSeconds
        
        
        
    else
        
        bins = ceil(correctedTime_noNans/timePerBin);      % bin 1 = first 25 sec of experiment
        bins_unique = unique(bins);
        
        % generalized for single shift experiments
        single_shiftBins = bins(bins*timePerBin > shiftTime);
        single_shiftBins_unique = unique(single_shiftBins);
        firstBin_single_shift = single_shiftBins_unique(1);    
        
        
    end
    
    
    
    % 15. choose which pre-shift data bins to plot
    if isnan(shiftTime) == 1
        
        if length(upshiftBins) >= 5 % upshift used here, but it is equal in length to downshift
            preShift_bins = 4;
        else
            preShift_bins = 2;
        end
        
        % shorter timescales (less bins) require pulling from Q4 growth data, in order to have 5 pre-shift points
        if lastBin_Q1 - preShift_bins <= 0
            first_preUPshiftBin = lastBin_Q1 - preShift_bins + lastBin_ofPeriod;
            pre_downshiftBins = [first_preUPshiftBin,lastBin_ofPeriod,1:lastBin_Q1];
        else
            % otherwise no need to tap into Q4 data
            pre_upshiftBins = lastBin_downshift - preShift_bins : lastBin_downshift;
            pre_downshiftBins = lastBin_Q1 - preShift_bins : lastBin_Q1;
        end
         
    else
        
        % single shift experiments don't have high/low phase interruptions
        preShift_bins = 10;
        
        index_single_shift = find(bins_unique == firstBin_single_shift);
        pre_upshiftBins = bins_unique(index_single_shift-preShift_bins-1 : index_single_shift-1); % same bins in both single down and upshifts
        pre_downshiftBins = bins_unique(index_single_shift-preShift_bins-1 : index_single_shift-1);
        
    end

    


    % 16. collect growth rate data into bins and calculate stats
    if isnan(shiftTime) == 1
        
        binned_growthRate = accumarray(bins,growthRt_noNaNs,[],@(x) {x});
        binned_mean = accumarray(bins,growthRt_noNaNs,[],@mean);
        binned_std = accumarray(bins,growthRt_noNaNs,[],@std);
        binned_counts = accumarray(bins,growthRt_noNaNs,[],@length);
        binned_sems = binned_std./sqrt(binned_counts);
        
    else
        
        binned_growthRate = accumarray(bins,growthRt_noNaNs,[],@(x) {x});
        binned_mean = accumarray(bins,growthRt_noNaNs,[],@mean);
        binned_std = accumarray(bins,growthRt_noNaNs,[],@std);
        binned_counts = accumarray(bins,growthRt_noNaNs,[],@length);
        binned_sems = binned_std./sqrt(binned_counts);
        
    end
    clear bins
    
    
    
    
    % 17. plot response in growth rate for all timescales over time
    if timescale == 300
        sp = 1;
        color_high = rgb('DarkSlateBlue');
        color_low = rgb('DarkMagenta');
    elseif timescale == 900
        sp = 2;
        color_high = rgb('Aquamarine');
        color_low = rgb('Teal');
    elseif timescale == 3600
        sp = 3;
        color_high = rgb('Chocolate');
        color_low = rgb('DodgerBlue');
    else
        color_high = rgb('MediumVioletRed');
        color_low = rgb('Pink');
    end
    
    
    if strcmp(shiftType,'upshift') == 1
        
        figure(1)   % upshift
        
        % pre upshift
        plot((preShift_bins*-1:0)*timePerBin,binned_mean(pre_upshiftBins),'Color',color_low,'LineWidth',1)
        hold on
        
        % post upshift
        if isnan(shiftTime) == 1
            plot((1:length(binned_mean(upshiftBins)))*timePerBin,binned_mean(upshiftBins),'Color',color_high,'LineWidth',1)
        else
            plot((1:length(binned_mean(single_shiftBins_unique)))*timePerBin,binned_mean(single_shiftBins_unique),'Color',color_high,'LineWidth',1)
        end
        grid on
        hold on
        title(strcat('response to upshift, binned every (',num2str(timePerBin),') sec'))
        
    else
        
        figure(2)    % downshift
        
        % pre downshift
        plot((preShift_bins*-1:0)*timePerBin,binned_mean(pre_downshiftBins),'Color',color_low,'LineWidth',1)
        hold on
        
        % post shift
        if isnan(shiftTime) == 1
            plot((1:length(binned_mean(downshiftBins)))*timePerBin,binned_mean(downshiftBins),'Color',color_high,'LineWidth',1)
        else
            plot((1:length(binned_mean(single_shiftBins_unique)))*timePerBin,binned_mean(single_shiftBins_unique),'Color',color_high,'LineWidth',1)
        end
        grid on
        hold on
        title(strcat('response to downshift, binned every (',num2str(timePerBin),') sec'))
        
    end

     
end


xlabel('time (sec)')
ylabel(strcat('growth rate: (', specificGrowthRate ,')'))
axis([preShift_bins*-1*timePerBin,3000,xmin,xmax])




