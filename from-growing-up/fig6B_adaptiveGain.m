%% figure 6B - growth rate response to nutrient upshift
%
%  Output: two plots of growth rate response to upshifts in real time.
%          (1) mean of each experimental replicate, one curve each
%          (2) mean + standard dev across experimental replicates
%
%          plot compares experiments  from three conditions:
%          (i) 15-min period, (ii) 60-min period, (iii) single upshift
%
%          Also outputs two vectors, mean single shift response data and time
%          into .mat workspace titled, response_singleUpshift.mat .
%          These two vectors are generated in section C2


%  General strategy:
%
%         Part A. initialize folder with stored meta data
%         Part B. curves for fluctuating data
%         Part C. curves for single shift data

%  Note: running parts B & C can be reversed in order!
%        current order with C first plots fluc over single shift data



%  last updated: jen, 2019 July 11

%  commit: add gain from downshift, but perhaps doesn't make sense?


% OK let's go!

%% A. Gain during upshift


clc
clear


% 0. load mean growth rate signal from upshift and downshift comparisons
cd('/Users/jen/Documents/StockerLab/Writing/manuscript 1/figures/figure5/')
load('response_singleUpshift.mat')
single_upshift_means = upshift_means_single;
single_upshift_times = upshift_times_single;
clear upshift_means_single upshift_times_single


% 0. fluctuating data, upshift
load('response_flucUpshift_2.mat','upshift_times_frep','upshift_growth_frep')
fluc_upshift_mean_60min = upshift_growth_frep;
fluc_upshift_times = upshift_times_frep;
clear replicate_mean_60min upshift_times_frep



% 1. define shift type, growth rate and time bin of interest
shiftType = 'upshift';


% 2. calculate mean growth rate between time zero and t = 30 min

% single
start_s = find(single_upshift_times == 0);
final_s = find(single_upshift_times == 60);
counter = 0;
for ss = start_s+1:final_s
    counter = counter + 1;
    growth_upshift_single(counter) = mean(single_upshift_means(start_s:ss));
end
clear counter ss

growth_upshift_single_times = single_upshift_times(start_s+1:final_s);

figure(1)
hold on
plot(growth_upshift_single_times,growth_upshift_single)




% fluctuating (60 min only)
start_f = find(fluc_upshift_times == 0);
final_f = find(fluc_upshift_times == 30);
ct = 0;
for ff = start_f+1:final_f
    ct = ct + 1;
    growth_upshift_fluc(ct) = mean(fluc_upshift_mean_60min(start_f:ff));
end
clear ct ff 

growth_upshift_fluc_times = fluc_upshift_times(start_f+1:final_f);

figure(1)
hold on
plot(growth_upshift_fluc_times,growth_upshift_fluc)

clear start_s final_s start_f final_f



% 3. calculate difference between means

% i. find timepoints in single that fit within range of fluctuating data
inBoth = ismember(fluc_upshift_times,single_upshift_times);
time_both = fluc_upshift_times(inBoth==1);
overlap = time_both(time_both > 0);
clear inBoth time_both



% ii. find means from both for these timepoints
for ot = 1:length(overlap)
    
    comparable_single(ot) = growth_upshift_single(growth_upshift_single_times == overlap(ot));
    comparable_fluc(ot) = growth_upshift_fluc(growth_upshift_fluc_times == overlap(ot));
    
end

figure(1)
hold on
plot(overlap,comparable_single,'o')
hold on
plot(overlap,comparable_fluc,'o')
xlim([0 65])
xlabel('mean growth rate since shift (1/h)')
ylabel('time since shift (min)')


% iii. find means from single data extending beyond overlap
overlap_final = find(growth_upshift_single_times == overlap(ot));
comparable_single_extended = growth_upshift_single(overlap_final+1:end);
comparable_single_extended_times = growth_upshift_single_times(overlap_final+1:end);
clear ot


% iv. generate vector of 'extended' data for fluctuating condition, of
%     equal length to extended single data
extension = ones(1,length(comparable_single_extended));
comparable_fluc_extended = extension * mean(comparable_fluc(11:end));

figure(1)
hold on
plot(comparable_single_extended_times,comparable_fluc_extended,'o')


% v. calculate % difference between means of measured and extended data
for comp = 1:length(comparable_fluc)
    
    sub = comparable_fluc(comp) - comparable_single(comp);
    compared_data(comp) = sub/comparable_single(comp) * 100;
    
end
clear comp sub

for ext = 1:length(comparable_fluc_extended)
    
    sb = comparable_fluc_extended(ext) - comparable_single_extended(ext);
    compared_extended(ext) = sb/comparable_single_extended(ext) * 100;
    
end
clear ext


% vi. plot comparison
compared_extended = [compared_data'; compared_extended'];
compared_extended_times = [overlap'; comparable_single_extended_times'];

figure(2)
bar(compared_extended_times,compared_extended)
xlim([0 65])
ylabel('percent difference from steady-state')
xlabel('time since shift (min)')


%% B. Gain during downshift

clear
clc

% 0. load single downshift data
load('response_singleDownshift.mat')
single_downshift_means = downshift_means_single;
single_downshift_times = downshift_times_single;
clear downshift_means_single downshift_times_single


% 0. load fluctuating data, downshift
load('response_flucDownshift_2.mat','downshift_times_frep','downshift_growth_frep')
%fluc_downshift_mean_15min = downshift_replicate_mean_15min;
fluc_downshift_mean_60min = downshift_growth_frep;
fluc_downshift_times = downshift_times_frep;
clear downshift_replicate_mean_60min downshift_times_frep


% 1. define shift type, growth rate and time bin of interest
shiftType = 'downshift';


% 2. calculate mean growth rate between time of most neg point in single and t = 30 min

% single
vel = diff(single_downshift_means);
start_sd = find(vel == min(vel)) + 1;
final_sd = length(single_downshift_times);
%final_sd = find(single_downshift_times == 60);
counter = 0;
for ss = start_sd+1:final_sd
    counter = counter + 1;
    growth_downshift_single(counter) = nanmean(single_downshift_means(start_sd:ss));
end
clear counter ss

growth_downshift_single_times = single_downshift_times(start_sd+1:final_sd);

figure(3)
hold on
plot(growth_downshift_single_times,growth_downshift_single)



% fluctuating (60 min only)
start_fd = find(fluc_downshift_times == growth_downshift_single_times(1));
final_fd = find(fluc_downshift_times == 30);
ct = 0;
for ff = start_fd+1:final_fd
    ct = ct + 1;
    growth_downshift_fluc(ct) = mean(fluc_downshift_mean_60min(start_fd:ff));
end
clear ct ff 

%figure(1)
%plot(fluc_downshift_times,fluc_downshift_mean_60min)
%hold on
%plot(start_fd,fluc_downshift_mean_60min(start_fd),'o')
%hold on
%plot(single_downshift_times,single_downshift_means)

growth_downshift_fluc_times = fluc_downshift_times(start_fd+1:final_fd);

figure(3)
hold on
plot(growth_downshift_fluc_times,growth_downshift_fluc)

clear start_sd final_sd start_fd final_fd




% 3. calculate difference between means

% i. find timepoints in single that fit within range of fluctuating data
inBoth = ismember(growth_downshift_fluc_times,growth_downshift_single_times);
time_both = growth_downshift_fluc_times(inBoth==1);
overlap = time_both(time_both > 0);
clear inBoth time_both


% ii. find means from both for these timepoints
for ot = 1:length(overlap)
    
    comparable_single(ot) = growth_downshift_single(growth_downshift_single_times == overlap(ot));
    comparable_fluc(ot) = growth_downshift_fluc(growth_downshift_fluc_times == overlap(ot));
    
end

figure(3)
hold on
plot(overlap,comparable_single,'o')
hold on
plot(overlap,comparable_fluc,'o')
xlim([0 185])
xlabel('mean growth rate since shift (1/h)')
ylabel('time since shift (min)')


% iii. find means from single data extending beyond overlap
overlap_final = find(growth_downshift_single_times == overlap(ot));
comparable_single_extended = growth_downshift_single(overlap_final+1:end);
comparable_single_extended_times = growth_downshift_single_times(overlap_final+1:end);
clear ot


% iv. generate vector of 'extended' data for fluctuating condition, of
%     equal length to extended single data
extension = ones(1,length(comparable_single_extended));
comparable_fluc_extended = extension * mean(comparable_fluc(11:end));

figure(3)
hold on
plot(comparable_single_extended_times,comparable_fluc_extended,'o')
ylim([0.1 0.5])

% v. calculate % difference between means of measured and extended data
for comp = 1:length(comparable_fluc)
    
    sub = comparable_fluc(comp) - comparable_single(comp);
    compared_data(comp) = sub/comparable_single(comp) * 100;
    
end
clear comp sub

for ext = 1:length(comparable_fluc_extended)
    
    sb = comparable_fluc_extended(ext) - comparable_single_extended(ext);
    compared_extended(ext) = sb/comparable_single_extended(ext) * 100;
    
end
clear ext



% vi. plot comparison
compared_extended = [compared_data'; compared_extended'];
compared_extended_times = [overlap'; comparable_single_extended_times'];

figure(4)
bar(compared_extended_times,compared_extended)
xlim([0 185])
ylabel('percent difference from steady-state')
xlabel('time since shift (min)')


%% supplement 1. compile upshift data from 60 min replicates 

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')

shiftType = 'upshift';
specificGrowthRate = 'log2';
specificColumn = 3; % log2 growth rate
xmin = -0.5;
xmax = 3.5;
timePerBin = 75;

% 0. initialize color designations
color900 = 'DeepSkyBlue'; 
color3600 = 'Navy'; 
color_single = 'DarkCyan'; 

% 1. create array of experiments of interest, then loop through each
exptArray = 13:15;


counter = 0; % because timescales will differ between experiments
for e = 1:length(exptArray)
    
    counter = counter + 1;
    
    % 2. initialize experiment meta data
    index = exptArray(e);                               % previous, dataIndex(e);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    timescale_vector(counter) = timescale;
    
    bubbletime = storedMetaData{index}.bubbletime;
    expType = storedMetaData{index}.experimentType;

    disp(strcat(date, ': analyze!'))
    
    
    
    % 3. load measured data    
    if ischar(timescale) == 0
        filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
    end
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    load(filename,'D5','T')
    clear experimentFolder
    
    
    
    % 4. specify condition of interest (fluc) and build data matrix
    condition = 1;                                      % 1 = fluctuating
    xy_start = storedMetaData{index}.xys(condition,1);
    xy_end = storedMetaData{index}.xys(condition,end);
    conditionData = buildDM(D5, T, xy_start, xy_end,index,expType);
    
    
    
    % 5. isolate volume (Va), timestamp, mu, drop and curveID data
    volumes = getGrowthParameter(conditionData,'volume');             % calculated va_vals (cubic um)
    timestamps_sec = getGrowthParameter(conditionData,'timestamp');   % timestamp in seconds
    isDrop = getGrowthParameter(conditionData,'isDrop');              % isDrop, 1 marks a birth event
    curveFinder = getGrowthParameter(conditionData,'curveFinder');    % curve finder (ID of curve in condition)
    trackNum = getGrowthParameter(conditionData,'trackNum');          % track number (not ID from particle tracking)
    clear expType xy_start xy_end
    
    
    % 6. calculate growth rate
    growthRates = calculateGrowthRate(volumes,timestamps_sec,isDrop,curveFinder,trackNum);
    clear curveFinder trackNum isDrop volumes
    
    
    
    % 7. isolate data to stabilized regions of growth
    %    NOTE: errors (excessive negative growth rates) occur at trimming
    %          point if growth rate calculation occurs AFTER time trim.
    minTime = 3; % hr
    
    maxTime = bubbletime(condition);
    timestamps_sec = conditionData(:,2); % time in seconds converted to hours
    timestamps_hr = timestamps_sec / 3600;
    
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
    clear growthRates conditionData maxTime minTime timestamps_sec

    
     
    % 8. isolate selected specific growth rate
    growthRt = growthRates_trim2(:,specificColumn);
    
    
    
    % 9. isolate corrected timestamp
    if strcmp(date, '2017-10-10') == 1
        correctedTime = conditionData_trim2(:,2);
    else
        correctedTime = conditionData_trim2(:,22); % col 22 = timestamps corrected for signal lag
    end
    clear D5 T isDrop timestamps_sec  
    
    
    
    % 10. remove nans from data analysis
    growthRt_noNaNs = growthRt(~isnan(growthRt),:);
    correctedTime_noNans = correctedTime(~isnan(growthRt),:);
    clear growthRt correctedTime times_trim1
    
    
    
    % 11. assign corrected timestamps to bins, by which to accumulate growth data
    timeInPeriods = correctedTime_noNans/timescale; % unit = sec/sec
    timeInPeriodFraction = timeInPeriods - floor(timeInPeriods);
    
    timeInPeriodFraction_inSeconds = timeInPeriodFraction * timescale;
    bins = ceil(timeInPeriodFraction_inSeconds/timePerBin);
    bins_unique = unique(bins);
    clear timeInPeriods timeInPeriodFraction
    
    
    
    % 12. find which bins are boundaries signal phases
    lastBin_Q1 = (timescale/timePerBin)/4;                      % last bin before downshift
    firstBin_downshift = (timescale/4)/timePerBin + 1;
    lastBin_downshift = (timescale*3/4)/timePerBin;             % last bin before upshift
    
    firstBin_upshift = (timescale*3/4)/timePerBin + 1;          % first bin of upshift
    lastBin_ofPeriod = timescale/timePerBin;                    % total bins in signal period
    
    
    
    
    % 13. list bins chronologically to combine broken up high nutrient phase
    %       i.e. start of upshift is Q4, concatenate Q1
    downshiftBins{counter} = firstBin_downshift:lastBin_downshift;
    upshiftBins{counter} = [firstBin_upshift:lastBin_ofPeriod, 1:lastBin_Q1];
    clear timeInPeriodFraction timeInPeriodFraction_inSeconds
    
    
    

    % 14. choose which pre-shift data bins to plot

    % decide how much data to plot prior to shift
    % upshift used here, but it is equal in length to downshift
    if length(upshiftBins{counter}) >= 5
        numPreshiftBins = 4;
    else
        numPreshiftBins = 2;
    end
    preShift_bins{counter} = numPreshiftBins;
    
    
    % for downshifts
    % shorter timescales (less bins) require pulling from Q4 growth data
    if lastBin_Q1 - numPreshiftBins <= 0
        
        % absolute value of (lastBin_Q1 - preShift_bins) = number of bins needed from Q4
        % flipping list of bins from last to first lets us use absolute value as index
        bins_unique_flipped = flipud(bins_unique);
        first_pre_downshiftBin = bins_unique_flipped( abs(lastBin_Q1 - numPreshiftBins) );
        
        % array of pre-downshift bin numbers in chronological order
        pre_downshiftBins{counter} = [first_pre_downshiftBin:lastBin_ofPeriod,1:lastBin_Q1];
        
        clear bins_unique_flipped
        
    else
        
        % if no need to tap into Q4 data...
        pre_downshiftBins{counter} = lastBin_Q1 - (numPreshiftBins-1) : lastBin_Q1;
        
    end
    
    % for upshifts
    % bins are already chronological in number, so the job is easy!
    pre_upshiftBins{counter} = lastBin_downshift - (numPreshiftBins-1) : lastBin_downshift;
    


    % 16. collect growth rate data into bins and calculate stats
    binned_growthRate{counter} = accumarray(bins,growthRt_noNaNs,[],@(x) {x});
    binned_mean{counter} = accumarray(bins,growthRt_noNaNs,[],@mean);
    clear bins
    
   
    
    
    % 17. plot response in growth rate for all timescales over time
    timePerBin_min = timePerBin/60; % time per bin in minutes
    
    if timescale == 900
        color = rgb(color900);
    elseif timescale == 3600
        color = rgb(color3600);
    end
    
    
    
    % plot!
    
    %concatenate pre and post upshift data
    
    % time
    preUpshift_times = ((numPreshiftBins-1)*-1:0)*timePerBin_min;
    postUpshift_times = (1:length( binned_mean{counter}( upshiftBins{counter} ) ) )*timePerBin_min;
    upshift_times_frep = [preUpshift_times,postUpshift_times];
    
    % growth rate
    preUpshift_growth = binned_mean{counter}(pre_upshiftBins{counter});
    postUpshift_growth = binned_mean{counter}(upshiftBins{counter});
    upshift_growth_frep = [preUpshift_growth;postUpshift_growth];
    
    
    figure(1)   % mean of each replicate
    plot(upshift_times_frep,upshift_growth_frep,'Color',color,'LineWidth',1,'Marker','.')
    hold on
    
    
    clear x y ycenter ybegin yend color
    clear conditionData_trim1 conditionData_trim2 growthRates_trim1 growthRates_trim2
    clear growthRt_noNaNs index lastBin_Q1 lastBin_ofPeriod lastBin_downshift firstBin_upshift firstBin_downshift
    clear first_pre_downshiftBin date condition e
    
end

figure(1)
title(strcat('mean of each replicate: response to upshift, binned every (',num2str(timePerBin),') sec'))
xlabel('time (min)')
ylabel(strcat('growth rate: (', specificGrowthRate,')'))
axis([numPreshiftBins*-1*timePerBin_min,80,xmin,xmax])

   
for t = 2 % loop through the two timescales, 15 min and 60 min
   
    clear replicate_means r
    
    if t == 1
        timescale = 900;
        color = rgb(color900);
    else
        timescale = 3600;
        color = rgb(color3600);
    end
    
    
    % determine which columns from collected B1 data go with current timescale
    columns = find(timescale_vector == timescale);
    
    
    % extract 2d matrix of mean replicate data, with columns replicate and rows being time (period bin)
    for r = 1:length(columns)
        replicate_means(:,r) = binned_mean{columns(r)};
    end
    
    
    % calculate mean and stdev across replicates (rows)
    replicate_means_mean = mean(replicate_means,2);
    if t == 1
        replicate_mean_15min = replicate_means_mean;
    else
        replicate_mean_60min = replicate_means_mean;
    end
    replicate_means_std = std(replicate_means,0,2); % w=0 normalizes by N-1, w=1 normalizes by N
    


    
    % plot
    clear preUpshift_growth postUpshift_growth
    
    % concatenate growth rate to match time relative to shift
    preUpshift_growth = replicate_means_mean(pre_upshiftBins{columns(1)});
    postUpshift_growth = replicate_means_mean(upshiftBins{columns(1)});
    upshift_growth_fluc_means = [preUpshift_growth;postUpshift_growth];
    
    
    figure(2) % mean of replicate means
    plot(upshift_times_frep(1:length(upshift_growth_fluc_means)),upshift_growth_fluc_means,'Color',color,'LineWidth',1,'Marker','.')
    hold on
    
    
    % prepare matrix of replicate data for errorbar plotting
    preUpshift_replicates = replicate_means(pre_upshiftBins{columns(1)},:);
    postUpshift_replicates = replicate_means(upshiftBins{columns(1)},:);
    upshift_replicates = [preUpshift_replicates;postUpshift_replicates];
    
    y = upshift_replicates';
    x = upshift_times_frep(1:length(upshift_replicates))';
    
    
    figure(3) % mean and std of replicate means
    hold on
    s = shadedErrorBar(x,y,{@mean,@std},'lineprops',{'Color',color},'patchSaturation',0.4);
    
    % set face and edge properties
    s.mainLine.LineWidth = 3;
    s.patch.FaceColor = color;

end

xlabel('time (min)')
ylabel(strcat('growth rate: (', specificGrowthRate,')'))
axis([-5,120,0,3.1])

save('response_flucUpshift_2','upshift_times_frep','upshift_growth_frep')
clear r t color x y timescale bubbletime 


%% supplement 2. compile downshift data from 60 min replicates

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')

shiftType = 'downshift';
specificGrowthRate = 'log2';
specificColumn = 3; % log2 growth rate
xmin = -0.5;
xmax = 3.5;
timePerBin = 75;

% 0. initialize color designations
color900 = 'DeepSkyBlue'; 
color3600 = 'Navy'; 
color_single = 'DarkCyan'; 

% 1. create array of experiments of interest, then loop through each
exptArray = 13:15;

counter = 0;
for e = 1:length(exptArray)
    
    counter = counter + 1;
    
    % 2. initialize experiment meta data
    index = exptArray(e);                               % previous, dataIndex(e);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    timescale_vector(counter) = timescale;
    
    bubbletime = storedMetaData{index}.bubbletime;
    expType = storedMetaData{index}.experimentType;

    disp(strcat(date, ': analyze!'))
    
    
    
    % 3. load measured data
    if ischar(timescale) == 0
        filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
    end
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    load(filename,'D5','T')
    clear experimentFolder
    
    
    
    % 4. specify condition of interest (fluc) and build data matrix
    condition = 1;                                      % 1 = fluctuating
    xy_start = storedMetaData{index}.xys(condition,1);
    xy_end = storedMetaData{index}.xys(condition,end);
    conditionData = buildDM(D5, T, xy_start, xy_end,index,expType);
    
    
    
    % 5. isolate volume (Va), timestamp, mu, drop and curveID data
    volumes = getGrowthParameter(conditionData,'volume');             % calculated va_vals (cubic um)
    timestamps_sec = getGrowthParameter(conditionData,'timestamp');   % timestamp in seconds
    isDrop = getGrowthParameter(conditionData,'isDrop');              % isDrop, 1 marks a birth event
    curveFinder = getGrowthParameter(conditionData,'curveFinder');    % curve finder (ID of curve in condition)
    trackNum = getGrowthParameter(conditionData,'trackNum');          % track number (not ID from particle tracking)
    clear expType xy_start xy_end
    
    
    % 6. calculate growth rate
    growthRates = calculateGrowthRate(volumes,timestamps_sec,isDrop,curveFinder,trackNum);
    clear curveFinder trackNum isDrop volumes
    
    
    
    % 7. isolate data to stabilized regions of growth
    %    NOTE: errors (excessive negative growth rates) occur at trimming
    %          point if growth rate calculation occurs AFTER time trim.
    minTime = 3; % hr
    
    maxTime = bubbletime(condition);
    timestamps_sec = conditionData(:,2); % time in seconds converted to hours
    timestamps_hr = timestamps_sec / 3600;
    
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
    clear growthRates conditionData maxTime minTime timestamps_sec

    
     
    % 8. isolate selected specific growth rate
    growthRt = growthRates_trim2(:,specificColumn);
    
    

    
    % 9. isolate corrected timestamp
    if strcmp(date, '2017-10-10') == 1
        correctedTime = conditionData_trim2(:,2);
    else
        correctedTime = conditionData_trim2(:,22); % col 22 = timestamps corrected for signal lag
    end
    clear D5 T isDrop timestamps_sec  
    
    
    
    
    % 10. remove nans from data analysis
    growthRt_noNaNs = growthRt(~isnan(growthRt),:);
    correctedTime_noNans = correctedTime(~isnan(growthRt),:);
    clear growthRt correctedTime times_trim1
    
    
    
    % 11. assign corrected timestamps to bins, by which to accumulate growth data
    timeInPeriods = correctedTime_noNans/timescale; % unit = sec/sec
    timeInPeriodFraction = timeInPeriods - floor(timeInPeriods);
    
    timeInPeriodFraction_inSeconds = timeInPeriodFraction * timescale;
    bins = ceil(timeInPeriodFraction_inSeconds/timePerBin);
    bins_unique = unique(bins);
    clear timeInPeriods timeInPeriodFraction
    
    
    
    % 12. find which bins are boundaries signal phases
    lastBin_Q1 = (timescale/timePerBin)/4;                      % last bin before downshift
    firstBin_downshift = (timescale/4)/timePerBin + 1;
    lastBin_downshift = (timescale*3/4)/timePerBin;             % last bin before upshift
    
    firstBin_upshift = (timescale*3/4)/timePerBin + 1;          % first bin of upshift
    lastBin_ofPeriod = timescale/timePerBin;                    % total bins in signal period
    
    
    
    
    % 13. list bins chronologically to combine broken up high nutrient phase
    %       i.e. start of upshift is Q4, concatenate Q1
    downshiftBins{counter} = firstBin_downshift:lastBin_downshift;
    upshiftBins{counter} = [firstBin_upshift:lastBin_ofPeriod, 1:lastBin_Q1];
    clear timeInPeriodFraction timeInPeriodFraction_inSeconds
    
    
    

    % 14. choose which pre-shift data bins to plot

    % decide how much data to plot prior to shift
    % upshift used here, but it is equal in length to downshift
    if length(upshiftBins{counter}) >= 5
        numPreshiftBins = 4;
    else
        numPreshiftBins = 2;
    end
    preShift_bins{counter} = numPreshiftBins;
    
    
    % for downshifts
    % shorter timescales (less bins) require pulling from Q4 growth data
    if lastBin_Q1 - numPreshiftBins <= 0
        
        % absolute value of (lastBin_Q1 - preShift_bins) = number of bins needed from Q4
        % flipping list of bins from last to first lets us use absolute value as index
        bins_unique_flipped = flipud(bins_unique);
        first_pre_downshiftBin = bins_unique_flipped( abs(lastBin_Q1 - numPreshiftBins) );
        
        % array of pre-downshift bin numbers in chronological order
        pre_downshiftBins{counter} = [first_pre_downshiftBin:lastBin_ofPeriod,1:lastBin_Q1];
        
        clear bins_unique_flipped
        
    else
        
        % if no need to tap into Q4 data...
        pre_downshiftBins{counter} = lastBin_Q1 - (numPreshiftBins-1) : lastBin_Q1;
        
    end
    
    % for upshifts
    % bins are already chronological in number, so the job is easy!
    % pre_upshiftBins{counter} = lastBin_downshift - (numPreshiftBins-1) : lastBin_downshift;
    


    % 16. collect growth rate data into bins and calculate stats
    binned_growthRate{counter} = accumarray(bins,growthRt_noNaNs,[],@(x) {x});
    binned_mean{counter} = accumarray(bins,growthRt_noNaNs,[],@mean);
    clear bins
    
   
    
    
    % 17. plot response in growth rate for all timescales over time
    timePerBin_min = timePerBin/60; % time per bin in minutes
    
    if timescale == 900
        color = rgb(color900);
    elseif timescale == 3600
        color = rgb(color3600);
    end
    
    
    
    % plot!
    if strcmp(shiftType,'downshift') == 1
        
        figure(1)  % binned downshift 
       
        %concatenate pre and post downshift data
        % time (same as upshift, just different variable names)
        preDownshift_times = ((numPreshiftBins-1)*-1:0)*timePerBin_min;
        postDownshift_times = (1:length( binned_mean{counter}( upshiftBins{counter} ) ) )*timePerBin_min;
        downshift_times_frep = [preDownshift_times,postDownshift_times];
        
        % growth rate
        preDownshift_growth = binned_mean{counter}(pre_downshiftBins{counter});
        postDownshift_growth = binned_mean{counter}(downshiftBins{counter});
        downshift_growth_frep = [preDownshift_growth;postDownshift_growth];
        
        plot(downshift_times_frep,downshift_growth_frep,'Color',color,'LineWidth',1,'Marker','.')
        hold on
        title(strcat('response to downshift, binned every (',num2str(timePerBin),') sec'))
        
        

    end
    clear x y ycenter ybegin yend color
    clear conditionData_trim1 conditionData_trim2 growthRates_trim1 growthRates_trim2
    clear growthRt_noNaNs index lastBin_Q1 lastBin_ofPeriod lastBin_downshift firstBin_upshift firstBin_downshift
    clear first_pre_downshiftBin date condition e
     
end

figure(1)
title(strcat('mean of each replicate: response to upshift, binned every (',num2str(timePerBin),') sec'))
xlabel('time (min)')
ylabel(strcat('growth rate: (', specificGrowthRate,')'))
axis([numPreshiftBins*-1*timePerBin_min,80,xmin,xmax])



for t = 2 % loop through the two timescales, 15 min and 60 min
   
    clear replicate_means r
    
    if t == 1
        timescale = 900;
        color = rgb(color900);
    else
        timescale = 3600;
        color = rgb(color3600);
    end
    
    
    % determine which columns from collected B1 data go with current timescale
    columns = find(timescale_vector == timescale);
    
    
    % extract 2d matrix of mean replicate data, with columns replicate and rows being time (period bin)
    for r = 1:length(columns)
        replicate_means(:,r) = binned_mean{columns(r)};
    end
    
    
    % calculate mean and stdev across replicates (rows)
    replicate_means_mean = mean(replicate_means,2);
    if t == 1
        downshift_replicate_mean_15min = replicate_means_mean;
    else
        downshift_replicate_mean_60min = replicate_means_mean;
    end
    replicate_means_std = std(replicate_means,0,2); % w=0 normalizes by N-1, w=1 normalizes by N
    
    
    

    
    % plot
    if strcmp(shiftType,'downshift') == 1
        
        figure(2) % binned downshift 
        
        % growth rate
        preDownshift_growth = replicate_means_mean(pre_downshiftBins{columns(1)});
        postDownshift_growth = replicate_means_mean(downshiftBins{columns(1)});
        downshift_growth_fluc_means = [preDownshift_growth;postDownshift_growth];
        
        %plot(downshift_times(1:length(downshift_growth_fluc_means)),downshift_growth_fluc_means,'Color',color,'LineWidth',1,'Marker','.')
        plot(downshift_times_frep(1:length(downshift_growth_fluc_means)),downshift_growth_fluc_means,'Color',color,'LineWidth',1,'Marker','.')
        hold on
        title(strcat('response to downshift, binned every (',num2str(timePerBin),') sec'))
        
        preDownshift_replicates = replicate_means(pre_downshiftBins{columns(1)},:);
        postDownshift_replicates = replicate_means(downshiftBins{columns(1)},:);
        downshift_replicates = [preDownshift_replicates;postDownshift_replicates];
        
        y = downshift_replicates';
        x = downshift_times_frep(1:length(downshift_replicates))';
        
        figure(3)
        hold on
        s = shadedErrorBar(x,y,{@mean,@std},'lineprops',{'Color',color},'patchSaturation',0.4);
        
        % set face and edge properties
        s.mainLine.LineWidth = 3;
        s.patch.FaceColor = color;
        
    end
    
end

xlabel('time (min)')
ylabel(strcat('growth rate: (', specificGrowthRate,')'))
axis([-5,120,-0.5,3.4])

save('response_flucDownshift_2','downshift_times_frep','downshift_growth_frep')

