%% figure 51
%
%  Goals: growth rate response to upshifts in real time.
%         and plots the mean of all means + standard dev between means


%  General strategy:
%
%         Part A. initialize folder with stored meta data
%         Part B. curves for fluctuating data
%         Part C. curves for single shift data

%  Note: running parts B & C can be reversed in order!
%        current order with C first plots fluc over single shift data



%  last updated: jen, 2019 Jan 4

%  commit: edit such that a save occurs before adding single shift condition,
%          for animation of figures in presentations

% OK let's go!

%% A. initialize

clc
clear

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')


% 0. define shift type, growth rate and time bin of interest
shiftType = 'upshift';
specificGrowthRate = 'log2';
specificColumn = 3; % log2 growth rate
xmin = -0.5;
xmax = 3.5;
timePerBin = 75;

counter = 0; % because timescales will differ between experiments

% color designations
color900 = 'DeepSkyBlue'; %'MediumSeaGreen'; %
color3600 = 'Navy'; % 'MediumPurple'; %
color_single = 'DarkCyan'; % 'DarkBlue'; %

%% C1. curves for single shift data


% 1. create array of experiments of interest, then loop through each
exptArray = [21,22]; % use corresponding dataIndex values for upshift


for e_shift = 1:length(exptArray)
    
    counter = counter + 1;
    
    % 2. initialize experiment meta data
    index = exptArray(e_shift); 
    date = storedMetaData{index}.date;
    
    % define which frames to ignore (noisy tracking)
    if strcmp(date,'2018-06-15') == 1
        ignoredFrames = [112,113,114];
    elseif strcmp(date,'2018-08-01') == 1
        ignoredFrames = [94,95];
    end
    
    
    timescale = storedMetaData{index}.timescale;
    bubbletime = storedMetaData{index}.bubbletime;
    expType = storedMetaData{index}.experimentType;
    shiftTime = storedMetaData{index}.shiftTime;        % sec

    disp(strcat(date, ': analyze!'))
    
    
    
    % 3. load measured data
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-width1p7-jiggle-0p5.mat');
    load(filename,'D5','T')
    
    
    
    
    % 4. specify condition of interest (fluc) and build data matrix
    condition = 1;                                      % 1 = fluctuating
    xy_start = storedMetaData{index}.xys(condition,1);
    xy_end = storedMetaData{index}.xys(condition,end);
    conditionData = buildDM(D5, T, xy_start, xy_end,index,expType);
    
    
    
    % 6. isolate volume (Va), timestamp, mu, drop and curveID data
    volumes = conditionData(:,11);            % col 11 = calculated va_vals (cubic um)
    timestamps_sec = conditionData(:,2);      % col 2  = timestamp in seconds
    isDrop = conditionData(:,4);              % col 4  = isDrop, 1 marks a birth event
    curveFinder = conditionData(:,5);         % col 5  = curve finder (ID of curve in condition)
    trackNum = conditionData(:,20);           % col 20 = total track number in condition
    clear xy_start xy_end
    
    
    % 7. calculate growth rate
    growthRates = calculateGrowthRate(volumes,timestamps_sec,isDrop,curveFinder,trackNum);
    clear curveFinder trackNum isDrop volumes
    
    
    
    % 8. isolate data to stabilized regions of growth
    %    NOTE: errors (excessive negative growth rates) occur at trimming
    %          point if growth rate calculation occurs AFTER time trim.
    
    minTime = 2.5; % for single shift data, unlike original fluc
    
    maxTime = bubbletime(condition);
    timestamps_sec = conditionData(:,2); % time in seconds converted to hours
    timestamps_hr = timestamps_sec / 3600;
    clear condition
    
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
    clear growthRates conditionData maxTime minTime timestamps_sec timestamps_hr

    
     
    % 9. isolate selected specific growth rate
    growthRt = growthRates_trim2(:,specificColumn);
    % specificColumn is already defined in Part B.
    % not re-defining it here ensures that we use the same metric between both
    clear growthRates_trim1 growthRates_trim2
    

    
    % 10. isolate corrected timestamp
    correctedTime = conditionData_trim2(:,22); % col 22 = timestamps corrected for signal lag
    clear D5 T isDrop conditionData_trim1
    
    
    
    
    % 11. assign NaN to all growth rates associated with frames to ignore
    frameNum = conditionData_trim2(:,16); % col 16 = original frame number
    growthRt_ignorant = growthRt;
    for fr = 1:length(ignoredFrames)
        growthRt_ignorant(frameNum == ignoredFrames(fr),1) = NaN;
    end
    clear fr
    
    
    
    % 12. remove nans from data analysis
    growthRt_noNaNs = growthRt_ignorant(~isnan(growthRt_ignorant),:);
    correctedTime_noNans = correctedTime(~isnan(growthRt_ignorant),:);
    clear growthRt growthRt_ignorant correctedTime frameNum
    
    
    
    
    % 13. assign corrected timestamps to bins, by which to accumulate growth data
    bins = ceil(correctedTime_noNans/timePerBin);      % bin 1 = first timePerBin sec of experiment
    bins_unique = (min(bins):max(bins))';              % avoids missing bins due to lack of data
    
    
    % find bins after shift
    % generalized for single shift experiments
    first_postshiftBin_single = ceil(shiftTime/timePerBin) + 1; % shift occurs at shiftTime/timePerBin, so first full bin with shifted data is the next one
    postshiftBins_single{counter} = (first_postshiftBin_single:max(bins))';
    
    
    
    
    % 14. choose which pre-shift data bins to plot
    
    % single shift experiments don't have high/low phase interruptions
    % however, they do have bins missing data!
    numPreshiftBins = 10;
    preShift_bins{counter} = numPreshiftBins;

    % determine pre-shift bins
    index_single_shift = find(bins_unique == first_postshiftBin_single);
    pre_upshiftBins{counter} = bins_unique(index_single_shift-numPreshiftBins : index_single_shift-1); % same bins in both single down and upshifts
    pre_downshiftBins{counter} = bins_unique(index_single_shift-numPreshiftBins : index_single_shift-1);
    

    
    
    % 15. collect growth rate data into bins and calculate stats
    %     WARNING: bins variable does not i
    binned_growthRate{counter} = accumarray(bins,growthRt_noNaNs,[],@(x) {x});
    binned_mean{counter} = accumarray(bins,growthRt_noNaNs,[],@mean);
    clear bins
    

    
    % 16. plot response in growth rate for all timescales over time
    timePerBin_min = timePerBin/60; % time per bin in minutes
    color = rgb(color_single);

    
    
    figure(1)   % binned upshift
    
    % plot!
    preUpshift_times = ((numPreshiftBins-1)*-1:0)*timePerBin_min;
    postUpshift_times = (1:length(binned_mean{counter}(postshiftBins_single{counter})))*timePerBin_min;
    upshift_times_gapped = [preUpshift_times, postUpshift_times];
    
    preUpshift_growth = binned_mean{counter}(pre_upshiftBins{counter});
    postUpshift_growth_single = binned_mean{counter}(postshiftBins_single{counter});
    upshift_growth_gapped = [preUpshift_growth; postUpshift_growth_single];
    
    
    % don't plot zeros that are place holders for gaps in data
    upshift_growth = upshift_growth_gapped(upshift_growth_gapped > 0);
    upshift_times = upshift_times_gapped(upshift_growth_gapped > 0);
    
    
    plot(upshift_times,upshift_growth,'Color',color,'LineWidth',1,'Marker','.')
    
    
    % store data for shaded error bars
    binned_singles{counter} = upshift_growth_gapped;
    binned_singles_times{counter} = upshift_times_gapped;
    
    
    grid on
    hold on
    title(strcat('response to upshift, binned every (',num2str(timePerBin),') sec'))
    

     
end


xlabel('time (min)')
ylabel(strcat('growth rate: (', specificGrowthRate,')'))
axis([numPreshiftBins*-1*timePerBin_min,160,xmin,xmax])

clear e_shift date
clear replicate_means replicate_means_mean replicate_means_std

%% C2. find mean and standard deviation between single-shift replicates


% extract 2d matrix of mean replicate data, with columns replicate and rows being time (period bin)
for rep = 1:2 % two replicates for upshift and downshift as of 2018-12-04
    col = rep;
    
    curr_rep_means = binned_singles{col}(1:149);
    curr_rep_means(curr_rep_means == 0) = NaN; % make zeros NaN prior to averaging between replicates
    replicate_single_means(rep,:) = curr_rep_means;
    replicate_single_times(rep,:) = binned_singles_times{col}(1:149); % demonstrates indeces are of the same timeline
    
end
clear rep    


% calculate mean and stdev across replicates (rows)
replicate_single_means_mean = mean(replicate_single_means);
replicate_single_means_std = std(replicate_single_means);



% don't plot zeros that are place holders for gaps in data
upshift_means = replicate_single_means_mean(replicate_single_means_mean > 0);
upshift_stds = replicate_single_means_std(replicate_single_means_mean > 0);
upshift_times = upshift_times_gapped(replicate_single_means_mean > 0);

upshift_replicate_means = replicate_single_means(:,replicate_single_means_mean > 0);



% plot
figure(2)
hold on
plot(upshift_times,upshift_means,'Color',color,'LineWidth',1,'Marker','.')
title('reponse to upshifts: mean of replicate means')
ylabel('growth rate: log2 (1/hr)')
xlabel('time (min)')


figure(3)
hold on
ss = shadedErrorBar(upshift_times,upshift_replicate_means,{@nanmean,@nanstd},'lineprops',{'Color',color},'patchSaturation',0.3);
title('response to upshift, mean and std of replicate means')

% set face and edge properties
ss.mainLine.LineWidth = 3;
ss.patch.FaceColor = color;
axis([-5,120,0,3.1])

% save only single shift data
figure(3)
plotName = strcat('figure51-upshift-',specificGrowthRate,'-mean&std-singleShiftOnly');
saveas(gcf,plotName,'epsc')

%% B1. compile and plot up or downshift data for each replicate 


% 1. create array of experiments of interest, then loop through each
exptArray = 9:15;


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
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    
    if strcmp(date,'2017-11-12') == 1
        filename = 'lb-fluc-2017-11-12-width1p4-jiggle-0p5.mat';
    elseif ischar(timescale) == 0
        filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
    end
    load(filename,'D5','T')
    clear experimentFolder
    
    
    
    % 4. specify condition of interest (fluc) and build data matrix
    condition = 1;                                      % 1 = fluctuating
    xy_start = storedMetaData{index}.xys(condition,1);
    xy_end = storedMetaData{index}.xys(condition,end);
    conditionData = buildDM(D5, T, xy_start, xy_end,index,expType);
    
    
    
    % 5. isolate volume (Va), timestamp, mu, drop and curveID data
    volumes = conditionData(:,11);            % col 11 = calculated va_vals (cubic um)
    timestamps_sec = conditionData(:,2);      % col 2  = timestamp in seconds
    isDrop = conditionData(:,4);              % col 4  = isDrop, 1 marks a birth event
    curveFinder = conditionData(:,5);         % col 5  = curve finder (ID of curve in condition)
    trackNum = conditionData(:,20);           % col 20 = total track number in condition
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
    upshift_times = [preUpshift_times,postUpshift_times];
    
    % growth rate
    preUpshift_growth = binned_mean{counter}(pre_upshiftBins{counter});
    postUpshift_growth = binned_mean{counter}(upshiftBins{counter});
    upshift_growth = [preUpshift_growth;postUpshift_growth];
    
    
    figure(1)   % mean of each replicate
    plot(upshift_times,upshift_growth,'Color',color,'LineWidth',1,'Marker','.')
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

%% B2. find mean and standard deviation between replicates
        
for t = 1:2 % loop through the two timescales, 15 min and 60 min
   
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
    replicate_means_std = std(replicate_means,0,2); % w=0 normalizes by N-1, w=1 normalizes by N
    


    
    % plot
    clear preUpshift_growth postUpshift_growth
    
    % concatenate growth rate to match time relative to shift
    preUpshift_growth = replicate_means_mean(pre_upshiftBins{columns(1)});
    postUpshift_growth = replicate_means_mean(upshiftBins{columns(1)});
    upshift_growth = [preUpshift_growth;postUpshift_growth];
    
    
    figure(2) % mean of replicate means
    plot(upshift_times(1:length(upshift_growth)),upshift_growth,'Color',color,'LineWidth',1,'Marker','.')
    hold on
    
    
    % prepare matrix of replicate data for errorbar plotting
    preUpshift_replicates = replicate_means(pre_upshiftBins{columns(1)},:);
    postUpshift_replicates = replicate_means(upshiftBins{columns(1)},:);
    upshift_replicates = [preUpshift_replicates;postUpshift_replicates];
    
    y = upshift_replicates';
    x = upshift_times(1:length(upshift_replicates))';
    
    
    figure(3) % mean and std of replicate means
    hold on
    s = shadedErrorBar(x,y,{@mean,@std},'lineprops',{'Color',color},'patchSaturation',0.4);
    
    % set face and edge properties
    s.mainLine.LineWidth = 3;
    s.patch.FaceColor = color;

end

xlabel('time (min)')
ylabel(strcat('growth rate: (', specificGrowthRate,')'))
%axis([numPreshiftBins*-1*timePerBin_min,80,0,3.1])
axis([-5,120,0,3.1])


clear r t color x y timescale bubbletime 

%% D1. save and close complete figures

cd('/Users/jen/Documents/StockerLab/Data_analysis/currentPlots/')


figure(1)
plotName = strcat('figure51-upshift-',specificGrowthRate,'-individualReplicates');
saveas(gcf,plotName,'epsc')
close(gcf)

figure(2)
plotName = strcat('figure51-upshift-',specificGrowthRate,'-meanOfReplicates');
saveas(gcf,plotName,'epsc')
close(gcf)

figure(3)
%plotName = strcat('figure51-upshift-',specificGrowthRate,'-mean&std');
plotName = strcat('figure51-upshift-',specificGrowthRate,'-mean&std-flucOnly');
%plotName = strcat('figure51-upshift-',specificGrowthRate,'-mean&std-flucAfterSingleShift');
saveas(gcf,plotName,'epsc')
close(gcf)



