%% figure 53
%
%  Goals: growth rate response to upshifts in real time.
%         plots growth rate over time, with t = 0 as time of shift

%         aggregates data from all replicates and plots raw and filtered mean



%  General strategy:
%
%         Part A. initialize folder with stored meta data
%         Part B. curves for fluctuating data
%         Part C. curve for single shift data



%  last updated: jen, 2018 Oct 23

%  commit: first commit, aggregates data from all replicates. plots raw and sgolay filtered mean


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

specificGrowthRate = 'log2';
specificColumn = 3; % log2 growth rate

prompt = 'Enter time per bin in seconds as double (i.e. 25): ';
timePerBin = input(prompt);

clear prompt




%% Part B. overlay up or downshift data

%  strategy:
%   
%   - loop through experiments to pull out:
%         i. instantaneous growth rates
%        ii. signal corrected timestamps
%
%   - within loop, concatenate data belonging to the same timescale
%   - outside of loop, once all data is compiled:
%         i. bin data by time pre- or post-shift
%        ii. smooth with filter
%       iii. plot



% 0. initialize arrays for growth rate concatenation
growthCat900 = [];
timeCat900 = [];

growthCat3600 = [];
timeCat3600 = [];


% 1. create array of experiments of interest, then loop through each
exptArray = 9:15;

%counter = 0; % because timescales will differ between experiments
for e = 1:length(exptArray)
    
    %counter = counter + 1;
    
    % 2. initialize experiment meta data
    index = exptArray(e);                               % previous, dataIndex(e);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
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
    clear filename
    
    
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
    clear growthRates conditionData maxTime minTime timestamps_sec condition bubbletime
    
    
    
    % 8. isolate selected specific growth rate
    if strcmp(specificGrowthRate,'raw') == 1
        specificColumn = 1;         % for selecting appropriate column in growthRates
        xmin = -5;                  % lower limit for plotting x axis
        xmax = 25;                  % upper limit for plotting x axis
    elseif strcmp(specificGrowthRate,'norm') == 1
        specificColumn = 2;
        xmin = -1;
        xmax = 3;
    elseif strcmp(specificGrowthRate,'log2') == 1
        specificColumn = 3;
        xmin = -0.5;
        xmax = 3.5;
    elseif strcmp(specificGrowthRate,'lognorm') == 1
        specificColumn = 4;
        xmin = -0.5;
        xmax = 1;
    end
    
    growthRt = growthRates_trim2(:,specificColumn);
    clear timestamps_hr growthRates_trim1 growthRates_trim2
    
    
    
    
    % 9. isolate corrected timestamp
    if strcmp(date, '2017-10-10') == 1
        correctedTime = conditionData_trim2(:,2);
    else
        correctedTime = conditionData_trim2(:,22); % col 22 = timestamps corrected for signal lag
    end
    clear D5 T isDrop conditionData_trim1 conditionData_trim2
    
    
    
    
    % 10. remove nans from data analysis
    growthRt_noNaNs = growthRt(~isnan(growthRt),:);
    correctedTime_noNans = correctedTime(~isnan(growthRt),:);
    clear growthRt correctedTime times_trim1
    
    
    
    % 11. concatenate growth rate and corrected time data base on timescale
    if timescale == 900
        
        growthCat900 = [growthCat900; growthRt_noNaNs];
        timeCat900 = [timeCat900; correctedTime_noNans];
        
    elseif timescale == 3600
        
        growthCat3600 =[growthCat3600; growthRt_noNaNs];
        timeCat3600 = [timeCat3600; correctedTime_noNans];
        
    end
    clear growthRt_noNaNs correctedTime_noNans
    
end
clear date e index


%% loop through concatenated timescale data to plot mean across replicates

% 12. loop through concatenated timescale data to plot mean across replicates

% initialize smoothing parameters
order = 1;        % Use sgolay to smooth the signal, with 5-sample frames and 2nd-order polynomials.
framelength = 3;
b = sgolay(order,framelength);

for t = 1:2
    
    if t == 1
        timescale = 900;
        color = rgb('Aquamarine');
        growthCat = growthCat900;
        timeCat = timeCat900;
        
    else
        timescale = 3600;
        color = rgb('Indigo');
        growthCat = growthCat3600;
        timeCat = timeCat3600;

    end
    
    
    % 13. assign corrected timestamps to bins, by which to accumulate growth data
    timeInPeriods = timeCat./timescale; % unit = sec/sec
    timeInPeriodFraction = timeInPeriods - floor(timeInPeriods);
    
    timeInPeriodFraction_inSeconds = timeInPeriodFraction * timescale;
    bins = ceil(timeInPeriodFraction_inSeconds/timePerBin);
    bins_unique = unique(bins);
    
    
    
    
    % 14. find which bins are boundaries signal phases
    lastBin_Q1 = (timescale/timePerBin)/4;                      % last bin before downshift
    firstBin_downshift = (timescale/4)/timePerBin + 1;
    lastBin_downshift = (timescale*3/4)/timePerBin;             % last bin before upshift
    
    firstBin_upshift = (timescale*3/4)/timePerBin + 1;          % first bin of upshift
    lastBin_ofPeriod = timescale/timePerBin;                    % total bins in signal period
    
    
    
    
    % 15. list bins chronologically to combine broken up high nutrient phase
    %       i.e. start of upshift is Q4, concatenate Q1
    downshiftBins = firstBin_downshift:lastBin_downshift;
    upshiftBins = [firstBin_upshift:lastBin_ofPeriod, 1:lastBin_Q1];
    clear timeInPeriods timeInPeriodFraction timeInPeriodFraction_inSeconds  
    
    
    
    % 16. choose which pre-shift data bins to plot
    % upshift used here to determine how many bins, but it is equal in length to downshift
    if length(upshiftBins) >= 5 
        numPreshiftBins = 4;
    else
        numPreshiftBins = 2;
    end
    
    
    % for downshifts
    % shorter timescales (less bins) require pulling from Q4 growth data
    if lastBin_Q1 - numPreshiftBins <= 0
        
        % absolute value of (lastBin_Q1 - preShift_bins) = number of bins needed from Q4
        % flipping list of bins from last to first lets us use absolute value as index
        bins_unique_flipped = flipud(bins_unique);
        first_pre_downshiftBin = bins_unique_flipped( abs(lastBin_Q1 - numPreshiftBins) ); % 1st for plotting, not in general
        
        % array of pre-downshift bin numbers in chronological order
        pre_downshiftBins = [first_pre_downshiftBin:lastBin_ofPeriod,1:lastBin_Q1];
        clear bins_unique_flipped first_pre_downshiftBin
        
    else
        
        % if no need to tap into Q4 data...
        pre_downshiftBins = lastBin_Q1 - (numPreshiftBins-1) : lastBin_Q1;
        
    end
    
    % for upshifts
    % bins are already chronological in number, so the job is easy!
    pre_upshiftBins = lastBin_downshift - (numPreshiftBins-1) : lastBin_downshift;


    
    % 16. collect growth rate data into bins and calculate stats
    binned_growthRate{t} = accumarray(bins,growthCat,[],@(x) {x});
    binned_mean{t} = accumarray(bins,growthCat,[],@mean);
    clear bins bins_unique
    
    
    
    % 17. plot response in growth rate for all timescales over time
    timePerBin_min = timePerBin/60; % time per bin in minutes
    

    % plot!
    if strcmp(shiftType,'upshift') == 1
        
        %concatenate pre and post upshift data
        
        % time
        preUpshift_times = ((numPreshiftBins-1)*-1:0)*timePerBin_min;
        postUpshift_times = (1:length( binned_mean{t}( upshiftBins ) ))*timePerBin_min;
        upshift_times = [preUpshift_times,postUpshift_times];
        
        % growth rate
        preUpshift_growth = binned_mean{t}(pre_upshiftBins);
        postUpshift_growth = binned_mean{t}(upshiftBins);
        upshift_growth = [preUpshift_growth;postUpshift_growth];
        
        
        figure(1)   % binned upshift
        plot(upshift_times,upshift_growth,'Color',color,'LineWidth',1)
        hold on
        title(strcat('response to upshift, binned every (',num2str(timePerBin),') sec'))
        
        
        % upshift, smoothed
        % pre shift growth signal
        x_pre = preUpshift_growth;
        ycenter_pre = conv(x_pre,b((framelength+1)/2,:),'valid');
        ybegin_pre = b(end:-1:(framelength+3)/2,:) * x_pre(framelength:-1:1);
        yend_pre = b((framelength-1)/2:-1:1,:) * x_pre(end:-1:end-(framelength-1));
        y_pre = [ybegin_pre; ycenter_pre; yend_pre];
       
        % post shift growth signal
        x_post = postUpshift_growth;
        ycenter_post = conv(x_post,b((framelength+1)/2,:),'valid');
        ybegin_post = b(end:-1:(framelength+3)/2,:) * x_post(framelength:-1:1);
        yend_post = b((framelength-1)/2:-1:1,:) * x_post(end:-1:end-(framelength-1));
        y_post = [ybegin_post; ycenter_post; yend_post];
        
        y = [y_pre;y_post];
        plot(upshift_times,y,'Color',color,'LineWidth',2)
        hold on
        title(strcat('response to downshift, binned every (',num2str(timePerBin),') sec - filtered'))
        
        clear upshift_times upshift_growth
        
    else % downshift
        
        figure(3)  % binned downshift
        
        %concatenate pre and post downshift data
        % time (same as upshift, just different variable names)
        preDownshift_times = ((numPreshiftBins-1)*-1:0)*timePerBin_min;
        postDownshift_times = (1:length( binned_mean{t}( upshiftBins ) ))*timePerBin_min;
        downshift_times = [preDownshift_times,postDownshift_times];
        
        % growth rate
        preDownshift_growth = binned_mean{t}(pre_downshiftBins);
        postDownshift_growth = binned_mean{t}(downshiftBins);
        downshift_growth = [preDownshift_growth;postDownshift_growth];
        
        plot(downshift_times,downshift_growth,'Color',color,'LineWidth',1)
        hold on
        title(strcat('response to downshift, binned every (',num2str(timePerBin),') sec'))

        
        % downshift, smoothed
        % pre shift growth signal
        x_pre = preDownshift_growth;
        ycenter_pre = conv(x_pre,b((framelength+1)/2,:),'valid');
        ybegin_pre = b(end:-1:(framelength+3)/2,:) * x_pre(framelength:-1:1);
        yend_pre = b((framelength-1)/2:-1:1,:) * x_pre(end:-1:end-(framelength-1));
        y_pre = [ybegin_pre; ycenter_pre; yend_pre];
       
        % post shift growth signal
        x_post = postDownshift_growth;
        ycenter_post = conv(x_post,b((framelength+1)/2,:),'valid');
        ybegin_post = b(end:-1:(framelength+3)/2,:) * x_post(framelength:-1:1);
        yend_post = b((framelength-1)/2:-1:1,:) * x_post(end:-1:end-(framelength-1));
        y_post = [ybegin_post; ycenter_post; yend_post];
        
        y = [y_pre;y_post];
        plot(downshift_times,y,'Color',color,'LineWidth',2)
        hold on
        title(strcat('response to downshift, binned every (',num2str(timePerBin),') sec - filtered'))
        
        clear downshift_times downshift_growth
        
    end
    clear preUpshift_growth postUpshift_growth preUpshift_times postUpshift_times
    
    x = [x_pre;x_post];
    
    figure(5)
    plot([x y])
    hold on
    
end
xlabel('bin number')
ylabel('growth rate')
title('+/- sgolay filter')

clear x y x_pre x_post y_post y_pre timescale 
clear ycenter_pre ybegin_pre yend_pre ycenter_post ybegin_post yend_post color
clear lastBin_Q1 lastBin_ofPeriod lastBin_downshift firstBin_upshift firstBin_downshift
clear growthCat timeCat
clear pre_downshiftBins pre_upshiftBins
clear numPreshiftBins order framelength b t
clear downshiftBins upshiftBins


%% Part C. curves for single shift data


% 1. create array of experiments of interest, then loop through each
if strcmp(shiftType,'upshift')
    exptArray = [21,22]; % use corresponding dataIndex values
else
    exptArray = [25,26,27];
end


growthCat_singleUP = [];
binCat_singleUP = [];

growthCat_singleDown = [];
binCat_singleDown = [];


for e_shift = 1:length(exptArray)
    
    
    % 2. initialize experiment meta data
    index = exptArray(e_shift);                               % previous, dataIndex(e);
    date = storedMetaData{index}.date;
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
    
    
    
    % 5. isolate volume (Va), timestamp, mu, drop and curveID data
    volumes = conditionData(:,11);            % col 11 = calculated va_vals (cubic um)
    timestamps_sec = conditionData(:,2);      % col 2  = timestamp in seconds
    isDrop = conditionData(:,4);              % col 4  = isDrop, 1 marks a birth event
    curveFinder = conditionData(:,5);         % col 5  = curve finder (ID of curve in condition)
    trackNum = conditionData(:,20);           % col 20 = total track number in condition
    clear xy_start xy_end
    
    
    % 6. calculate growth rate
    growthRates = calculateGrowthRate(volumes,timestamps_sec,isDrop,curveFinder,trackNum);
    clear curveFinder trackNum isDrop volumes
    
    
    
    % 7. isolate data to stabilized regions of growth
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
    clear growthRates conditionData maxTime minTime timestamps_sec timestamps_hr times_trim1
    
    
    
    % 8. isolate selected specific growth rate
    growthRt = growthRates_trim2(:,specificColumn);
    % specificColumn is already defined in Part B.
    % not re-defining it here ensures that we use the same metric between both
    
    
    
    
    % 9. isolate corrected timestamp
    correctedTime = conditionData_trim2(:,22); % col 22 = timestamps corrected for signal lag
    clear D5 T isDrop conditionData_trim1 conditionData_trim2
    
    
    
    
    % 10. remove nans from data analysis
    growthRt_noNaNs = growthRt(~isnan(growthRt),:);
    correctedTime_noNans = correctedTime(~isnan(growthRt),:);
    clear growthRt correctedTime growthRates_trim1 growthRates_trim2
    
    
    
    % 11. assign corrected timestamps to bins, by which to accumulate growth data
    bins = ceil(correctedTime_noNans/timePerBin);      % bin 1 = first timePerBin sec of experiment
    bins_unique = (min(bins):max(bins))';              % avoids missing bins due to lack of data
    
    % find bins after shift
    % generalized for single shift experiments
    first_postshiftBin_single = ceil(shiftTime/timePerBin) + 1; % shift occurs at shiftTime/timePerBin, so first full bin with shifted data is the next one
    
    
    % 14. choose which pre-shift data bins to plot
    % single shift experiments don't have high/low phase interruptions
    % however, they do have bins missing data!
    numPreshiftBins = 10;
    
    % normalize bins such that 1:numPreshiftBins are the bins immediately before shift
    bins_normalizedByShiftTime = bins - (first_postshiftBin_single-1) + numPreshiftBins;
    
    % remove data earlier preshift range of interest
    bins_normalizedByShiftTime_trimmed = bins_normalizedByShiftTime(bins_normalizedByShiftTime > 0);
    growthRt_normalizedByShiftTime_trimmed = growthRt_noNaNs(bins_normalizedByShiftTime > 0);
    
    % define pre-shift bins
    index_single_shift = numPreshiftBins+1;
    pre_upshiftBins = 1:numPreshiftBins;
    pre_downshiftBins = 1:numPreshiftBins;
    
    

    % 11. concatenate growth rate and corrected time data based on shift type
    if strcmp(expType,'upshift') == 1
        
        growthCat_singleUP = [growthCat_singleUP; growthRt_normalizedByShiftTime_trimmed];
        binCat_singleUP = [binCat_singleUP; bins_normalizedByShiftTime_trimmed];
        
    elseif strcmp(expType,'downshift') == 1
        
        growthCat_singleDown =[growthCat_singleDown; growthRt_normalizedByShiftTime_trimmed];
        binCat_singleDown = [binCat_singleDown; bins_normalizedByShiftTime_trimmed];
        
    end
    clear growthRt_noNaNs correctedTime_noNans
    
    
end
clear date e index specificColumn timescale bubbletime

%% plot mean across single shift replicates

t = 3;
if strcmp(shiftType,'upshift') == 1
    
    growthCat = growthCat_singleUP;
    binCat = binCat_singleUP;
    
else
    growthCat = growthCat_singleDown;
    binCat = binCat_singleDown;
    
end

    
% 15. collect growth rate data into bins and calculate stats
%     WARNING: bins variable does not i
binned_growthRate{t} = accumarray(binCat,growthCat,[],@(x) {x});
binned_mean{t} = accumarray(binCat,growthCat,[],@mean);
clear bins



% 17. plot response in growth rate for all timescales over time
timePerBin_min = timePerBin/60; % time per bin in minutes
color = rgb('DarkOliveGreen');


if strcmp(shiftType,'upshift') == 1
    
    figure(1)   % binned upshift
    
    % plot!
    preUpshift_times = flip( (pre_upshiftBins*timePerBin_min)*-1 );
    postUpshift_times = ((index_single_shift: max(binCat))-index_single_shift)*timePerBin_min;
    upshift_times_gapped = [preUpshift_times, postUpshift_times];
    
    preUpshift_growth = binned_mean{t}([pre_upshiftBins,11,12]);
    postUpshift_growth_single = binned_mean{t}((index_single_shift+2: max(binCat)));
    upshift_growth_gapped = [preUpshift_growth; postUpshift_growth_single];
    
    
    % don't plot zeros that are place holders for gaps in data
    upshift_growth = upshift_growth_gapped(upshift_growth_gapped > 0);
    upshift_times = upshift_times_gapped(upshift_growth_gapped > 0);
    
    
    plot(upshift_times,upshift_growth,'Color',color,'LineWidth',1)
    hold on
    title(strcat('response to upshift, binned every (',num2str(timePerBin),') sec'))
    
    
    % upshift, smoothed
    
    order = 2;        % Use sgolay to smooth the signal
    framelength = 9;
    b = sgolay(order,framelength);
    
    x_pre = preUpshift_growth(preUpshift_growth > 0);
    ycenter_pre = conv(x_pre,b((framelength+1)/2,:),'valid');
    ybegin_pre = b(end:-1:(framelength+3)/2,:) * x_pre(framelength:-1:1);
    yend_pre = b((framelength-1)/2:-1:1,:) * x_pre(end:-1:end-(framelength-1));
    y_pre = [ybegin_pre; ycenter_pre; yend_pre];
    
    order = 3;        % Use sgolay to smooth the signal, with 35-sample frames and 2nd-order polynomials.
    framelength = 19;
    b = sgolay(order,framelength);
    
    x_post = postUpshift_growth_single(postUpshift_growth_single > 0);
    ycenter_post = conv(x_post,b((framelength+1)/2,:),'valid');
    ybegin_post = b(end:-1:(framelength+3)/2,:) * x_post(framelength:-1:1);
    yend_post = b((framelength-1)/2:-1:1,:) * x_post(end:-1:end-(framelength-1));
    y_post = [ybegin_post; ycenter_post; yend_post];
    
    y = [y_pre;y_post];
    plot(upshift_times,y,'Color',color,'LineWidth',2)
    hold on
   
else
    
    figure(3)    % downshift
    
    % plot!
    preDownshift_times = flip( (pre_downshiftBins*timePerBin_min)*-1 );
    postDownshift_times = ((index_single_shift: max(binCat))-index_single_shift)*timePerBin_min;
    downshift_times_gapped = [preDownshift_times, postDownshift_times];
    
    preDownshift_growth_gapped = binned_mean{t}(pre_downshiftBins);
    postDownshift_growth_single = binned_mean{t}((index_single_shift: max(binCat)));
    downshift_growth_gapped = [preDownshift_growth_gapped; postDownshift_growth_single];
    
    % don't plot zeros that are place holders for gaps in data
    downshift_growth = downshift_growth_gapped(downshift_growth_gapped > 0);
    downshift_times = downshift_times_gapped(downshift_growth_gapped > 0);
    
    plot(downshift_times,downshift_growth,'Color',color,'LineWidth',1)
    grid on
    hold on
    title(strcat('response to downshift, binned every (',num2str(timePerBin),') sec'))
    
    
    
    
    % downshift, smoothed
    preDownshift_growth = preDownshift_growth_gapped(preDownshift_growth_gapped > 0);
    postDownshift_growth = postDownshift_growth_single(postDownshift_growth_single > 0);
    
    order = 2;        % Use sgolay to smooth the signal
    framelength = 9;
    b = sgolay(order,framelength);
    
    x_pre = preDownshift_growth;
    ycenter_pre = conv(x_pre,b((framelength+1)/2,:),'valid');
    ybegin_pre = b(end:-1:(framelength+3)/2,:) * x_pre(framelength:-1:1);
    yend_pre = b((framelength-1)/2:-1:1,:) * x_pre(end:-1:end-(framelength-1));
    y_pre = [ybegin_pre; ycenter_pre; yend_pre];
    
    
    order = 3;        % Use sgolay to smooth the signal, with 35-sample frames and 2nd-order polynomials.
    framelength = 29;
    b = sgolay(order,framelength);
    
    x_post = postDownshift_growth;
    ycenter_post = conv(x_post,b((framelength+1)/2,:),'valid');
    ybegin_post = b(end:-1:(framelength+3)/2,:) * x_post(framelength:-1:1);
    yend_post = b((framelength-1)/2:-1:1,:) * x_post(end:-1:end-(framelength-1));
    y_post = [ybegin_post; ycenter_post; yend_post];
    
    
    y = [y_pre;y_post];
    plot(downshift_times,y,'Color',color,'LineWidth',2)
    hold on
    
    
end

x = [x_pre; x_post];
figure(6)
plot([x y])
xlabel('bin number')
ylabel(strcat('growth rate: (', specificGrowthRate,')'))
title('+/- sgolay filter')



%% test smoothing with sgolay filter
%data_to_smooth = binned_mean{counter}(single_shiftBins_unique{counter});

% Use sgolay to smooth the signal. Use 21-sample frames and 4th-order polynomials.
order = 2;
framelength = 27;

b = sgolay(order,framelength);


% Compute the steady-state portion of the signal by convolving it with the center row of b.
ycenter = conv(x,b((framelength+1)/2,:),'valid');


% Compute the transients.
% Use the last rows of b for the startup and the first rows of b for the terminal.
ybegin = b(end:-1:(framelength+3)/2,:) * x(framelength:-1:1);
yend = b((framelength-1)/2:-1:1,:) * x(end:-1:end-(framelength-1));


% Concatenate the transients and the steady-state portion to generate the complete smoothed signal.
% Plot the original signal and the Savitzky-Golay estimate.
figure(2)
y = [ybegin; ycenter; yend];
plot([x y])
%

%%


% plot average and standard dev of experimental means

% define color scheme
color300 = rgb('Chocolate');
color900 = rgb('ForestGreen');
color3600 = rgb('Amethyst');
color_single = rgb('MidnightBlue');


% brute force
% isolate experiments based on type (as indicated in notes)
binnedMean_300 = [binned_mean{1}, binned_mean{2}, binned_mean{3}];
binnedMean_900 = [binned_mean{4}, binned_mean{5}, binned_mean{6}];
binnedMean_3600 = [binned_mean{7}, binned_mean{8}, binned_mean{9}];

binnedMean_shift = [binned_mean{10}, binned_mean{11}(1:length(binned_mean{10}))];
binnedMean_shift(binnedMean_shift == 0) = NaN;


% collect mean from each row
m300 = mean(binnedMean_300,2)/log(2);
m900 = mean(binnedMean_900,2)/log(2);
m3600 = mean(binnedMean_3600,2)/log(2);
m_shift = mean(binnedMean_shift,2)/log(2);

% collect standard dev from each row
sd300 = std(binnedMean_300,0,2)/log(2);
sd900 = std(binnedMean_900,0,2)/log(2);
sd3600 = std(binnedMean_3600,0,2)/log(2);
sd_shift = std(binnedMean_shift,0,2)/log(2);


%% upshift response, average between experimental replicates

figure(3)

% post single upshift
tsingle = (-10:length(m_shift(single_shiftBins_unique{10})))*timePerBin_min;
m_single = [m_shift(pre_upshiftBins{10}); m_shift(single_shiftBins_unique{10})];
sd_single = [sd_shift(pre_upshiftBins{10}); sd_shift(single_shiftBins_unique{10})];
errorbar(tsingle,m_single,sd_single,'Color',color_single,'LineWidth',1)
hold on

% 3600
% concatenate pre-upshift and post shift
t3600 = (-4:length(m3600(upshiftBins{9})))*timePerBin_min;
m3600_2plot = [m3600(pre_upshiftBins{9}); m3600(upshiftBins{9})];
sd3600_2plot = [sd3600(pre_upshiftBins{9}); sd3600(upshiftBins{9})];
errorbar(t3600,m3600_2plot,sd3600_2plot,'Color',color3600,'LineWidth',1)
hold on

% 900
% concatenate pre-upshift and post shift
t900 = (-4:length(m900(upshiftBins{6})))*timePerBin_min;
m900_2plot = [m900(pre_upshiftBins{6}); m900(upshiftBins{6})];
sd900_2plot = [sd900(pre_upshiftBins{6}); sd900(upshiftBins{6})];
errorbar(t900,m900_2plot,sd900_2plot,'Color',color900,'LineWidth',1)
hold on

% 300
% concatenate pre-upshift and post shift
t300 = (-2:length(m300(upshiftBins{3})))*timePerBin_min;
m300_2plot = [m300(pre_upshiftBins{3}); m300(upshiftBins{3})];
sd300_2plot = [sd300(pre_upshiftBins{3}); sd300(upshiftBins{3})];
errorbar(t300,m300_2plot,sd300_2plot,'Color',color300,'LineWidth',1)

hold on

title(strcat('response to upshift, binned every (',num2str(timePerBin),') sec'))
xlabel('time (sec)')
ylabel(strcat('growth rate: (', specificGrowthRate ,'/ln(2))'))
axis([numPreshiftBins*-1*timePerBin_min,60,0,4])




%% downshift response, average between experimental replicates
 
figure(4)


% post single upshift
tsingle = (-10:length(m_shift(single_shiftBins_unique{10})))*timePerBin_min;
m_single = [m_shift(pre_downshiftBins{10}); m_shift(single_shiftBins_unique{10})];
sd_single = [sd_shift(pre_downshiftBins{10}); sd_shift(single_shiftBins_unique{10})];
errorbar(tsingle,m_single,sd_single,'Color',color_single,'LineWidth',1)
hold on



% 3600
% concatenate pre-upshift and post shift
t3600 = (-4:length(m3600(downshiftBins{9})))*timePerBin_min;
m3600_2plot = [m3600(pre_downshiftBins{9}); m3600(downshiftBins{9})];
sd3600_2plot = [sd3600(pre_downshiftBins{9}); sd3600(downshiftBins{9})];
errorbar(t3600,m3600_2plot,sd3600_2plot,'Color',color3600,'LineWidth',1)
hold on


% 900
% concatenate pre-upshift and post shift
t900 = (-4:length(m900(downshiftBins{6})))*timePerBin_min;
m900_2plot = [m900(pre_downshiftBins{6}); m900(downshiftBins{6})];
sd900_2plot = [sd900(pre_downshiftBins{6}); sd900(downshiftBins{6})];
errorbar(t900,m900_2plot,sd900_2plot,'Color',color900,'LineWidth',1)
hold on

% 300
% concatenate pre-upshift and post shift
t300 = (-2:length(m300(downshiftBins{3})))*timePerBin_min;
m300_2plot = [m300(pre_downshiftBins{3}); m300(downshiftBins{3})];
sd300_2plot = [sd300(pre_downshiftBins{3}); sd300(downshiftBins{3})];
errorbar(t300,m300_2plot,sd300_2plot,'Color',color300,'LineWidth',1)
hold on

title(strcat('response to downshift, binned every (',num2str(timePerBin),') sec'))
xlabel('time (sec)')
ylabel(strcat('growth rate: (', specificGrowthRate ,')'))
axis([numPreshiftBins*-1*timePerBin_min,60,-1,3.5])





