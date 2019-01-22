%% figure 55
%
%  Goals: from downshift data, calculate
%                   (1) time to stabilization and
%                   (2) mean growth rate during stable region


%  Strategy:
%
%         1. find window of stabilized region that gives fit slope closest to zero
%         2. time until window of zero-most slope = time to stabilization 
%         3. mean growth rate within window = mean growth rate of stabilized signal


%  Note: or single downshift experiments (class 3 below), only one replicate reaches G_low
%        thus, (a) calculate "time to stable" as "time to G_low" 
%              (b) confirm that between "time to G_low" and end that mean gr is G_low
        


%  Summary of code sections:
%
%         Part A. initialize folder with stored meta data
%         Part B. curves for fluctuating data
%         Part C. curve for single shift data
%         Part D. fitting for quantifications
%         Part E. plot quantifications


%  last updated: jen, 2018 Jan 22

%  commit: first commit, quantification of gr responses to downshift
%          

% OK let's go!


%% Part A. initialize

clc
clear

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')


% 0. define shift type, growth rate and time bin of interest
shiftType = 'downshift';

specificGrowthRate = 'log2';
specificColumn = 3; % log2 growth rate
xmin = -0.5;
xmax = 3.5;

timePerBin = 75; % matches binning in shift response plots


% 0. define value of G_low
G_low = 1.07;
std_low = 0.23;


%% Part B. accumulate downshift data


% 1. create array of experiments of interest, then loop through each
exptArray = 9:15;


counter = 0; % because timescales will differ between experiments
for e = 1:length(exptArray)
    
    counter = counter + 1;
    
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
    lastBin_ofPeriod = timescale/timePerBin;                    % total bins in signal period
    
    
    
    
    % 13. list bins chronologically to combine broken up high nutrient phase
    %       i.e. start of upshift is Q4, concatenate Q1
    downshiftBins{counter} = firstBin_downshift:lastBin_downshift;
    clear timeInPeriodFraction timeInPeriodFraction_inSeconds
    
    
    

    % 14. choose which pre-shift data bins to plot
    % decide how much data to plot prior to shift
    if length(downshiftBins{counter}) >= 5
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
        first_pre_downshiftBin = bins_unique_flipped( abs(lastBin_Q1 - numPreshiftBins) );
        
        % array of pre-downshift bin numbers in chronological order
        pre_downshiftBins{counter} = [first_pre_downshiftBin:lastBin_ofPeriod,1:lastBin_Q1];
        
        clear bins_unique_flipped
        
    else
        
        % if no need to tap into Q4 data...
        pre_downshiftBins{counter} = lastBin_Q1 - (numPreshiftBins-1) : lastBin_Q1;
        
    end
    


    % 16. collect growth rate data into bins and calculate stats
    binned_growthRate{counter} = accumarray(bins,growthRt_noNaNs,[],@(x) {x});
    binned_mean{counter} = accumarray(bins,growthRt_noNaNs,[],@mean);
    clear bins
    
   
    
    
    % 17. plot response in growth rate for all timescales over time
    timePerBin_min = timePerBin/60; % time per bin in minutes
    
    if timescale == 900
        color = rgb('Aquamarine');
    elseif timescale == 3600
        color = rgb('Indigo');
    end
    
    
    
    figure(1)  % binned downshift
    
    %concatenate pre and post downshift data
    % time (same as upshift, just different variable names)
    preDownshift_times = ((numPreshiftBins-1)*-1:0)*timePerBin_min;
    postDownshift_times = (1:length( binned_mean{counter}( downshiftBins{counter} ) ) )*timePerBin_min;
    downshift_times = [preDownshift_times,postDownshift_times];
    
    % growth rate
    preDownshift_growth = binned_mean{counter}(pre_downshiftBins{counter});
    postDownshift_growth = binned_mean{counter}(downshiftBins{counter});
    downshift_growth = [preDownshift_growth;postDownshift_growth];
    
    plot(downshift_times,downshift_growth,'Color',color,'LineWidth',1,'Marker','.')
    hold on
    title(strcat('response to downshift, binned every (',num2str(timePerBin),') sec'))
    
    clear x y ycenter ybegin yend color
    clear conditionData_trim1 conditionData_trim2 growthRates_trim1 growthRates_trim2
    clear growthRt_noNaNs index lastBin_Q1 lastBin_ofPeriod lastBin_downshift firstBin_upshift firstBin_downshift
    clear first_pre_downshiftBin date condition e
    
end

xlabel('time (min)')
ylabel(strcat('growth rate: (', specificGrowthRate,')'))
axis([numPreshiftBins*-1*timePerBin_min,80,xmin,xmax])


clear bins_unique bubbletime correctedTime_noNans filename numPreshiftBins
clear timestamps_hr 
clear preDownshift_growth preDownshift_times downshift_growth downshift_times
clear postDownshift_growth postDownshift_times
clear timePerBin_min preShift_bins


%% Part C. accumulate single shift data


% 1. create array of experiments of interest, then loop through each
exptArray = [26,27];  % dataIndex values of single downshift experiments


%counter = 0;  % keep counter value from part B and continue
for e_shift = 1:length(exptArray)
    
    counter = counter + 1;
    
    % 2. initialize experiment meta data
    index = exptArray(e_shift); 
    date = storedMetaData{index}.date;
    
    % define which frames to ignore (noisy tracking)
    if strcmp(date,'2018-08-09') == 1
        ignoredFrames = [115,116,117];
    elseif strcmp(date,'2018-08-08') == 1
        ignoredFrames = [112,113,114];
    end
    
    
    timescale = storedMetaData{index}.timescale;
    bubbletime = storedMetaData{index}.bubbletime;
    expType = storedMetaData{index}.experimentType;
    shiftTime(e_shift) = storedMetaData{index}.shiftTime;        % sec

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
    clear growthRates conditionData maxTime minTime timestamps_sec timestamps_hr

    
     
    % 8. isolate selected specific growth rate
    growthRt = growthRates_trim2(:,specificColumn);
    % specificColumn is already defined in Part B.
    % not re-defining it here ensures that we use the same metric between both
    clear growthRates_trim1 growthRates_trim2
    

    
    % 9. isolate corrected timestamp
    correctedTime = conditionData_trim2(:,22); % col 22 = timestamps corrected for signal lag
    clear D5 T isDrop conditionData_trim1
    
    
    
    % 10. assign NaN to all growth rates associated with frames to ignore
    frameNum = conditionData_trim2(:,16); % col 16 = original frame number
    growthRt_ignorant = growthRt;
    for fr = 1:length(ignoredFrames)
        growthRt_ignorant(frameNum == ignoredFrames(fr),1) = NaN;
    end
    
    
    
    % 11. remove nans from data analysis
    growthRt_noNaNs = growthRt_ignorant(~isnan(growthRt_ignorant),:);
    correctedTime_noNans = correctedTime(~isnan(growthRt_ignorant),:);
    clear growthRt growthRt_ignorant correctedTime frameNum
    
    
    
    
    % 12. assign corrected timestamps to bins, by which to accumulate growth data
    bins = ceil(correctedTime_noNans/timePerBin);      % bin 1 = first timePerBin sec of experiment
    bins_unique = (min(bins):max(bins))';              % avoids missing bins due to lack of data
    
    
    % find bins after shift
    % generalized for single shift experiments
    first_postshiftBin_single = ceil(shiftTime(e_shift)/timePerBin) + 1; % shift occurs at shiftTime/timePerBin, so first full bin with shifted data is the next one
    postshiftBins_single{counter} = (first_postshiftBin_single:max(bins))';
    
    
    
    
    % 13. choose which pre-shift data bins to plot
    
    % single shift experiments don't have high/low phase interruptions
    % however, they do have bins missing data!
    numPreshiftBins = 10;
    %preShift_bins{counter} = numPreshiftBins;

    % determine pre-shift bins
    index_single_shift = find(bins_unique == first_postshiftBin_single);
    %pre_upshiftBins{counter} = bins_unique(index_single_shift-numPreshiftBins : index_single_shift-1); % same bins in both single down and upshifts
    pre_downshiftBins{counter} = bins_unique(index_single_shift-numPreshiftBins : index_single_shift-1);
    

    
    
    % 14. collect growth rate data into bins and calculate stats
    %     WARNING: bins variable does not i
    binned_growthRate{counter} = accumarray(bins,growthRt_noNaNs,[],@(x) {x});
    binned_mean{counter} = accumarray(bins,growthRt_noNaNs,[],@mean);
    clear bins
    
   
    
    
    % 15. plot response in growth rate for all timescales over time
    timePerBin_min = timePerBin/60; % time per bin in minutes
    color = rgb('DarkOliveGreen');
    
    
    figure(1)    % downshift
    
    % plot!
    preDownshift_times = ((numPreshiftBins-1)*-1:0)*timePerBin_min;
    postDownshift_times = (1:length(binned_mean{counter}(postshiftBins_single{counter})))*timePerBin_min;
    downshift_times_gapped = [preDownshift_times, postDownshift_times];
    
    preDownshift_growth_gapped = binned_mean{counter}(pre_downshiftBins{counter});
    postDownshift_growth_single = binned_mean{counter}(postshiftBins_single{counter});
    downshift_growth_gapped = [preDownshift_growth_gapped; postDownshift_growth_single];
    
    % don't plot zeros that are place holders for gaps in data
    downshift_growth = downshift_growth_gapped(downshift_growth_gapped > 0);
    downshift_times = downshift_times_gapped(downshift_growth_gapped > 0);
    
    plot(downshift_times,downshift_growth,'Color',color,'LineWidth',1)
    grid on
    hold on
    title(strcat('response to downshift, binned every (',num2str(timePerBin),') sec'))
    
    
     
end


xlabel('time (min)')
ylabel(strcat('growth rate: (', specificGrowthRate,')'))
axis([numPreshiftBins*-1*timePerBin_min,160,xmin,xmax])


clear bubbletime color conditionData_trim2 date correctedTime_noNans
clear bins_unique e_shift experimentFolder exptArray expType filename
clear fr growthRt_noNaNs ignoredFrames index index_single_shift
clear numPreshiftBins postshiftBins_single preDownshift_times preDownshift_growth
clear shiftType specificColumn specificGrowthRate timePerBin_min
clear timescale times_trim1 downshift_growth downshift_growth_gapped
clear xmax xmin first_postshiftBin_single postDownshift_growth_single
clear postDownshift_times downshift_times_gapped downshift_times


%% Part D. fitting slopes to find time until stable
%       ...and once stable region is determined, calculate mean steady growth rate

clear steadinessTimescale steadinessValue r

expClass = [1,1,1,1,2,2,2,3,3];
cl = unique(expClass);

for currentClass = 1:length(cl)
    
    clear steadiedGR time2steady 
    if currentClass == 1 % 15 min timescale
        
        % isolate 15 min experiment data
        % entire signal
        currentData = cell2mat( binned_mean(expClass==currentClass) );
        finalBin = length(currentData);
        
        % order signal such that 1st half is downshift (start to end of low phase),
        %                        2nd half is upshift (start to end of high phase)
        binOrder = [4:finalBin, 1:3];                   % order for sorting
        orderedData = currentData(binOrder,:);   % sort data so that first low after high starts the array
        
        % generate generic timesignal
        genericTime = ((1:finalBin)*timePerBin)';
        orderedTime = [genericTime, genericTime, genericTime, genericTime];
        
        % from tpt zero to end, find fit a linear function and save slope
        clear slopes steadiedGR section_compiled
        for r = 1:6-2
            
            % select section to compute slope over
            section = orderedData(r:6,:);
            section_time = orderedTime(r:6,:);
            
            for col = 1:4
                section_slope = polyfit(section_time(:,col),section(:,col),1);
                slopes(r,col) = section_slope(1);
                section_compiled{r,col} = section(:,col);
            end
            
        end
        clear r section_slope col
        
        % find which slope is closest to zero
        %slopes_dataOnly = slopes(3:end,:);
        %section_compiled_dataOnly = section_compiled(3:end,:);
        slopes_abs = abs(slopes);
        flats = min(slopes_abs,[],1);
        
        % find times associated with flattest slope
        orderedTime_trim = orderedTime(1:6,:);
        for col = 1:4
            time2steady(col) = orderedTime_trim(slopes_abs(:,col) == flats(col));
            steadiedGR(col) = cellfun(@mean,section_compiled(slopes_abs(:,col) == flats(col),col));
        end
        
        steadinessTimescale{currentClass} = time2steady;
        steadinessValue{currentClass} = steadiedGR;
        
    elseif currentClass == 2 % 60 min timescale
        
        clear steadiedGR time2steady section_compiled section_compiled_dataOnly slopes_dataOnly
        
        % isolate 60min experiment data
        % entire signal
        currentData = cell2mat( binned_mean(expClass==currentClass) );
        finalBin = length(currentData);
        
        
        % order signal such that 1st half is downshift (start to end of low phase),
        %                        2nd half is upshift (start to end of high phase)
        binOrder = [13:finalBin, 1:12];         % order for sorting
        orderedData = currentData(binOrder,:);
        
        % generate generic timesignal
        genericTime = ((1:finalBin)*75)';
        orderedTime = [genericTime, genericTime, genericTime];
        
        % from tpt zero to end, find fit a linear function and save slope
        clear slopes section_compiled
        for r = 1:24-5
            
            % select section to compute slope over
            section = orderedData(r:24,:);
            section_time = orderedTime(r:24,:);
            
            for col = 1:3
                section_slope = polyfit(section_time(:,col),section(:,col),1);
                slopes(r,col) = section_slope(1);
                section_compiled{r,col} = section(:,col);
            end
        end
        clear r section_slope col
        
        % find which slope is closest to zero
        slopes_abs = abs(slopes);
        flats = min(slopes_abs,[],1);
        
        % find times associated with flattest slope
        orderedTime_trim = orderedTime(1:23,:);
        for col = 1:3
            time2steady(col) = orderedTime_trim(slopes_abs(:,col) == flats(col));
            steadiedGR(col) = cellfun(@mean,section_compiled(slopes_abs(:,col) == flats(col),col));
        end
        
        steadinessTimescale{currentClass} = time2steady;
        steadinessValue{currentClass} = steadiedGR;
    
        
    elseif currentClass == 3
        
        %%
        % for downshifts, only one replicate reaches G_low
        % thus, only calculate "time to stable" as "time to G_low" 
        % confirm that between "time to G_low" and end that mean gr is G_low
        %
        
        clear steadiedGR section_compiled orderedData time2steady slopes slopes_abs orderedTime currentData section section_time orderedTime_trim
        
        
        % isolate single shift experiment data
        currentData_temp = binned_mean(expClass==currentClass); % entire signal
        

        % trim to start at upshift
        shiftTime_bin = ceil(shiftTime./timePerBin); % bin at t=zero. 75 seconds is first true bin after shift.
        for rep = 1:2
            
            % data for current replicate
            d = currentData_temp{rep};
            
            % bin of upshift in current replicate
            zeroBin = shiftTime_bin(rep);
            
            % trim to start at upshift
            d_zeroed = d(zeroBin:end);
            currentData_zeroed{rep} = d_zeroed;
            
        end
        clear d d_zeroed zeroBin rep
        
        
        % trim so datasets have equal lengths (though for single downshift they already do)
        numBins = cellfun(@length,currentData_zeroed);
        shorty = min(numBins);
        for rep = 1:2
            
            dz = currentData_zeroed{rep};
            dz(dz==0) = NaN;
            orderedData(:,rep) = dz(1:shorty);
            
        end
        clear currentData_temp currentData_zeroed rep dz
        
        
        
        % generate generic timesignal
        genericTime = ((1:shorty)*75)';
        orderedTime = [genericTime, genericTime];
        
        
        
        % from tpt zero to end, find growth rates within std of G_low
        orderedData(1,:) = NaN; % ignore data from 1st timepoint of each, as they are still coming down from high N
        near_G = find(orderedData > (G_low-std_low) );  
        
        
        % as expected near_G only finds frames within range of G_low in first single downshift, 2018-08-08
        % because NaNs do not appear consecutively, consider "steady-state low" 
        % as reached when indeces between near_G are consistently less than 3
        % 1st index allowed for NaN, 2nd for noise in growth signal
        
        idx_diffs = [NaN; diff(near_G)];
        allowables = idx_diffs < 3;
        last_jump = max(find(allowables == 0)); % last large gap in growth rates that hit G_low range
        stable_start = last_jump + 1;
        
        
        % from first "stable" data index in near_G to end, find mean slope
        clear slopes
        
        finalIndex = length(near_G);
        col = 1;

        % select section to compute slope over
        tpt_first = near_G(stable_start);
        gr = orderedData(tpt_first:end,col);
        t = orderedTime(tpt_first:end,col);
        
        % compute slope 
        idx = isnan(gr);
        section_slope = polyfit(t(~idx),gr(~idx),1);
        
        % store slope, mean growth rate, first timepoint of steaded range, and growth rate data within section
        slopes = section_slope(1);
        steadiedGR = nanmean(gr);
        time2steady = t(1);
        section_compiled{1,1} = gr;
        
        clear section_slope col idx gr t idx_diffs
        
        % store mean growth rate and first timepoint of steaded range across classes
        steadinessTimescale{currentClass} = time2steady;
        steadinessValue{currentClass} = steadiedGR;
        
    end
    
end

clear ans binOrder binOrder_sorted col counter slopes slopes_abs fullPeriod
clear flats orderedTime orderedTime_trim orderedData cl currentClass idx
clear shorty t section section_time numBins gr finalBin initialBin
clear genericTime steadiedGR


%% plot bar graphs of time to saturation and final saturation value

% time to saturation
t_sat_mean = cellfun(@mean,steadinessTimescale)./60;
t_sat_std = cellfun(@std,steadinessTimescale)./60;

sat_gr_mean = cellfun(@mean,steadinessValue);
sat_gr_std = cellfun(@std,steadinessValue);

figure(2)
hold on
bar(t_sat_mean)
errorbar(1:3,t_sat_mean,t_sat_std,'.')
ylabel('time (min)')

t_sat_mean % display values
t_sat_std


figure(3)
hold on
bar(sat_gr_mean)
errorbar(1:3,sat_gr_mean,sat_gr_std,'.')
ylabel('growth rate (1/hr)')

sat_gr_mean % display values
sat_gr_std







