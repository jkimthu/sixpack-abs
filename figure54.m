%% figure 54
%
%  Goals: from upshift data, calculate
%                   (1) time to stabilization and
%                   (2) mean growth rate during stable region


%  Strategy:
%
%         1. find window of stabilized region that gives fit slope closest to zero
%         2. time until window of zero-most slope = time to stabilization 
%         3. mean growth rate within window = mean growth rate of stabilized signal



%  Summary of code sections:
%
%         Part A. initialize folder with stored meta data
%         Part B. curves for fluctuating data
%         Part C. curve for single shift data
%         Part D. fitting for quantifications
%         Part E. plot quantifications


%  last updated: jen, 2019 Jan 22

%  commit: edit to streamline for upshifts
%          

% OK let's go!


%% Part A. initialize

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

timePerBin = 75; % matches binning in shift response plots

clear prompt


%% Part B. accumulate upshift data


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
        filename = 'lb-fluc-2017-11-12-width1p4-jiggle-0p5.mat'; % says only 1.4 but it is both
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
    clear timeInPeriods timeInPeriodFraction
    
    
    % 12. find which bins are boundaries signal phases
    lastBin_Q1 = (timescale/timePerBin)/4;                      % last bin before downshift
    lastBin_downshift = (timescale*3/4)/timePerBin;             % last bin before upshift
    
    firstBin_upshift = (timescale*3/4)/timePerBin + 1;          % first bin of upshift
    lastBin_ofPeriod = timescale/timePerBin;                    % total bins in signal period
    
    
    
    
    % 13. list bins chronologically to combine broken up high nutrient phase
    %       i.e. start of upshift is Q4, concatenate Q1
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
        color = rgb('Aquamarine');
    elseif timescale == 3600
        color = rgb('Indigo');
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
    
    
    figure(1)   % binned upshift
    plot(upshift_times,upshift_growth,'Color',color,'LineWidth',1,'Marker','.')
    hold on
    title(strcat('response to upshift, binned every (',num2str(timePerBin),') sec'))
    
    clear x y ycenter ybegin yend color
    clear conditionData_trim1 conditionData_trim2 growthRates_trim1 growthRates_trim2
    clear growthRt_noNaNs index lastBin_Q1 lastBin_ofPeriod lastBin_downshift firstBin_upshift firstBin_downshift
    clear first_pre_downshiftBin date condition e
    
end

xlabel('time (min)')
ylabel(strcat('growth rate: (', specificGrowthRate,')'))
axis([numPreshiftBins*-1*timePerBin_min,80,xmin,xmax])


clear bins_unique bubbletime correctedTime_noNans filename numPreshiftBins
clear timestamps_hr postUpshift_growth postUpshift_times
clear preUpshift_growth preUpshift_times preShift_bins upshift_growth upshift_times
clear timePerBin_min


%% Part C. accumulate single shift data


% 1. create array of experiments of interest, then loop through each
exptArray = [21,22]; % dataIndex values of single upshift experiments


%counter = 0;  % keep counter value from part B and continue
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
    
    
    % 12. remove nans from data analysis
    growthRt_noNaNs = growthRt_ignorant(~isnan(growthRt_ignorant),:);
    correctedTime_noNans = correctedTime(~isnan(growthRt_ignorant),:);
    clear growthRt growthRt_ignorant correctedTime frameNum
    
    
    
    
    % 13. assign corrected timestamps to bins, by which to accumulate growth data
    bins = ceil(correctedTime_noNans/timePerBin);      % bin 1 = first timePerBin sec of experiment
    bins_unique = (min(bins):max(bins))';              % avoids missing bins due to lack of data
    
    
    % find bins after shift
    % generalized for single shift experiments
    first_postshiftBin_single = ceil(shiftTime(e_shift)/timePerBin) + 1; % shift occurs at shiftTime/timePerBin, so first full bin with shifted data is the next one
    postshiftBins_single{counter} = (first_postshiftBin_single:max(bins))';
    
    
    
    
    % 14. choose which pre-shift data bins to plot
    
    % single shift experiments don't have high/low phase interruptions
    % however, they do have bins missing data!
    numPreshiftBins = 10;


    % determine pre-shift bins
    index_single_shift = find(bins_unique == first_postshiftBin_single);
    pre_upshiftBins{counter} = bins_unique(index_single_shift-numPreshiftBins : index_single_shift-1); % same bins in both single down and upshifts
    %pre_downshiftBins{counter} = bins_unique(index_single_shift-numPreshiftBins : index_single_shift-1);
    

    
    
    % 15. collect growth rate data into bins and calculate stats
    %     WARNING: bins variable does not i
    binned_growthRate{counter} = accumarray(bins,growthRt_noNaNs,[],@(x) {x});
    binned_mean{counter} = accumarray(bins,growthRt_noNaNs,[],@mean);
    clear bins
    
   
    
    
    % 17. plot response in growth rate for all timescales over time
    timePerBin_min = timePerBin/60; % time per bin in minutes
    color = rgb('DarkOliveGreen');

    
%     if strcmp(shiftType,'upshift') == 1
        
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
        
        
        grid on
        hold on
        title(strcat('response to upshift, binned every (',num2str(timePerBin),') sec'))

        
        
%     else
%         
%         figure(3)    % downshift
%         
%         % plot!
%         preDownshift_times = ((numPreshiftBins-1)*-1:0)*timePerBin_min;
%         postDownshift_times = (1:length(binned_mean{counter}(postshiftBins_single{counter})))*timePerBin_min;
%         downshift_times_gapped = [preDownshift_times, postDownshift_times];
%         
%         preDownshift_growth_gapped = binned_mean{counter}(pre_downshiftBins{counter});
%         postDownshift_growth_single = binned_mean{counter}(postshiftBins_single{counter});
%         downshift_growth_gapped = [preDownshift_growth_gapped; postDownshift_growth_single];
%         
%         % don't plot zeros that are place holders for gaps in data
%         downshift_growth = downshift_growth_gapped(downshift_growth_gapped > 0);
%         downshift_times = downshift_times_gapped(downshift_growth_gapped > 0);
%         
%         plot(downshift_times,downshift_growth,'Color',color,'LineWidth',1)
%         grid on
%         hold on
%         title(strcat('response to downshift, binned every (',num2str(timePerBin),') sec'))
%              
        
%     end

     
end


xlabel('time (min)')
ylabel(strcat('growth rate: (', specificGrowthRate,')'))
axis([numPreshiftBins*-1*timePerBin_min,160,xmin,xmax])


clear bubbletime color conditionData_trim2 date correctedTime_noNans
clear bins_unique e_shift experimentFolder exptArray expType filename
clear fr growthRt_noNaNs ignoredFrames index index_single_shift
clear numPreshiftBins postshiftBins_single preUpshift_times preUpshift_growth
clear shiftType specificColumn specificGrowthRate timePerBin_min
clear timescale times_trim1 upshift_growth upshift_growth_gapped
clear xmax xmin first_postshiftBin_single postUpshift_growth_single
clear postUpshift_times upshift_times_gapped upshift_times


%% Part D. fitting slopes to find time until stable
%
% .. and once stable region is determined, calculate mean steady growth rate

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
        
        % order signal such that 1st half is upshift (start to end of high phase),
        %                        2nd half is downshift (start to end of low phase)
        binOrder = [4:finalBin, 1:3];
        [~,binOrder_sorted] = sort(binOrder);           % get the order for sorting
        orderedData = currentData(binOrder_sorted,:);   % sort data so that first high after low starts the array
        
        % generate generic timesignal
        genericTime = ((1:finalBin)*75)';
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
        
        % order signal such that 1st half is upshift, 2nd half is downshift
        binOrder = [13:finalBin, 1:12];
        [~,binOrder_sorted] = sort(binOrder); % get the order of binOrder
        orderedData = currentData(binOrder_sorted,:);
        
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
        
        clear steadiedGR section_compiled orderedData time2steady slopes slopes_abs orderedTime currentData section section_time orderedTime_trim
        
        shiftTime_bin = ceil(shiftTime./75); % bin at t=zero. 75 seconds is first true bin after shift.
        
        
        % isolate single shift experiment data
        % entire signal
        currentData_temp = binned_mean(expClass==currentClass);
        finalBin = 128; % 75 sec/bin * 128 bins =  9600 sec  = 160 min
                        % reflects final of analysis which starts at shift
        initialBin = 64; % 80 min * 60 sec/min / 75 sec/bin = 64 bin
                        
                        
        % trim to start at upshift
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
        
        
        % trim so datasets have equal lengths
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
        
        
        % from tpt zero to end, find fit a linear function and save slope
        clear slopes
        for r = initialBin:finalBin-5
            
            % select section to compute slope over
            section = orderedData(r:finalBin,:);
            section_time = orderedTime(r:finalBin,:);
            
            for col = 1:2
                
                gr = section(:,col);
                t = section_time(:,col);
                idx = isnan(gr);
                section_slope = polyfit(t(~idx),gr(~idx),1);
                
                %section_slope = polyfit(section_time(:,col),section(:,col),1);
                slopes(r,col) = section_slope(1);
                section_compiled{r,col} = gr;
            end
            
        end
        clear r section_slope col
        
        % find which slope is closest to zero
        slopes_abs = abs(slopes(initialBin:end,:));
        section_compiled_abs = section_compiled(initialBin:end,:);
        flats = min(slopes_abs,[],1);
        
        % find times associated with flattest slope
        orderedTime_trim = orderedTime(initialBin:end,:);
        for col = 1:2
            time2steady(col) = orderedTime_trim(slopes_abs(:,col) == flats(col));
            steadiedGR(col) = cellfun(@nanmean,section_compiled_abs(slopes_abs(:,col) == flats(col),col));
        end
        
        steadinessTimescale{currentClass} = time2steady;
        steadinessValue{currentClass} = steadiedGR;
        
    end
    
end

clear ans binOrder binOrder_sorted col counter slopes slopes_abs fullPeriod
clear flats orderedTime orderedTime_trim orderedData cl currentClass idx
clear shorty t section section_time numBins gr finalBin initialBin
clear genericTime steadiedGR


%% Part E. plot bar graphs of time to saturation and final saturation value

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



