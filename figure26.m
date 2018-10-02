%% figure 26

%  Goal: plot pdfs of growth rates, raw and normalized.

%        (A) trim data prior to three hours, raw
%        (B) trim data prior to three hours, divide by population mean
%        (c) trim data prior to three hours, subtract population mean


%  Strategy:
%
%       0. initialize data & specify target concentration
%       1. create a directory of experiments with target concentration
%       2. for all experiments in target directory... accumulate cell size and curve duration data
%               3. move to experiment folder and build data matrix
%               4. for each condition with target concentration...
%                       5. build data maxtrix from data for current condition
%                       6. isolate volume(Va), curve duration, added volume per cell cycle, drop and time data
%                       7. isolate only data during which drop == 1 (birth event)
%                               - birth size = volume at birth event
%                               - one value for duration and added volume is gathered per cell cycle
%                                 (durations and added volume are assembled
%                                 in data matrix as final values, repeated
%                                 for all timepoints in a curve)
%                       8. remove data from cell cycles with added volume > 0
%                       9. remove data from cell cycles with durations shorter than 10 min (zero values are incomplete)
%                      10. trim data to stabilized / non-bubble timestamps
%                      11. trim outliers (those 3 std dev away from media) from final dataset
%                      12. calculate count number of data points per bin
%                      13. bin data and normalize bin counts by total counts
%                      14. plot pdf and histograms per experiment
%              15. repeat for each condition, storing fluc and stable data in their own variables
%      16.  after assembling data from all experiments, plot violin plots
%           comparing distrubtions between fluc and stable for each experiment
%


%  Last edit: jen, 2018 Oct 1
%  commit: pdfs of growth rate for stable and 60 min replicates, include
%          subtracted mean



% OK! Lez go!


%%
% 0. initialize data & specify target concentration

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));


% 0. define growth rates of interest, see comments below for details
prompt = 'Enter specific growth rate definition as string (raw / norm / log2 / lognorm): ';
specificGrowthRate = input(prompt);


% 0. define bin size
binSize = 0.1;          % 1/hr
extremes = [-10,12];    % 1/hr

%%
% 1. create array of experiments of interest, then loop through each:
exptArray = [13,14,15]; % use corresponding dataIndex values


for e = 1:length(exptArray)
    
    
    % 2. initialize experiment meta data
    index = exptArray(e);
    date = storedMetaData{index}.date;
    expType = storedMetaData{index}.experimentType;
    bubbletime = storedMetaData{index}.bubbletime;
    
    timescale = storedMetaData{index}.timescale;
    disp(strcat(date, ': analyze!'))
    
    
    
    % 3. load measured experiment data
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-width1p7-jiggle-0p5.mat');
    load(filename,'D5','T');
    
    
    
    % 4. build data matrix from specified condition
    for condition = 1:length(bubbletime)
        
        xy_start = storedMetaData{index}.xys(condition,1);
        xy_end = storedMetaData{index}.xys(condition,end);
        conditionData = buildDM(D5, T, xy_start, xy_end,index,expType);
        
        
        
        % 5. isolate condition data to those with full cell cycles
        curveIDs = conditionData(:,5);           % col 5 = curve ID
        conditionData_fullOnly = conditionData(curveIDs > 0,:);
        clear curveFinder
        
        
        
        % 6. isolate volume (Va), timestamp, mu, drop and curveID data
        volumes = conditionData_fullOnly(:,11);        % col 11 = calculated va_vals (cubic um)
        timestamps_sec = conditionData_fullOnly(:,2);  % col 2  = timestamp in seconds
        isDrop = conditionData_fullOnly(:,4);          % col 4  = isDrop, 1 marks a birth event
        curveFinder = conditionData_fullOnly(:,5);     % col 5  = curve finder (ID of curve in condition)
        trackNum = conditionData_fullOnly(:,20);       % col 20 = track number (not ID from particle tracking)
        
        
        
        % 7. calculate growth rate
        growthRates = calculateGrowthRate(volumes,timestamps_sec,isDrop,curveFinder,trackNum);
        
        
        
        % 8. truncate data to non-erroneous (e.g. bubbles) timestamps
        maxTime = bubbletime(condition);
        timestamps_hr = conditionData_fullOnly(:,2)/3600; % time in seconds converted to hours
        
        if maxTime > 0
            conditionData_bubbleTrimmed = conditionData_fullOnly(timestamps_hr <= maxTime,:);
            growthRates_bubbleTrimmed = growthRates(timestamps_hr <= maxTime,:);
        else
            conditionData_bubbleTrimmed = conditionData_fullOnly;
            growthRates_bubbleTrimmed = growthRates;
        end
        clear timestamps_sec timestamps_hr maxTime
        
        
        
        % 9. truncate data to stabilized regions
        minTime = 3;
        timestamps_hr = conditionData_bubbleTrimmed(:,2)/3600; % time in seconds converted to hours
        
        conditionData_fullyTrimmed = conditionData_bubbleTrimmed(timestamps_hr >= minTime,:);
        growthRates_fullyTrimmed = growthRates_bubbleTrimmed(timestamps_hr >= minTime,:);
        clear timestamps_hr
        
        
        % 10. isolate selected specific growth rate and remove nans from data analysis
        if strcmp(specificGrowthRate,'raw') == 1
            specificColumn = 1;         % for selecting appropriate column in growthRates
        elseif strcmp(specificGrowthRate,'norm') == 1
            specificColumn = 2;
        elseif strcmp(specificGrowthRate,'log2') == 1
            specificColumn = 3;
        elseif strcmp(specificGrowthRate,'lognorm') == 1
            specificColumn = 4;
        end
        
        growthRt = growthRates_fullyTrimmed(:,specificColumn);
        growthRt_noNaNs = growthRt(~isnan(growthRt),:);
        clear gr_median gr_std gr_temp

        
        
        
        % 11. trim outliers (those 3 std dev away from median) from final dataset
        gr_median = median(growthRt_noNaNs);
        gr_std = std(growthRt_noNaNs);
        gr_temp = growthRt_noNaNs(growthRt_noNaNs <= (gr_median+gr_std*3)); % cut smallest vals, over 3 std out
        growthRt_final = gr_temp(gr_temp >= (gr_median-gr_std*3));          % cut largest vals, over 3 std out
        
        
        
        
        % 12. remove growth rate data that occurs during a switch
        %
%         if condition == 1
%             
%             % extract and trim corrected timestamp in sec
%             correctedTime = conditionData_fullyTrimmed(:,22);
%             correctedTime_noNans = correctedTime(~isnan(growthRt),:);
%             ct_temp = correctedTime_noNans(growthRt_noNaNs <= (gr_median+gr_std*3));
%             correctedTime_final = ct_temp(gr_temp >= (gr_median-gr_std*3)); 
%             
%             
%             % compute nutrient signal, where 1 = high and 0 = low (2 steps)
%             
%             % step 1 of 1: translate timestamps into quarters of nutrient signal
%             timeInPeriods = correctedTime_final/timescale;        % unit = sec/sec
%             timeInPeriodFraction = timeInPeriods - floor(timeInPeriods);
%             timeInPeriodQuarters = ceil(timeInPeriodFraction * 4);
%             
%             % step 2 of 2: from nutrient signal quarters, generate a binary nutrient signal
%             binaryNutrientSignal = zeros(length(timeInPeriodQuarters),1);
%             binaryNutrientSignal(timeInPeriodQuarters == 1) = 1;
%             binaryNutrientSignal(timeInPeriodQuarters == 4) = 1;
%             
%             
%             % locate indeces with a nutrient switch (up or down)
%             switches = diff(binaryNutrientSignal);
%             
%             %  assign corrected timestamps to bins, by which to accumulate growth rate data
%             timeInPeriodFraction_inSeconds = timeInPeriodFraction * timescale;
%             timeInPeriodFraction_inBins = ceil(timeInPeriodFraction_inSeconds/timePerBin);
%             
%             
%         end
        

        
        % 13. calculate final populuation mean and count
        gr_mean = mean(growthRt_final);
        gr_count = length(growthRt_final);
        
        
        
        % 14. bin data and normalize bin counts by total counts
        assignedBins_raw = ( ceil(growthRt_final * (1/binSize)) );
        assignedBins_raw_shiftup = assignedBins_raw + 1000;    % shift right by 1000,
        % because accumarray cannot deal with negative valued bins
        % correct for this shift later, in step _____
        
        assignedBins_norm = ( ceil( (growthRt_final./gr_mean) * (1/binSize)) );
        assignedBins_norm_shiftup = assignedBins_norm + 1000;
        
        assignedBins_sub = ( ceil( (growthRt_final-gr_mean) * (1/binSize)) );
        assignedBins_sub_shiftup = assignedBins_sub + 1000;
        
        
        
        % 15. calculate pdf, still shifted up 1000
        binned_gr_raw = accumarray(assignedBins_raw_shiftup, growthRt_final, [], @(x) {x});
        binCounts_raw = cellfun(@length,binned_gr_raw);
        pdf_gr_raw = binCounts_raw/gr_count;
        
        binned_gr_norm = accumarray(assignedBins_norm_shiftup, assignedBins_norm, [], @(x) {x});
        binCounts_norm = cellfun(@length,binned_gr_norm);
        pdf_gr_norm = binCounts_norm/gr_count;
        
        binned_gr_sub = accumarray(assignedBins_sub_shiftup, assignedBins_sub, [], @(x) {x});
        binCounts_sub = cellfun(@length,binned_gr_sub);
        pdf_gr_sub = binCounts_sub/gr_count;
        
        
        
        % 16. plot PDF
        palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
        color = rgb(palette(condition));
        
        % raw growth rates
        gr_raw = ((1:length(binned_gr_raw))-1000).*binSize;
        
        figure(1)
        plot(gr_raw,pdf_gr_raw,'Color',color,'LineWidth',1)
        hold on
        title('raw pdf')
        legend('fluc','low','ave','high')
        xlabel('growth rate (1/hr)')
        ylabel('pdf')
        axis([extremes(1) extremes(2) 0 .05])
        
        
        % divided by population mean
        gr_norm = ((1:length(binned_gr_norm))-1000).*binSize;
        
        figure(2)
        plot(gr_norm,pdf_gr_norm,'Color',color,'LineWidth',1)
        hold on
        title('norm pdf, div')
        legend('fluc','low','ave','high')
        xlabel('gr/<gr>')
        ylabel('pdf')
        axis([extremes(1) extremes(2) 0 .15])
        
        
        
        % subtracted by population mean
        gr_sub = ((1:length(binned_gr_sub))-1000).*binSize;
        
        figure(3)
        plot(gr_sub,pdf_gr_sub,'Color',color,'LineWidth',1)
        hold on
        title('norm pdf, sub')
        legend('fluc','low','ave','high')
        xlabel('gr-<gr>')
        ylabel('pdf')
        axis([extremes(1) extremes(2) 0 .05])
        

    end
    
    
end

% 17. save plots in active folder
cd('/Users/jen/Documents/StockerLab/Data_analysis/currentPlots/')

figure(1)
plotName = strcat('figure26-',num2str(timescale),'-',specificGrowthRate,'-rawPDF');
saveas(gcf,plotName,'epsc')
%close(gcf)

figure(2)
plotName = strcat('figure26-',num2str(timescale),'-',specificGrowthRate,'-divPDF');
saveas(gcf,plotName,'epsc')
%close(gcf)

figure(3)
plotName = strcat('figure26-',num2str(timescale),'-',specificGrowthRate,'-subPDF');
saveas(gcf,plotName,'epsc')
%close(gcf)


        

    

