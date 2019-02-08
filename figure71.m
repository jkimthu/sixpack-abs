%% figure 71


%  Goals: plot growth rate over time for sensitivity analysis
%         
%         1) one plot of all four conditions per varied width threshold
%         2) one four-panel subplot: one panel per condition with all thresholds overlaid


%  Strategy:
%
%       0. initialize directory and meta data
%       0. define time binning and growth rates of interest, see comments below for details 
%       1. create array of experiments of interest, for each:
%               2. initialize experiment meta data
%               3. load measured experiment data    
%               4. for single shift experiments, define which frames to ignore (noisy tracking)
%               5. for each condition in current experiment, build data matrix from specified condition
%                       6. isolate volume (Va), timestamp, drop, curve, and trackNum data
%                       7. calculate growth rate
%                       8. truncate data to non-erroneous (e.g. bubbles) timestamps
%                       9. isolate selected specific growth rate and timestamp
%                      10. if appropriate, assign NaN to all growth rates associated with frames to ignore
%                          else simply remove existing nans from analysis
%                      11. bin growth rate into time bins based on timestamp
%                      12. calculate mean, standard dev, counts, and standard error
%                      13. plot growth rate over time
%              14. save plot with experiment #, specific growth rate definition, and binning          
%      15. repeat for all experiments 


%  last updated: jen, 2018 December 10

%  commit: first commit, sensitivity analysis for width threshold with 2018-02-01 data


% OK let's go!

%% A1. initialize

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
cd('/Users/jen/Documents/StockerLab/Data_analysis/sensitivity analysis')

% 0. define growth rate metric and binning
specificGrowthRate = 'log2';
specificColumn = 3;
xmin = -1;
xmax = 5;
ymax = 0;

specificBinning = 30; % mins per bin
binsPerHour = 60/specificBinning;




%% B1. loop through width thresholds and plot all four conditions together

% 1. create array of experiments of interest, then loop through each:
widthArray = [4,5,6,7]; % use corresponding dataIndex values

for w = 1:length(widthArray)
    
    
    % 2. initialize experiment meta data
    index = 15; % stored meta data index for experiment 2018-02-01
    date = storedMetaData{index}.date;
    expType = storedMetaData{index}.experimentType;
    bubbletime = storedMetaData{index}.bubbletime;
    
    timescale = storedMetaData{index}.timescale;
    disp(strcat(date, ': analyze!'))
    
    
    
    % 3. load measured experiment data
    filename = strcat('lb-fluc-',date,'-width1p',num2str(widthArray(w)),'-jiggle-0p5.mat');
    load(filename,'D5','T');
    
    
    
    % 4. build data matrix from specified condition
    for condition = 1:length(bubbletime)
        
        xy_start = storedMetaData{index}.xys(condition,1);
        xy_end = storedMetaData{index}.xys(condition,end);
        conditionData = buildDM(D5, T, xy_start, xy_end,index,expType);
        
        
        
        % 5. isolate volume (Va), timestamp, drop, curve, and trackNum data
        volumes = conditionData(:,11);        % col 11 = calculated va_vals (cubic um)
        timestamps_sec = conditionData(:,2);  % col 2  = timestamp in seconds
        isDrop = conditionData(:,4);          % col 4  = isDrop, 1 marks a birth event
        curveFinder = conditionData(:,5);     % col 5  = curve finder (ID of curve in condition)
        trackNum = conditionData(:,20);       % col 20 = track number (not ID from particle tracking)
        
        
        
        % 6. calculate growth rate
        growthRates = calculateGrowthRate(volumes,timestamps_sec,isDrop,curveFinder,trackNum);
        
        

        % 7. truncate data to non-erroneous (e.g. bubbles) timestamps
        maxTime = bubbletime(condition);
        timestamps_hr = conditionData(:,2)/3600; % time in seconds converted to hours
        
        if maxTime > 0
            conditionData_bubbleTrimmed = conditionData(timestamps_hr <= maxTime,:);
            growthRates_bubbleTrimmed = growthRates(timestamps_hr <= maxTime,:);
        else
            conditionData_bubbleTrimmed = conditionData;
            growthRates_bubbleTrimmed = growthRates;
        end
        clear timestamps_hr timestamps_sec
        
        
        
        % 8. isolate selected specific growth rate and timestamp
        growthRt = growthRates_bubbleTrimmed(:,specificColumn);
        timeInHours = conditionData_bubbleTrimmed(:,2)/3600;   % col 2 = raw timestamps
        clear isDrop trackNum volumes curveFinder conditionData
        
        
        
        % 9. remove existing nans from analysis
        growthRt_noNaNs = growthRt(~isnan(growthRt),:);
        timeInHours_noNans = timeInHours(~isnan(growthRt),:);
        
        clear frameNum conditionData_trim2 growthRates fr
        clear growthRt growthRates_bubbleTrimmed
        
        
        
        
        
        % 10. bin growth rate into time bins based on timestamp
        bins = ceil(timeInHours_noNans*binsPerHour);

        
        
        
        % 11. calculate mean, standard dev, counts, and standard error
        binned_growthRt = accumarray(bins,growthRt_noNaNs,[],@(x) {x});
        bin_means = cellfun(@mean,binned_growthRt);
        bin_stds = cellfun(@std,binned_growthRt);
        bin_counts = cellfun(@length,binned_growthRt);
        bin_sems = bin_stds./sqrt(bin_counts);
        
        
        
        % 12. plot growth rate over time
        palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick','LimeGreen','MediumPurple'};
        palette_width = {'SeaGreen','SlateBlue','Indigo','BlueViolet'};
        condition_legend = {'fluc','low','ave','high'};

        color = rgb(palette(condition));
        width_color = rgb(palette_width(w));
        xmark = '.';
        
        % adjust ymax in plot to longest bubble time
        if ymax < maxTime
            ymax = maxTime;
        end
        
        
        % four conditions with single width threshold 
        figure(1)
        plot((1:length(bin_means))/binsPerHour,bin_means,'Color',color,'Marker',xmark)
        hold on
        grid on
        legend('high,untreated','high,treated','low,untreated','low,treated')
        axis([0,10,xmin,xmax])
        xlabel('Time (hr)')
        ylabel('Growth rate')
        title(strcat(date,': (',specificGrowthRate,')'))
        
        

    end

    % 14. save plots in active folder
    cd('/Users/jen/Documents/StockerLab/Data_analysis/sensitivity analysis')
    figure(1)
    plotName = strcat('figure71-',specificGrowthRate,'-',date,'-width1p',num2str(widthArray(w)),'-2minbins');
    saveas(gcf,plotName,'epsc')
    close(gcf)
    clc
    
    % 15. repeat for all experiments
end


%% C1. loop through conditions and plot all four width thresholds together

% 1. create array of experiments of interest, then loop through each:
widthArray = [4,5,6,7]; % use corresponding dataIndex values


% 2. initialize experiment meta data
index = 15; % stored meta data index for experiment 2018-02-01
date = storedMetaData{index}.date;
expType = storedMetaData{index}.experimentType;
bubbletime = storedMetaData{index}.bubbletime;

timescale = storedMetaData{index}.timescale;
disp(strcat(date, ': analyze!'))

% 3. loop through conditions and plot growth rate data
for condition = 1:4
    
    for w = 1:length(widthArray)
    
        
    % 4. load measured experiment data
    filename = strcat('lb-fluc-',date,'-width1p',num2str(widthArray(w)),'-jiggle-0p5.mat');
    cd('/Users/jen/Documents/StockerLab/Data_analysis/sensitivity analysis')
    load(filename,'D5','T');
    
    
    % 5. build data matrix from specified condition
    xy_start = storedMetaData{index}.xys(condition,1);
    xy_end = storedMetaData{index}.xys(condition,end);
    conditionData = buildDM(D5, T, xy_start, xy_end,index,expType);
    
    
    % 6. isolate volume (Va), timestamp, drop, curve, and trackNum data
    volumes = conditionData(:,11);        % col 11 = calculated va_vals (cubic um)
    timestamps_sec = conditionData(:,2);  % col 2  = timestamp in seconds
    isDrop = conditionData(:,4);          % col 4  = isDrop, 1 marks a birth event
    curveFinder = conditionData(:,5);     % col 5  = curve finder (ID of curve in condition)
    trackNum = conditionData(:,20);       % col 20 = track number (not ID from particle tracking)
    
    
    % 7. calculate growth rate
    growthRates = calculateGrowthRate(volumes,timestamps_sec,isDrop,curveFinder,trackNum);
        
    
    % 8. truncate data to non-erroneous (e.g. bubbles) timestamps
    maxTime = bubbletime(condition);
    timestamps_hr = conditionData(:,2)/3600; % time in seconds converted to hours
    
    if maxTime > 0
        conditionData_bubbleTrimmed = conditionData(timestamps_hr <= maxTime,:);
        growthRates_bubbleTrimmed = growthRates(timestamps_hr <= maxTime,:);
    else
        conditionData_bubbleTrimmed = conditionData;
        growthRates_bubbleTrimmed = growthRates;
    end
    clear timestamps_hr timestamps_sec
    
    
    
    % 9. isolate selected specific growth rate and timestamp
    growthRt = growthRates_bubbleTrimmed(:,specificColumn);
    timeInHours = conditionData_bubbleTrimmed(:,2)/3600;   % col 2 = raw timestamps
    clear isDrop trackNum volumes curveFinder conditionData
    
    
    % 10. remove existing nans from analysis
    growthRt_noNaNs = growthRt(~isnan(growthRt),:);
    timeInHours_noNans = timeInHours(~isnan(growthRt),:);
    
    clear frameNum conditionData_trim2 growthRates fr
    clear growthRt growthRates_bubbleTrimmed
    
    
    % 11. bin growth rate into time bins based on timestamp
    bins = ceil(timeInHours_noNans*binsPerHour);
    
    
    
    % 12. calculate mean, standard dev, counts, and standard error
    binned_growthRt = accumarray(bins,growthRt_noNaNs,[],@(x) {x});
    bin_means = cellfun(@mean,binned_growthRt);
    bin_stds = cellfun(@std,binned_growthRt);
    bin_counts = cellfun(@length,binned_growthRt);
    bin_sems = bin_stds./sqrt(bin_counts);
    
    
    
    % 13. plot growth rate over time
    palette_width = {'GoldenRod','SeaGreen','Indigo','Orchid'};
    condition_legend = {'fluc','low','ave','high'};
    
    width_color = rgb(palette_width(w));
    xmark = '.';
    
    % four conditions with single width threshold
    figure(condition)
    plot((1:length(bin_means))/binsPerHour,bin_means,'Color',width_color,'Marker',xmark)
    hold on
    grid on
    legend('1.4','1.5','1.6','1.7')
    axis([0,10,xmin,xmax])
    xlabel('Time (hr)')
    ylabel('Growth rate')
    title(condition_legend(condition))
    
    
    end
    
end
    
% 14. save plots in active folder
cd('/Users/jen/Documents/StockerLab/Data_analysis/sensitivity analysis')
figure(1)
axis([0,9.1,-0.8,3])
plotName = strcat('figure71-compiled-fluc-',num2str(specificBinning),'minbins');
saveas(gcf,plotName,'epsc')
close(gcf)

figure(2)
axis([0,9.1,0,1.4])
plotName = strcat('figure71-compiled-low-',num2str(specificBinning),'minbins');
saveas(gcf,plotName,'epsc')
close(gcf)

figure(3)
axis([0,9.1,1,4.5])
plotName = strcat('figure71-compiled-ave-',num2str(specificBinning),'minbins');
saveas(gcf,plotName,'epsc')
close(gcf)

figure(4)
axis([0,9.1,0.5,5.5])
plotName = strcat('figure71-compiled-high-',num2str(specificBinning),'minbins');
saveas(gcf,plotName,'epsc')
close(gcf)

clc

%%

    
    % 15. repeat for all experiments



%%
% four conditions with all four thresholds
figure(2)
hold on
subplot(2,2,condition)
plot((1:length(bin_means))/binsPerHour,bin_means,'Color',width_color,'Marker',xmark)
xlabel('Time (hr)')
ylabel('Growth rate')
title(condition_legend(condition))


        
cd('/Users/jen/Documents/StockerLab/Data_analysis/sensitivity analysis')
figure(2)
plotName = strcat('figure71-',specificGrowthRate,'-',date,'-rangeOfWidths');
saveas(gcf,plotName,'epsc')

