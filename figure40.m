%% figure 40

%  Goal: plot pdfs of interdivision time, raw and normalized.

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


%  Last edit: jen, 2018 November 22
%  commit: raw pdf of interdivision time for stable and 60 min replicates


% OK! Lez go!


%%
% 0. initialize data & specify target concentration

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));

% 
% % 0. define growth rates of interest, see comments below for details
% prompt = 'Enter specific growth rate definition as string (raw / norm / log2 / lognorm): ';
% specificGrowthRate = input(prompt);


% 0. define bin size
binSize = 5;          % min
extremes = [0,120];  % 1/hr

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
    filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
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
        
        
        
        % 6. isolate isDrop, timestamp, and timeSinceBirth data
        timestamps_hr = conditionData_fullOnly(:,2)./3600;    % col 2  = timestamp in seconds converted to hours
        isDrop = conditionData_fullOnly(:,4);                 % col 4  = isDrop, 1 marks a birth event
        timeSinceBirth_min = conditionData_fullOnly(:,6)./60; % col 6  = timeSinceBirth in seconds converted to min
        
        
        % 7. extract only final times from each growth curve
        birthIndeces = find(isDrop==1); 
        finalDivTimes = timeSinceBirth_min(birthIndeces(2:end)-1); % final timeSince birth is indexed right before every drop event
        finalTimestamps = timestamps_hr(birthIndeces(2:end)-1);

        
        % 8. remove zeros, which occur if no full track data exists at a drop
        divTimes = finalDivTimes(finalDivTimes>0);
        divTimestamps = finalTimestamps(finalDivTimes>0);
        
        
        % 9. truncate data to non-erroneous (e.g. bubbles) timestamps
        maxTime = bubbletime(condition);
        
        if maxTime > 0
            divTimes_bubbleTrimmed = divTimes(divTimestamps <= maxTime,:);
            divTimestamps_bubbleTrimmed = divTimestamps(divTimestamps <= maxTime,:);
        else
            divTimes_bubbleTrimmed = divTimes;
            divTimestamps_bubbleTrimmed = divTimestamps;
        end
        clear timestamps_hr maxTime isDrop birthIndeces
        
        
        
        % 10. truncate data to stabilized regions
        minTime = 3;
        divTimes_fullyTrimmed = divTimes_bubbleTrimmed(divTimestamps_bubbleTrimmed >= minTime,:);
        divTimestamps_fullyTrimmed = divTimestamps_bubbleTrimmed(divTimestamps_bubbleTrimmed >= minTime,:);
        
        
        
        
        % 11. trim outliers (those 3 std dev away from median) from final dataset
        divT_median = median(divTimes_fullyTrimmed);
        divT_std = std(divTimes_fullyTrimmed);
        divT_temp = divTimes_fullyTrimmed(divTimes_fullyTrimmed <= (divT_median+divT_std*3)); % cut smallest vals, over 3 std out
        divT_final = divT_temp(divT_temp >= (divT_median-divT_std*3));          % cut largest vals, over 3 std out
        
        
        
        
        % 12. trim interdivision times < 5 min
        divT_final = divT_final(divT_final > 5);
        
        
        
        
        % 13. calculate final populuation mean and count
        divT_mean = mean(divT_final);
        divT_count = length(divT_final);
        
        
        
        
        % 14. bin data and normalize bin counts by total counts
        assignedBins_raw = ceil(divT_final/binSize);

        divT_norm = divT_final./divT_mean;
        assignedBins_norm = ceil(divT_norm/(1/binSize));
        
        divT_sub = divT_final-divT_mean;
        assignedBins_sub = ceil(divT_sub/(1/binSize));
        assignedBins_sub_shiftup = assignedBins_sub + 1000;
        
        
        
        % 14. calculate pdf
        binned_divT_raw = accumarray(assignedBins_raw, divT_final, [], @(x) {x});
        binned_divT_norm = accumarray(assignedBins_norm, divT_norm, [], @(x) {x});
        binned_divT_sub = accumarray(assignedBins_sub_shiftup, divT_sub, [], @(x) {x}); % shift bins up to avoid negative vals
        
        binCounts_raw = cellfun(@length,binned_divT_raw);
        binCounts_norm = cellfun(@length,binned_divT_norm);
        binCounts_sub = cellfun(@length,binned_divT_sub);
        
        pdf_divT_raw = binCounts_raw/divT_count;
        pdf_divT_norm = binCounts_norm/divT_count;
        pdf_divT_sub = binCounts_sub/divT_count;
        
        
        
        % 16. plot PDF
        palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
        color = rgb(palette(condition));
        
        % raw growth rates
        raw_vals = (1:max(assignedBins_raw)).*binSize;
        
        figure(1)
        plot(raw_vals,pdf_divT_raw,'Color',color,'LineWidth',1)
        hold on
        title('raw pdf')
        legend('fluc','low','ave','high')
        xlabel('interdivision time (min)')
        ylabel('pdf')
        %axis([extremes(1) extremes(2) 0 .1])
        
        
        % divided by population mean
%         norm_vals = (min(assignedBins_norm):max(assignedBins_norm)).*(1/binSize);
%         
%         figure(2)
%         plot(norm_vals,pdf_divT_norm,'Color',color,'LineWidth',1)
%         hold on
%         title('norm pdf, div')
%         legend('fluc','low','ave','high')
%         xlabel('divT/<divT>')
%         ylabel('pdf')
%         %axis([extremes(1) extremes(2) 0 .15])
%         
%         
%         
%         % subtracted by population mean
%         gr_sub = ((1:length(binned_gr_sub))-1000).*(1/binSize);
%         
%         figure(3)
%         plot(gr_sub,pdf_gr_sub,'Color',color,'LineWidth',1)
%         hold on
%         title('norm pdf, sub')
%         legend('fluc','low','ave','high')
%         xlabel('gr-<gr>')
%         ylabel('pdf')
%         axis([extremes(1) extremes(2) 0 .05])
        

    end
    
    
end

% 17. save plots in active folder
cd('/Users/jen/Documents/StockerLab/Data_analysis/currentPlots/')

figure(1)
plotName = strcat('figure40-',num2str(timescale),'interdivT-rawPDF');
saveas(gcf,plotName,'epsc')
%close(gcf)
% 
% figure(2)
% plotName = strcat('figure26-',num2str(timescale),'-',specificGrowthRate,'-divPDF');
% saveas(gcf,plotName,'epsc')
% %close(gcf)
% 
% figure(3)
% plotName = strcat('figure26-',num2str(timescale),'-',specificGrowthRate,'-subPDF');
% saveas(gcf,plotName,'epsc')
% %close(gcf)
% 
% 
        

    

