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
%                       6. isolate isDrop, timestamp, and timeSinceBirth data
%                       7. extract only final times from each growth curve
%                       8. remove zeros, which occur if no full track data exists at a drop        
%                       9. truncate data to non-erroneous (e.g. bubbles) timestamps
%                      10. truncate data to stabilized regions
%                      11. trim outliers (those 3 std dev away from median) from final dataset
%                      12. trim interdivision times < 5 min
%                      13. calculate final populuation mean and count
%                      14. bin data and normalize bin counts by total counts
%                      15. calculate pdf
%                      16. plot PDF
%             17. repeat for each condition
%      18. save plot in active folder




%  Last edit: jen, 2018 November 22
%  commit: raw pdf of interdivision time for stable and 30 sec replicates,
%          edit comments


% OK! Lez go!


%%
% 0. initialize data & specify target concentration

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));


% 0. define bin size
binSize = 5;          % min


%%
% 1. create array of experiments of interest, then loop through each:
exptArray = [2,3,4]; % use corresponding dataIndex values


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
    if strcmp(date,'2017-11-12') == 1
        filename = strcat('lb-fluc-',date,'-width1p4-jiggle-0p5.mat');
    else
        filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
    end
    load(filename,'D5','T');
    
    
    
    % 4. build data matrix from specified condition
    for condition = 1:length(bubbletime)
        
        xy_start = storedMetaData{index}.xys(condition,1);
        xy_end = storedMetaData{index}.xys(condition,end);
        conditionData = buildDM(D5, T, xy_start, xy_end,index,expType);
        
        
        
        % 5. isolate condition data to those with full cell cycles
        curveIDs = conditionData(:,5);           % col 5 = curve ID
        conditionData_fullOnly = conditionData(curveIDs > 0,:);
        clear curveFinder conditionData
        
        
        
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
        clear timestamps_hr maxTime isDrop birthIndeces finalDivTimes
        
        
        
        % 10. truncate data to stabilized regions
        minTime = 3;
        divTimes_fullyTrimmed = divTimes_bubbleTrimmed(divTimestamps_bubbleTrimmed >= minTime,:);
        divTimestamps_fullyTrimmed = divTimestamps_bubbleTrimmed(divTimestamps_bubbleTrimmed >= minTime,:);
        
        
        % if no div data in steady-state, skip condition
        if isempty(divTimes_fullyTrimmed) == 1
            continue
        else
            
            
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
            
            
            
            % 15. calculate pdf
            binned_divT_raw = accumarray(assignedBins_raw, divT_final, [], @(x) {x});
            binCounts_raw = cellfun(@length,binned_divT_raw);
            pdf_divT_raw = binCounts_raw/divT_count;
            
            
            
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
            
        end
        
        % 17. repeat for each condition
    end
    
    
end

% 18. save plots in active folder
cd('/Users/jen/Documents/StockerLab/Data_analysis/currentPlots/')

figure(1)
plotName = strcat('figure40-',num2str(timescale),'interdivT-rawPDF');
saveas(gcf,plotName,'epsc')


    

