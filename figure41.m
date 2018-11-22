%% figure 40

%  Goal: plot pdfs of birth size, raw and normalized.

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
%                       6. isolate isDrop, timestamp, volume and timeSinceBirth data
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
%  commit: raw pdf of birth volume for stable and 60 min replicates


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
binSize = 0.5;          % cubic microns


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
        
        
        
        % 6. isolate isDrop, timestamp, and volume data
        timestamps_hr = conditionData_fullOnly(:,2)./3600;    % col 2  = timestamp in seconds converted to hours
        isDrop = conditionData_fullOnly(:,4);                 % col 4  = isDrop, 1 marks a birth event
        volumes = conditionData_fullOnly(:,11);               % col 11 = Va volumes (cubic microns)
        clear curveIDs
        
        
        % 7. extract only final times from each growth curve
        birthSizes = volumes(isDrop==1);
        birthTimestamps = timestamps_hr(isDrop==1);
        

        
        % 8. truncate data to non-erroneous (e.g. bubbles) timestamps
        maxTime = bubbletime(condition);
        
        if maxTime > 0
            birthSize_bubbleTrimmed = birthSizes(birthTimestamps <= maxTime,:);
            birthTimestamps_bubbleTrimmed = birthTimestamps(birthTimestamps <= maxTime,:);
        else
            birthSize_bubbleTrimmed = divTimes;
            birthTimestamps_bubbleTrimmed = divTimestamps;
        end
        clear timestamps_hr maxTime isDrop
        
        
        
        % 10. truncate data to stabilized regions
        minTime = 3;
        birthSize_fullyTrimmed = birthSize_bubbleTrimmed(birthTimestamps_bubbleTrimmed >= minTime,:);
        clear volumes 
        
        
        % if no div data in steady-state, skip condition
        if isempty(birthSize_fullyTrimmed) == 1
            continue
        else
            
            
            % 11. trim outliers (those 3 std dev away from median) from final dataset
            birthS_median = median(birthSize_fullyTrimmed);
            birthS_std = std(birthSize_fullyTrimmed);
            birthS_temp = birthSize_fullyTrimmed(birthSize_fullyTrimmed <= (birthS_median+birthS_std*3)); % cut smallest vals, over 3 std out
            birthS_final = birthS_temp(birthS_temp >= (birthS_median-birthS_std*3));          % cut largest vals, over 3 std out
            

            
            % 12. calculate final populuation mean and count
            birthT_mean(e) = mean(birthS_final);
            birthT_count(e) = length(birthS_final);
            
            
            
            
            % 14. bin data and normalize bin counts by total counts
            assignedBins_raw = ceil(birthS_final/binSize);
            
            
            
            % 15. calculate pdf
            binned_birthS_raw = accumarray(assignedBins_raw, birthS_final, [], @(x) {x});
            binCounts_raw = cellfun(@length,binned_birthS_raw);
            pdf_birthS_raw = binCounts_raw/birthT_count(e);
            
            
            
            % 16. plot PDF
            palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
            color = rgb(palette(condition));
            
            % raw growth rates
            raw_vals = (1:max(assignedBins_raw)).*binSize;
            
            figure(1)
            plot(raw_vals,pdf_birthS_raw,'Color',color,'LineWidth',1)
            hold on
            title('raw pdf')
            legend('fluc','low','ave','high')
            xlabel('birth volume (cubic um)')
            ylabel('pdf')
            %axis([extremes(1) extremes(2) 0 .1])
            
        end
        
        % 17. repeat for each condition
    end
    
    
end

% 18. save plots in active folder
cd('/Users/jen/Documents/StockerLab/Data_analysis/currentPlots/')

figure(1)
plotName = strcat('figure41-',num2str(timescale),'-birthSize-rawPDF');
saveas(gcf,plotName,'epsc')

birthT_mean
birthT_count    

