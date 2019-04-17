%% figure 47 - mean width over time


%  Goals: plot width over time, like growth rate

%  Strategy:
%
% 0. initialize complete meta data
% 0. define growth rates and bin size (time) of interest
% 1. create array of experiments of interest, then loop through each:
%       2. initialize experiment meta data
%       3. load measured experiment data
%       4. for single shift experiments, define which frames to ignore (noisy tracking)
%          these frames to be ignored only apply to condition 1.
%       5. build data matrix from specified condition
%       6. truncate data to non-erroneous (e.g. bubbles) timestamps
%       7. isolate width and timestamp
%       8. if appropriate, assign NaN to all growth rates associated with frames to ignore
%          else simply remove existing nans from analysis
%       9. bin growth rate into time bins based on timestamp
%      10. calculate mean, standard dev, counts, and standard error
%      11. plot growth rate over time
%      12. save plots in active folder
% 13. repeat for all experiments



%  last updated: jen, 2019 April 17

%  commit: first commit, population width vs time with 2 min bins


% OK let's go!

%%

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));


% 0. define growth rates and bin size (time) of interest
specificBinning = 2;
binsPerHour = 60/specificBinning;


%%
% 1. create array of experiments of interest, then loop through each:
exptArray = [2:4,5:7,9:12,13:15]; % use corresponding dataIndex values

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
    if strcmp(date,'2017-11-09') == 1
        filename = strcat('lb-control-',date,'-width1p4-jiggle-0p5.mat');
    elseif strcmp(date,'2017-09-26') == 1
        filename = 'lb-monod-2017-09-26-jiggle-c12-0p1-c3456-0p5-bigger1p8.mat';
    elseif strcmp(date,'2018-12-04') == 1
        filename = 'lb-monod-2018-12-04-c12-width1p7-c34-width1p4-jiggle-0p5.mat';
    elseif strcmp(expType,'origFluc') == 1
        filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
    elseif strcmp(expType,'fluc2stable') == 1
        filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
    else
        filename = strcat('lb-fluc-',date,'-width1p7-jiggle-0p5.mat');
        % single upshift and downshift data only uses larger width thresh
    end
    load(filename,'D5','T');
    
    
    
    % 4. for single shift experiments, define which frames to ignore (noisy tracking)
    %    these frames to be ignored only apply to condition 1.
    if strcmp(date,'2018-06-15') == 1
        ignoredFrames = [112,113,114];
    elseif strcmp(date,'2018-08-01') == 1
        ignoredFrames = [94,95];
    elseif strcmp(date,'2018-08-09') == 1
        ignoredFrames = [115,116,117];
    elseif strcmp(date,'2018-08-08') == 1
        ignoredFrames = [112,113,114];
    else
        ignoredFrames = [];
    end
    
    
    
    % 5. build data matrix from specified condition
    for condition = 1:length(bubbletime)
        
        
        xy_start = storedMetaData{index}.xys(condition,1);
        xy_end = storedMetaData{index}.xys(condition,end);
        conditionData = buildDM(D5, T, xy_start, xy_end,index,expType);
        
        
        % 6. truncate data to non-erroneous (e.g. bubbles) timestamps
        maxTime = bubbletime(condition);
        timestamps_hr = conditionData(:,2)/3600; % time in seconds converted to hours
        
        if maxTime > 0
            conditionData_bubbleTrimmed = conditionData(timestamps_hr <= maxTime,:);
        else
            conditionData_bubbleTrimmed = conditionData;
        end
        clear timestamps_hr
        
        
        
        
        % 7. isolate selected width and timestamp
        width = getGrowthParameter(conditionData_bubbleTrimmed,'width');                % width (um)
        timestamps_sec = getGrowthParameter(conditionData_bubbleTrimmed,'timestamp');   % ND2 file timestamp in seconds
        timeInHours = timestamps_sec./3600;
        clear conditionData timestamps_sec
        
        
        
        
        % 8. if appropriate, assign NaN to all growth rates associated with frames to ignore
        %     else simply remove existing nans from analysis
        if condition == 1 && isempty(ignoredFrames) == 0
            
            
            % assign NaN to frames to ignore
            frameNum = conditionData_bubbleTrimmed(:,16); % col 16 = original frame number
            width_ignorant = width;
            for fr = 1:length(ignoredFrames)
                width_ignorant(frameNum == ignoredFrames(fr),1) = NaN;
            end
            
            
            % remove nans from data analysis
            width_noNaNs = width_ignorant(~isnan(width_ignorant),:);
            timeInHours_noNans = timeInHours(~isnan(width_ignorant),:);
            
            
            
        else
            
            % remove nans from data analysis
            width_noNaNs = width(~isnan(width),:);
            timeInHours_noNans = timeInHours(~isnan(width),:);
            
        end
        clear frameNum conditionData_trim2 width fr
        
        
        
        
        
        % 9. bin growth rate into time bins based on timestamp
        bins = ceil(timeInHours_noNans*binsPerHour);

        
        
        
        % 10. calculate mean, standard dev, counts, and standard error
        binned_width = accumarray(bins,width_noNaNs,[],@(x) {x});
        bin_means = cellfun(@mean,binned_width);
        bin_stds = cellfun(@std,binned_width);
        bin_counts = cellfun(@length,binned_width);
        bin_sems = bin_stds./sqrt(bin_counts);
        
        
        
        % 11. plot growth rate over time
        palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick','LimeGreen','MediumPurple'};
        

        color = rgb(palette(condition));
        xmark = '.';
        

        
        figure(e)
        plot((1:length(bin_means))/binsPerHour,bin_means,'Color',color,'Marker',xmark)
        hold on
        legend('fluc','low','ave','high')
        axis([0,10,1.1,1.5])
        xlabel('Time (hr)')
        ylabel('Width (um)')
        title(strcat(date))
        
        bin_counts(bin_counts == 0) = NaN;
        min(bin_counts)
        max(bin_counts)

    end
    
    
    % 12. save plots in active folder
    cd('/Users/jen/Documents/StockerLab/Data_analysis/currentPlots/')
    plotName = strcat('figure47-widthvtime-',date,'-',num2str(specificBinning),'minbins');
    saveas(gcf,plotName,'epsc')
    
    close(gcf)
    
    
    % 13. repeat for all experiments
end


%%

