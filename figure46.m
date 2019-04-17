%% figure 46 - mean length over time


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
%       7. isolate length, volume and timestamp
%       8. if appropriate, assign NaN to all growth rates associated with frames to ignore
%          else simply remove existing nans from analysis
%       9. bin growth rate into time bins based on timestamp
%      10. calculate mean, standard dev, counts, and standard error
%      11. plot growth rate over time
%      12. save plots in active folder
% 13. repeat for all experiments


%  last updated: jen, 2019 April 17

%  commit: first commit, length and volume vs time with 2 min bins


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
        
        
        
        
        % 7. isolate selected length, volume and timestamp
        majorAx = getGrowthParameter(conditionData_bubbleTrimmed,'length');                % length (um)
        volume = getGrowthParameter(conditionData_bubbleTrimmed,'volume');                 % volume (va) cubic microns
        timestamps_sec = getGrowthParameter(conditionData_bubbleTrimmed,'timestamp');      % ND2 file timestamp in seconds
        timeInHours = timestamps_sec./3600;
        clear conditionData timestamps_sec
        
        
        
        
        % 8. if appropriate, assign NaN to all growth rates associated with frames to ignore
        %     else simply remove existing nans from analysis
        if condition == 1 && isempty(ignoredFrames) == 0
            
            
            % assign NaN to frames to ignore
            frameNum = conditionData_bubbleTrimmed(:,16); % col 16 = original frame number
            length_ignorant = majorAx;
            volume_ignorant = volume;
            for fr = 1:length(ignoredFrames)
                length_ignorant(frameNum == ignoredFrames(fr),1) = NaN;
                volume_ignorant(frameNum == ignoredFrames(fr),1) = NaN;
            end
            
            
            % remove nans from data analysis
            length_noNaNs = length_ignorant(~isnan(length_ignorant),:);
            volume_noNaNs = volume_ignorant(~isnan(volume_ignorant),:);
            timeInHours_noNans = timeInHours(~isnan(length_ignorant),:);
            
            
            
        else
            
            % remove nans from data analysis
            length_noNaNs = majorAx(~isnan(majorAx),:);
            volume_noNaNs = volume(~isnan(volume),:);
            timeInHours_noNans = timeInHours(~isnan(majorAx),:);
            
        end
        clear frameNum conditionData_trim2 majorAx volume fr
        
        
        
        
        
        % 9. bin size into time bins based on timestamp
        bins = ceil(timeInHours_noNans*binsPerHour);
        binned_length = accumarray(bins,length_noNaNs,[],@(x) {x});
        binned_volume = accumarray(bins,volume_noNaNs,[],@(x) {x});

        
        
        
        % 10. calculate mean, standard dev, counts, and standard error
        L_means = cellfun(@mean,binned_length);
        L_stds = cellfun(@std,binned_length);
        L_counts = cellfun(@length,binned_length);
        
        V_means = cellfun(@mean,binned_volume);
        V_stds = cellfun(@std,binned_volume);
        V_counts = cellfun(@length,binned_volume);

        
        
        
        % 11. plot growth rate over time
        palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick','LimeGreen','MediumPurple'};
        

        color = rgb(palette(condition));
        xmark = '.';
        

        
        figure(1)
        plot((1:length(L_means))/binsPerHour,L_means,'Color',color,'Marker',xmark)
        hold on
        legend('fluc','low','ave','high')
        xlabel('Time (hr)')
        ylabel('Length (um)')
        title(strcat(date))
        
        figure(2)
        plot((1:length(V_means))/binsPerHour,V_means,'Color',color,'Marker',xmark)
        hold on
        legend('fluc','low','ave','high')
        xlabel('Time (hr)')
        ylabel('Volume (cubic um)')
        title(strcat(date))
        

    end
    
    
    % 12. save plots in active folder
    cd('/Users/jen/Documents/StockerLab/Data_analysis/currentPlots/')
    
    figure(1)
    plotName = strcat('figure46-lengthvtime-',date,'-',num2str(specificBinning),'minbins');
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
    figure(2)
    plotName = strcat('figure46-volumevtime-',date,'-',num2str(specificBinning),'minbins');
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
    
    
    % 13. repeat for all experiments
end


%%

