%% figure 7


%  Goals: plot growth rate over time

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


%  last updated: jen, 2019 July 1

%  commit: first growth rate over time for 2019-06-26 experiment


% OK let's go!

%%

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));


% 0. define growth rates and bin size (time) of interest
specificGrowthRate = 'log2';
specificColumn = 3;

specificBinning = 2;
binsPerHour = 60/specificBinning;


%%
% 1. create array of experiments of interest, then loop through each:
exptArray = 38; % use corresponding dataIndex values

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
    elseif strcmp(date,'2017-11-09') == 1
        filename = strcat('lb-control-',date,'-width1p4-jiggle-0p5.mat');
    elseif strcmp(date,'2017-09-26') == 1
        filename = 'lb-monod-2017-09-26-jiggle-c12-0p1-c3456-0p5-bigger1p8.mat';
    elseif strcmp(date,'2018-12-04') == 1
        filename = 'lb-monod-2018-12-04-c12-width1p7-c34-width1p4-jiggle-0p5.mat';
    elseif strcmp(expType,'origFluc') == 1
        filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
    elseif strcmp(expType,'steady2fluc') == 1
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
        
        if strcmp(expType,'steady2fluc') == 1
            xy_start = storedMetaData{index}.xys{condition}(1);
            xy_end = storedMetaData{index}.xys{condition}(end);
        else
            xy_start = storedMetaData{index}.xys(condition,1);
            xy_end = storedMetaData{index}.xys(condition,end);
        end
        conditionData = buildDM(D5, T, xy_start, xy_end,index,expType);
        
        
        
        % 6. isolate volume (Va), timestamp, drop, curve, and trackNum data
        volumes = getGrowthParameter(conditionData,'volume');             % volume = calculated va_vals (cubic um)
        timestamps_sec = getGrowthParameter(conditionData,'timestamp');   % ND2 file timestamp in seconds
        isDrop = getGrowthParameter(conditionData,'isDrop');              % isDrop == 1 marks a birth event
        curveFinder = getGrowthParameter(conditionData,'curveFinder');    % col 5  = curve finder (ID of curve in condition)
        trackNum = getGrowthParameter(conditionData,'trackNum');          % track number, not ID from particle tracking
        
        
        
        % 7. calculate growth rate
        growthRates = calculateGrowthRate(volumes,timestamps_sec,isDrop,curveFinder,trackNum);
        clear volumes isDrop curveFinder trackNum
        
        

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
        timestamps_sec = getGrowthParameter(conditionData_bubbleTrimmed,'timestamp'); % ND2 file timestamp in seconds
        timeInHours = timestamps_sec./3600;
        clear conditionData timestamps_sec
        
        
        
        % 10. if appropriate, assign NaN to all growth rates associated with frames to ignore
        %     else simply remove existing nans from analysis
        if condition == 1 && isempty(ignoredFrames) == 0
            
            
            % assign NaN to frames to ignore
            frameNum = conditionData_bubbleTrimmed(:,16); % col 16 = original frame number
            growthRt_ignorant = growthRt;
            for fr = 1:length(ignoredFrames)
                growthRt_ignorant(frameNum == ignoredFrames(fr),1) = NaN;
            end
            
            
            % remove nans from data analysis
            growthRt_noNaNs = growthRt_ignorant(~isnan(growthRt_ignorant),:);
            timeInHours_noNans = timeInHours(~isnan(growthRt_ignorant),:);
            
            
            
        else
            
            % remove nans from data analysis
            growthRt_noNaNs = growthRt(~isnan(growthRt),:);
            timeInHours_noNans = timeInHours(~isnan(growthRt),:);
            
        end
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
        palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick','LimeGreen','MediumPurple'};
        

        color = rgb(palette(condition));
        xmark = '.';
        
%         % adjust ymax in plot to longest bubble time
%         if ymax < maxTime
%             ymax = maxTime;
%         end
        
        figure(e)
        plot((1:length(bin_means))/binsPerHour,bin_means,'Color',color,'Marker',xmark)
        hold on
        grid on
        legend('high,untreated','high,treated','low,untreated','low,treated')
        %axis([2,ymax+0.1,xmin,xmax])
        axis([0,10,-0.5,3.5])
        xlabel('Time (hr)')
        ylabel('Growth rate')
        title(strcat(date,': (',specificGrowthRate,')'))
        
        bin_counts(bin_counts == 0) = NaN;
        min(bin_counts)
        max(bin_counts)
        %end
    end
    
    
    % 14. save plots in active folder
    %cd('/Users/jen/Documents/StockerLab/Data_analysis/currentPlots/')
    %plotName = strcat('figure7-',specificGrowthRate,'-',date,'-',num2str(specificBinning),'minbins');
    %saveas(gcf,plotName,'epsc')
    
    %close(gcf)
    %clc
    
    % 15. repeat for all experiments
end


%%

