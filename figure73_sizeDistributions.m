%%  Figure 73. size distributions from a single timepoint

%  Goal: what is the distribution of cell size at any given time?


%  Strategy:
%
%  Part A:
%     0. initialize complete meta data

%  Part B:
%     1. for all experiments in dataset:
%           2. collect experiment data
%           3. load measured data
%           4. for each condition,
%                   5. gather specified condition data
%                   6. isolate timestamp and truncate data to non-erroneous (e.g. bubbles) timestamps
%                   7. isolate selected specific growth rate
%                   8. bin growth rate into time bins based on timestamp
%                   9. plot volume distribution at three timepoints (early, mid, late)
%          10. repeat for all conditions
%    11. repeat for all experiments


%  Last edit: jen, 2019 July 25

%  commit: first commit, volume distributions at various timepoints across
%          experiment



%% (A) initialize analysis
clc
clear

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')

dataIndex = find(~cellfun(@isempty,storedMetaData));
experimentCount = length(dataIndex);

%% (B) plot volume over time, overlaying nutrient signal

% 1. for all experiments in dataset
exptArray = 15; % 2019-02-01

for e = 1:length(exptArray)
    
    % 2. collect experiment data
    index = exptArray(e);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    expType = storedMetaData{index}.experimentType;
    xys = storedMetaData{index}.xys;
    bubbletime = storedMetaData{index}.bubbletime;
    disp(strcat(date, ': analyze!'))
    
    
    % 3. load measured data
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
    load(filename,'D5','T');
    xy_start = min(min(xys));
    xy_end = max(max(xys));
    exptData = buildDM(D5, T, xy_start, xy_end,index,expType);

    
    % 4. for each condition, specify:
    for condition = 1:4 % 1 = fluctuating; 3 = ave nutrient condition
        
       
        % 5. gather specified condition data
        conditionData = exptData(exptData(:,21) == condition,:);
        
        
        % 6. isolate timestamp and truncate data to non-erroneous (e.g. bubbles) timestamps
        maxTime = bubbletime(condition);
        timestamps_sec = getGrowthParameter(conditionData,'timestamp');   % ND2 file timestamp in seconds
        timestamps_hr = timestamps_sec/3600; % time in seconds converted to hours
        
        if maxTime > 0
            conditionData_bubbleTrimmed = conditionData(timestamps_hr <= maxTime,:);
        else
            conditionData_bubbleTrimmed = conditionData;
        end
        clear timestamps_sec timestamps_hr
        
        
        % 7. isolate selected specific growth rate
        volumes = getGrowthParameter(conditionData_bubbleTrimmed,'volume');  % volume = calculated va_vals (cubic um)
        
        

        % 8. bin growth rate into time bins based on timestamp
        binsPerHour = 12;
        timestamps_sec = getGrowthParameter(conditionData_bubbleTrimmed,'timestamp');   % ND2 file timestamp in seconds
        timestamps_hr = timestamps_sec/3600; % time in seconds converted to hours
        bins = ceil(timestamps_hr*binsPerHour);
        binned_volume = accumarray(bins,volumes,[],@(x) {x});

        
        
        % 9. plot volume distribution at three timepoints (early, mid, late)
        tpt = [min(bins); floor(mean( [min(bins),max(bins)] ));max(bins)];
        palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick','LimeGreen','MediumPurple'};
        color = rgb(palette(condition));

        
        figure(condition)
        for tt = 1:length(tpt)
            subplot(1,3,tt)
            histogram(binned_volume{tt},'FaceColor',color)
            title(num2str(tpt(tt)))
            hold on
        end
        xlabel('Volume (cubic microns)')
        
        % 10. repeat for all conditions
    end
    % 11. repeat for all experiments
end



