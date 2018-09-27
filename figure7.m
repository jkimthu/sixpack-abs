%% figure 7


%  Goals: plot growth rate over time

%  Strategy:
%
%       0. initialize directory and meta data
%       0. define time binning and growth rates of interest, see comments below for details 
%       1. create array of experiments of interest, for each:
%               2. initialize experiment meta data
%               3. 
%               4. load measured data
%               5. for each condition in current experiment, build data matrix from specified condition
%                       6. isolate condition data to those with full cell cycles
%                       7. truncate data to non-erroneous (e.g. bubbles) timestamps
%                       8. isolate volume (Va), timestamp, mu, drop and curveID data
%                       9. calculate growth rate
%                      10. bin growth rate into time bins based on timestamp
%                      11. isolate selected specific growth rate and remove nans from data analysis
%                      12. calculate mean, standard dev, counts, and standard error
%              13. plot growth rate over time
%              14. save plot with experiment #, specific growth rate definition, and binning          
%      13. repeat for all experiments 


%  last updated: jen, 2018 September 26

%  commit: base 2 on fresh re-do of 2018-01-29 analysis, finish 60 min reps 

         

% OK let's go!

%%

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));


% 0. define growth rates of interest, see comments below for details
prompt = 'Enter specific growth rate definition as string (raw / norm / log2 / lognorm): ';
specificGrowthRate = input(prompt);

prompt = 'Enter specific binning in minutes as double: ';
specificBinning = input(prompt);
binsPerHour = 60/specificBinning;


%%
% 1. create array of experiments of interest, then loop through each:
exptArray = 13; % use corresponding dataIndex values

for e = 1:length(exptArray)
    
    
    % 2. initialize experiment meta data
    index = exptArray(e); %dataIndex(e);
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
    
    
    
    % 5. build data matrix from specified condition
    for condition = 1:length(bubbletime)
    
        xy_start = storedMetaData{index}.xys(condition,1);
        xy_end = storedMetaData{index}.xys(condition,end);
        conditionData = buildDM(D5, T, xy_start, xy_end,index,expType);
        
        
        
        % 6. isolate condition data to those with full cell cycles
        curveIDs = conditionData(:,5);           % col 5 = curve ID
        conditionData_fullOnly = conditionData(curveIDs > 0,:);
        clear curveFinder
        
        
  
        
        % 7. isolate volume (Va), timestamp, mu, drop and curveID data
        volumes = conditionData_fullOnly(:,11);        % col 11 = calculated va_vals (cubic um)
        timestamps_sec = conditionData_fullOnly(:,2);  % col 2  = timestamp in seconds
        isDrop = conditionData_fullOnly(:,4);          % col 4  = isDrop, 1 marks a birth event
        curveFinder = conditionData_fullOnly(:,5);     % col 5  = curve finder (ID of curve in condition)
        trackNum = conditionData_fullOnly(:,20);       % col 20 = track number (not ID from particle tracking)
        
        
        
        % 8. calculate growth rate
        growthRates = calculateGrowthRate(volumes,timestamps_sec,isDrop,curveFinder,trackNum);


        
        
        % 9. truncate data to non-erroneous (e.g. bubbles) timestamps
        maxTime = bubbletime(condition);
        timestamps_hr = conditionData_fullOnly(:,2)/3600; % time in seconds converted to hours
        
        if maxTime > 0
            conditionData_bubbleTrimmed = conditionData_fullOnly(timestamps_hr <= maxTime,:);
            growthRates_bubbleTrimmed = growthRates(timestamps_hr <= maxTime,:);
        else
            conditionData_bubbleTrimmed = conditionData_fullOnly;
            growthRates_bubbleTrimmed = growthRates;
        end
        clear timestamps_hr timestamps_sec maxTime
        
        
        
        
        % 10. bin growth rate into time bins based on timestamp
        timeInHours = conditionData_bubbleTrimmed(:,2)/3600;
        bins = ceil(timeInHours*binsPerHour);
        %binVector = linspace(1,binsPerHour*10,binsPerHour*10);
        
        
        
        % 11. isolate selected specific growth rate and remove nans from data analysis
        if strcmp(specificGrowthRate,'raw') == 1
            specificColumn = 1;         % for selecting appropriate column in growthRates
            xmin = -5;                  % lower limit for plotting x axis
            xmax = 25;                  % upper limit for plotting x axis
        elseif strcmp(specificGrowthRate,'norm') == 1
            specificColumn = 2;
            xmin = -1;
            xmax = 5;
        elseif strcmp(specificGrowthRate,'log2') == 1
            specificColumn = 3;
            xmin = -1;
            xmax = 5;
        elseif strcmp(specificGrowthRate,'lognorm') == 1
            specificColumn = 4;
            xmin = -0.5;
            xmax = 1;
        end
        
        %growthRt = growthRates(:,specificColumn);
        growthRt = growthRates_bubbleTrimmed(:,specificColumn);
        
        growthRt_noNaNs = growthRt(~isnan(growthRt),:);
        bins_noNaNs = bins(~isnan(growthRt),:);
        
        
        
        % 12. calculate mean, standard dev, counts, and standard error
        binned_growthRt = accumarray(bins_noNaNs,growthRt_noNaNs,[],@(x) {x});
        bin_means = cellfun(@mean,binned_growthRt);
        bin_stds = cellfun(@std,binned_growthRt);
        bin_counts = cellfun(@length,binned_growthRt);
        bin_sems = bin_stds./sqrt(bin_counts);
       
       
        
        % 13. plot growth rate over time
        palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
        
        color = rgb(palette(condition));
        xmark = 'o';
        
        figure(e)
        errorbar((1:length(bin_means))/binsPerHour,bin_means,bin_sems,'Color',color)
        hold on
        plot((1:length(bin_means))/binsPerHour,bin_means,'Color',color,'Marker',xmark)
        hold on
        grid on
        axis([0,10.1,xmin,xmax])
        xlabel('Time (hr)')
        ylabel('Growth rate')
        title(strcat(date,': (',specificGrowthRate,')'))
        
        
    end
    
    % 14. save plots in active folder
    cd('/Users/jen/Documents/StockerLab/Data_analysis/currentPlots/')
    plotName = strcat('figure7-',specificGrowthRate,'-',date,'-',num2str(specificBinning),'minbins');
    saveas(gcf,plotName,'epsc')
    
    close(gcf)
    
% 15. repeat for all experiments 
end


%%

