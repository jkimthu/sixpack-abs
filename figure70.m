%% figure 70


%  Goals: plot growth rate over time, for each xy position separately

%  Strategy:
%
%       0. initialize directory and meta data
%       0. define time binning and growth rates of interest, see comments below for details 
%       0. define condition of interest
%       
%       1. determine experiments of interest, by index in storedMetaData
%       2. initialize experiment meta data
%       3. load measured experiment data
%       4. determine number of positions in specified condition
%       5. loop through xys and plot growth rate vs time
%          note: this procedure is from figure7.m



%  last updated: jen, 2019 Aug 6

%  commit: plot data from promixal, middle and distal xy individually

% OK let's go!

%%

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));


% 0. define growth rates of interest
specificGrowthRate = 'log2';
specificColumn = 3;
xmin = -1.5;
xmax = 4;


% 0. define binning parameters
specificBinning = 5;
binsPerHour = 60/specificBinning;


% 0. define condition of interest
palette = {'MediumBlue','DodgerBlue','LightSkyBlue','Indigo','Purple','Plum','DarkGoldenRod','GoldenRod','PaleGoldenRod','DarkRed','FireBrick','Crimson'};

        
counter = 0;
for condition = 1:4
 

    % 1. determine experiments of interest
    index = 15; % 2018-02-01 data
    
    
    
    % 2. initialize experiment meta data
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
    else
        filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
        %filename = strcat('lb-fluc-',date,'-width1p7-jiggle-0p5.mat');
        % single upshift and downshift data only uses larger width thresh
    end
    load(filename,'D5','T');
    
    
    
    % 4. determine number of positions in specified condition
    xy_start = storedMetaData{index}.xys(condition,1);
    xy_end = storedMetaData{index}.xys(condition,end);
    xy_middle = storedMetaData{index}.xys(condition,5);
    
    xy_list = [xy_start,xy_middle,xy_end];
    %xy_legend = cell(3,1);
    
    
    
    % 5. loop through xys and plot growth rate vs time
    for ii = 1:length(xy_list)
        
        counter = counter + 1;
        xy = xy_list(ii);
        %xy_legend{xy} = num2str(xy);
        
        % build data matrix from current xy position
        xyData = buildDM(D5, T, xy, xy,index,expType);
        
        
        
        % isolate volume (Va), timestamp, drop, curve, and trackNum data
        volumes = getGrowthParameter(xyData,'volume');             % volume = calculated va_vals (cubic um)
        timestamps_sec = getGrowthParameter(xyData,'timestamp');   % ND2 file timestamp in seconds
        isDrop = getGrowthParameter(xyData,'isDrop');              % isDrop == 1 marks a birth event
        curveFinder = getGrowthParameter(xyData,'curveFinder');    % col 5  = curve finder (ID of curve in condition)
        trackNum = getGrowthParameter(xyData,'trackNum');          % track number, not ID from particle tracking
        
        
        
        
        
        % calculate growth rate
        growthRates = calculateGrowthRate(volumes,timestamps_sec,isDrop,curveFinder,trackNum);
        
        
        
        
        % truncate data to non-erroneous (e.g. bubbles) timestamps
        maxTime = bubbletime(condition);
        timestamps_hr = timestamps_sec/3600; % time in seconds converted to hours
        
        if maxTime > 0
            xyData_bubbleTrimmed = xyData(timestamps_hr <= maxTime,:);
            growthRates_bubbleTrimmed = growthRates(timestamps_hr <= maxTime,:);
        else
            xyData_bubbleTrimmed = xyData;
            growthRates_bubbleTrimmed = growthRates;
        end
        clear timestamps_hr timestamps_sec maxTime
        
        
        
        
        % bin growth rate into time bins based on timestamp
        timeInHours = xyData_bubbleTrimmed(:,2)/3600;
        bins = ceil(timeInHours*binsPerHour);
        
        
        
        
        % isolate selected specific growth rate and remove nans from data analysis
        growthRt = growthRates_bubbleTrimmed(:,specificColumn);
        growthRt_noNaNs = growthRt(~isnan(growthRt),:);
        bins_noNaNs = bins(~isnan(growthRt),:);
        
        
        
        % calculate mean, standard dev, counts, and standard error
        binned_growthRt = accumarray(bins_noNaNs,growthRt_noNaNs,[],@(x) {x});
        bin_means = cellfun(@mean,binned_growthRt);
        bin_stds = cellfun(@std,binned_growthRt);
        bin_counts = cellfun(@length,binned_growthRt);
        bin_sems = bin_stds./sqrt(bin_counts);
        
%         if condition == 1
%             bin_means_compiled{xy} = bin_means; % for determine where noise(?) peaks are in time
%             bin_times{xy} = (1:length(bin_means))/binsPerHour;
%         end
%         
        
        % 13. plot growth rate over time
        color = rgb(palette(counter));
        xmark = '.';
        
        figure(1)
        plot((1:length(bin_means))/binsPerHour,bin_means,'Color',color)
        %plot((1:length(bin_means))/binsPerHour,bin_means,'Marker',xmark)
        hold on
        grid on
        axis([0,10.1,xmin,xmax])
        xlabel('Time (h)')
        ylabel('Growth rate')
        title(strcat(date,': (',specificGrowthRate,')'))
        %legend(xy_legend{1:xy})
        
    end

end
% 14. save plots in active folder
%cd('/Users/jen/Documents/StockerLab/Data_analysis/currentPlots/')
%plotName = strcat('figure7-',specificGrowthRate,'-',date,'-',num2str(specificBinning),'minbins-variedWidth-constjiggle-allNOTonlyFull');
%saveas(gcf,plotName,'epsc')

%close(gcf)

clc
% 15. repeat for all experiments




