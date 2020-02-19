% relativeRates.m

% goal: testing the relative rates model (Harris & Theriot 2016) by
%       asking plotting growth rates of SA and V.

%       the model hypothesizes that biosynthesis rates of volume (alpha) and
%       surface components (beta) are proprotional to cell volume, V(t).
%       more formally:
%
%               beta/alpha = relatively constant
%
%       where,           
%               alpha = (dV/dt)/V(t)
%               beta  = (dA/dt)/V(t)
%
%
%       compellingly, this relationship brings about another hypothesis:
%       throughout the cell cycle, cells drop in SA/V as they grow UNTIL
%       septation begins to occur, which then brings SA/V back up. with
%       relative rates remaining equal, this suggests an accumulation of
%       cell surface materials until surface constriction is initiated,
%       plausibly due to some threshold accumulation.


% open questions:

%       1. is the relative rates model true in our system?

%       2. why/how might bacteria regulate their width (wider or thinner)
%          to acheive SA/V ratios consistent with relative rates model?

%       3. if division trigger is due to surface material accumulation,
%          what sets this threshold for each steady-state?


% strategy:
%
%   Part I:  SA/V across nutrient concentrations and timescales (population-level)
%                   (A)  collect and store data
%                   (B)  plot!
%
%   Part II: SA/V response to nutrient shifts (population-level)
%
%   Part III: comparison of volume biosynthesis rate (alpha) and SA biosynthesis rate (beta)
%             responses to nutrient shifts
%                   - is beta/alpha ratio is constant before and after shifts?
%                     if yes, this supports "relative rates" model
%                   - is beta/alpha ratio immediate the beta/alpha ratio of high and low? 
%                     what would this mean??


% last update: jen, 2018 March 26

% commit: plot SA/V across nutrient concentration, timescale and it's
% response to nutrient shifts. for comparison to Harris and Theriot, Cell
% 2016


% OK let's go!


%% Part I: SA/V across nutrient concentrations and timescales
%     (A): compiled and store data 

% strategy:

%       0. initialize experiment data
%       1. for each experiment...load data
%               -  identify date
%               -  exclude outlier: 2017-10-31
%               -  build data matrix of current experimental data
%               2. for each condition in experiment...
%                       3. isolate all data from current condition
%                       4. isolate condition specific raw timestamps for bubble trimming
%                       5. remove data with timestamps prior to and after stabilization
%                       6. isolate size at birth data
%                       7. isolate condition specific length, width, V and SA
%                       8. calculate average and s.e.m. of stablized datas
%                       9. accumulate data for storage / plotting
%              10. store data from all conditions into measured data structure
%      11.  save stored data into data structure
%      12.  plot average biovolume production rate over time


% 0. initialize experiment data
clear
clc

% initialize summary vectors for calculated data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
sizeData = cell(size(storedMetaData));
dataIndex = find(~cellfun(@isempty,storedMetaData));
experimentCount = length(dataIndex);


% 1. for each experiment, move to folder and load data
for e = 1:experimentCount
    
    % identify experiment by date
    index = dataIndex(e);
    date = storedMetaData{index}.date;
    
    % move directory to experiment data
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    
    % load data
    timescale = storedMetaData{index}.timescale;
    if ischar(timescale) == 0
        filename = strcat('lb-fluc-',date,'-window5-width1p4-1p7-jiggle-0p5.mat');
    elseif strcmp(date,'2017-09-26') == 1
        filename = 'lb-monod-2017-09-26-window5-va-jiggle-c12-0p1-c3456-0p5-bigger1p8.mat';
    elseif strcmp(date, '2017-11-09') == 1
        filename = 'lb-control-2017-11-09-window5-width1p4-jiggle-0p5.mat';
    end
    load(filename,'D5','M','M_va','T')
    
    % build experiment data matrix
    display(strcat('Experiment (', num2str(e),') of (', num2str(length(dataIndex)),')'))
    xy_start = 1;
    xy_end = length(D5);
    exptData = buildDM(D5,M,M_va,T,xy_start,xy_end,e);
    
    clear D5 M M_va T filename experimentFolder
   
    
    % 2. for each condition in experiment...
    xys = storedMetaData{index}.xys;
    xy_dimensions = size(xys);
    totalConditions = xy_dimensions(1);
    bubbletime = storedMetaData{index}.bubbletime;

    for c = 1:totalConditions
        
        % 3. isolate all data from current condition
        conditionData = exptData(exptData(:,23) == c,:); % col 23 = condition values
        
        
        % 4. isolate condition specific raw timestamps for bubble trimming
        timestamps = conditionData(:,2)/3600;       % time in seconds converted to hours
        
        
        % 5. remove data with timestamps prior to and after stabilization
        minTime = 3;  % hr
        maxTime = storedMetaData{index}.bubbletime(c);
        
        times_trim1 = timestamps(timestamps >= minTime);
        conditionData_trim1 = conditionData(timestamps >= minTime,:);

        if maxTime > 0
            conditionData_trim2 = conditionData_trim1(times_trim1 <= maxTime,:);
        else
            conditionData_trim2 = conditionData_trim1;
        end
        clear timestamps times_trim1 conditionData conditionData_trim1 minTime maxTime
        
        
        % 6. isolate size at birth data
        isDrop = conditionData_trim2(:,5);          % col 5 = is drop (1=birth event)
        conditionData_atBirth = conditionData_trim2(isDrop == 1,:);
        clear isDrop conditionData_trim2
        
        % 7. isolate condition specific length, width, V and SA data
        lengthVals = conditionData_atBirth(:,3);        % col 3 = length values
        width = conditionData_atBirth(:,11);            % col 11 = width values
        volumes = conditionData_atBirth(:,12);          % col 12 = calculated va_vals (cubic um)
        surfaceArea = conditionData_atBirth(:,13);      % col 13 = surfaceArea
        
        
        % 8. calculate average and s.e.m. of stabilized data
        SAtoV = surfaceArea./volumes;
        stabilizedData = [lengthVals width SAtoV];
        
        means = mean(stabilizedData);
        counts = length(stabilizedData);
        stds = std(stabilizedData);
        sems = stds./sqrt(counts);
        

        % 9. accumulate data for storage / plotting
        compiledSizeData{c}.mean = means;
        compiledSizeData{c}.std = stds;
        compiledSizeData{c}.count = counts;
        compiledSizeData{c}.sem = sems;
        
        clear means counts stds sems stabilizedData xys xy_dimensions
        clear lengthVals width volumes surfaceArea
    
    end
    
    
    % 10. store data from all conditions into measured data structure
    sizeData{index} = compiledSizeData;     
    clear compiledSizeData
    
    
end

% 10.  save stored data into data structure
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
save('sizeData.mat','sizeData')


%% Part I: SA/V across nutrient concentrations and timescales
%     (B): plot!

% 11.  plot average biovolume production rate over time
clc
clear

cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
load('sizeData.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));
experimentCount = length(dataIndex);

% initialize summary stats for fitting
counter = 0;
summaryMeans = zeros(1,(experimentCount-1)*3 + 6);
summaryConcentrations = zeros(1,(experimentCount-1)*3 + 6);


yAxis = {'length at birth (um)', 'width at birth (um)', 'SA/V at birth (1/um)'};

for e = 1:experimentCount
    
    % identify experiment by date
    index = dataIndex(e);
    date = storedMetaData{index}.date;
    
    % exclude outlier from analysis
    if strcmp(date, '2017-10-31') == 1
        disp(strcat(date,': excluded from analysis'))
        continue
    end
    disp(strcat(date, ': analyze!'))
    
    %     % collect analyzed dates for legend
    %     analyzedDates{e} = date;
    
    % load timescale
    timescale = storedMetaData{index}.timescale;
    
    % isolate biomass prod data for current experiment
    experimentData = sizeData{index};
    
    % isolate concentration data for current experiment
    concentration = storedMetaData{index}.concentrations;
    
    % plot for figures: (1) length, (2) width, (3) SA/volume
    for s = 1:3
        
        figure(s)
        for c = 1:length(concentration)
            
            % if monod experiment
            if ischar(timescale)
                errorbar(log(concentration(c)), experimentData{c}.mean(s), experimentData{c}.sem(s),'o','Color','k','MarkerSize',10);
                hold on
                
                % for stable conditions, accumulate data into summary vector
                counter = counter + 1;
                summaryMeans(counter) = experimentData{c}.mean(s);
                summaryConcentrations(counter) = concentration(c);
                
                % if fluc condition in flux experiment
            elseif timescale == 30 && c == 1
                errorbar(log(concentration(c)), experimentData{c}.mean(s), experimentData{c}.sem(s),'o','Color',[0.25 0.25 0.9],'MarkerSize',10);
                hold on
            elseif timescale == 300 && c == 1
                errorbar(log(concentration(c)), experimentData{c}.mean(s), experimentData{c}.sem(s),'o','Color',[0 .7 .7],'MarkerSize',10);
                hold on
                legend('5 min')
            elseif timescale == 900 && c == 1
                errorbar(log(concentration(c)), experimentData{c}.mean(s), experimentData{c}.sem(s),'o','Color',[1 0.6 0],'MarkerSize',10);
                hold on
            elseif timescale == 3600 && c == 1
                errorbar(log(concentration(c)), experimentData{c}.mean(s), experimentData{c}.sem(s),'o','Color',[1 0.5 0.5],'MarkerSize',10);
                hold on
                
            else % if stable condition in flux experiment
                errorbar(log(concentration(c)), experimentData{c}.mean(s), experimentData{c}.sem(s),'o','Color','k','MarkerSize',10);
                hold on
                
                % for stable conditions, accumulate data into summary vector
                counter = counter + 1;
                summaryMeans(counter) = experimentData{c}.mean(s);
                summaryConcentrations(counter) = concentration(c);
            end
        end
        legend('30 sec','5 min','15 min','60 min','stable')
        ylabel(yAxis{s})
        xlabel('log fold LB dilution')
        title(strcat(yAxis{s},' vs nutrient'))
    end
end


%% Part II: SA/V response to nutrient shifts

%  Strategy:
%
%     0. initialize analysis parameters
%     0. initialize complete meta data

%     1. for all experiments in dataset:
%           2. collect experiment date and exclude outliers (2017-10-31)
%           3. initialize experiment meta data
%           4. load measured data
%           5. gather data for specified condition
%           6. isolate volume and timestamp (corrected for signal lag) data of interest
%           7. remove data not in stabilized region
%           8. remove zeros from mu data (always bounding start and end of tracks)
%           9. bin volumes by 20th of period
%          10. calculate average volume and s.e.m. per timebin
%          11. plot
%    12. repeat for all experiments

%  Part C:
%    13. save volume stats into stored data structure

%  Part D:
%    Same as A-C, except volumes are first separated by cell cycle quarter.
%    See section for complete strategy.

clc
clear

% 0. initialize analysis parameters
binsPerPeriod = 20;

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')

dataIndex = find(~cellfun(@isempty,storedMetaData));
experimentCount = length(dataIndex);

% 1. for all experiments in dataset
exptCounter = 0;
for e = 1:experimentCount
       
    % 2. collect experiment date
    index = dataIndex(e);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    
    % exclude outlier from analysis
    if strcmp(date, '2017-10-31') == 1 || strcmp (timescale, 'monod') == 1
        disp(strcat(date,': excluded from analysis'))
        continue
    end
    disp(strcat(date, ': analyze!'))
    exptCounter = exptCounter + 1;
    datesForLegend{exptCounter} = date;

    
    % 3. initialize experiment meta data
    xys = storedMetaData{index}.xys;
    bubbletime = storedMetaData{index}.bubbletime;
    
    
    % 4. load measured data
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-window5-width1p4-1p7-jiggle-0p5.mat');
    load(filename,'D','D5','M','M_va','T');
    
    
    % 5. gather specified condition data
    condition = 1; % 1 = fluctuating; 3 = ave nutrient condition
    xy_start = storedMetaData{index}.xys(condition,1);
    xy_end = storedMetaData{index}.xys(condition,end);
    flucData = buildDM(D5, M, M_va, T, xy_start, xy_end,e);
    
    
    % 6. remove data not in stabilized region
    minTime = 3;  % hr converted to min
    timestamps = flucData(:,2)/3600;
    
    data_trim1 = flucData(timestamps >= minTime,:);
    time_trim1 = timestamps(timestamps >= minTime);
    
    if bubbletime(condition) == 0
        data_trim2 = data_trim1;
        %time_trim2 = time_trim1;
    else
        maxTime = bubbletime(condition);
        data_trim2 = data_trim1(time_trim1 <= maxTime,:);
        %time_trim2 = time_trim1(time_trim1 <= maxTime);
    end
    clear minTime maxTime bubbletime timestamps time_trim1 time_trim2 data_trim1 
    
    
    % 7. isolate volume, surface area, and timestamp (corrected for signal lag) data of interest
    volumes = data_trim2(:,12);          % col 12 = calculated va_vals (cubic um)
    surfaceArea = data_trim2(:,13);      % col 13 = surfaceArea
    if strcmp(date, '2017-10-10') == 1
        correctedTime = data_trim2(:,2)/3600;
    else
        correctedTime = data_trim2(:,25)/3600; % col 30 = timestamps corrected for signal lag
    end
    clear D D5 M M_va T xy_start xy_end xys
    
    
    % 9. bin SA/V ratio by 20th of period
    SA2V = surfaceArea./volumes;
    timeInSeconds = correctedTime*3600;
    timeInPeriods = timeInSeconds/timescale; % units = sec/sec
    timeInPeriods_floors = floor(timeInPeriods);
    timeInPeriodFraction = timeInPeriods - timeInPeriods_floors;
    assignedBin = timeInPeriodFraction * binsPerPeriod;
    assignedBin = ceil(assignedBin);
    
    binnedRatios = accumarray(assignedBin,SA2V,[],@(x) {x});
    clear timeInSeconds timeInPeriods timeInPeriodFraction assignedBin
    
    
    % 10.  calculate average volume and s.e.m. per timebin
    meanRatio(exptCounter,:) = cellfun(@mean,binnedRatios);
    countRatio(exptCounter,:) = cellfun(@length,binnedRatios);
    stdRatio(exptCounter,:) = cellfun(@std,binnedRatios);
    semRatio(exptCounter,:) = stdRatio(exptCounter,:)./sqrt(countRatio(exptCounter,:));
    
    
    % 11. plot
    if timescale == 30
        color = rgb('FireBrick');
        shapeNum = index-1;
    elseif timescale == 300
        color = rgb('Gold');
        shapeNum = index-4;
    elseif timescale == 900
        color = rgb('MediumSeaGreen');
        shapeNum = index-8;
    elseif timescale == 3600
        color = rgb('MediumSlateBlue');
        shapeNum = index-12;
    end
    
    if shapeNum == 1
        shape = 'x';
    elseif shapeNum == 2
        shape = 'o';
    elseif shapeNum == 3
        shape = 'square';
    else
        shape = '+';
    end

    figure(4)
    errorbar(meanRatio(exptCounter,:),semRatio(exptCounter,:),'Color',color,'Marker',shape)
    hold on
    grid on
    title('surface area:volume across nutrient shifts')
    xlabel('period bin (1/20)')
    ylabel('SA/V, unsynchronized')
    axis([0 21 2.4 3.8])
    legend(datesForLegend)
    

% 12. repeat for all experiments
end



