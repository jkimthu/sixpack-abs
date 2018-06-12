%% figure 34


% Goal: Monod plot of growth rate vs nutrient concentration,
%       where growth rate is calculated in two ways:
%
%       Like figure 21, except width instead of volume
%
%       A. mean of all instantaneous dW/dt in experiment, from ful curves ONLY
%       B. mean of all instantaneous dW/dt, from full curves ONLY,
%          normalized by instantaneous width


% Strategy:
%
%       1. collect instantaneous dW/dt and normalized dW/dt from full curves only, after 3 hours
%       2. calculate mean, std, counts, and sem for each condition of each experiment
%       3. store stats into a structure
%       4. save two structures: one for dW/dt, one for normalized dW/dt
%       5. call data from structures for plotting


% Last edit: jen, 2018 Jun 12

% Commit: plot monod of growth rate in width over LB dilution

% OK let's go!


%% Build compiled data structures with all experiments

% 0. initialize experiment data
clear
clc

cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')

dataIndex = find(~cellfun(@isempty,storedMetaData));
dWdtData_fullOnly = cell(size(storedMetaData));
dWdtData_fullOnly_normalized = cell(size(storedMetaData));

% initialize summary vectors for calculated data
experimentCount = length(dataIndex);

% 1. for each experiment, move to folder and load data
for e = 1:experimentCount
    
    % identify experiment by date
    index = dataIndex(e);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    
    
    % move directory to experiment data
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    
    
    % load data
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
   
    
    % 2. for each condition, calculate mean biovolume production rate per condition
    xys = storedMetaData{index}.xys;
    xy_dimensions = size(xys);
    totalConditions = xy_dimensions(1);
    
    for c = 1:totalConditions
        
        % 3. isolate all data from current condition
        conditionData = exptData(exptData(:,23) == c,:);
        
        % 4. trim data to full cell cycles ONLY
        curveFinder = conditionData(:,6);        % col 6 = curveFinder, ID of full cell cycles
        conditionData_fullOnly = conditionData(curveFinder > 0,:);
        clear curveFinder
        
        % 5. isolate data to stabilized regions of growth
        minTime = 3;  % hr
        maxTime = storedMetaData{index}.bubbletime(c);
        timestamps = conditionData_fullOnly(:,2)/3600; % time in seconds converted to hours
        
        times_trim1 = timestamps(timestamps >= minTime);
        conditionData_trim1 = conditionData_fullOnly(timestamps >= minTime,:);
        
        if maxTime > 0
            conditionData_trim2 = conditionData_trim1(times_trim1 <= maxTime,:);
        else
            conditionData_trim2 = conditionData_trim1;
        end
        clear times_trim1 timestamps minTime maxTime
        
        
        % 6. isolate volume (Va), dVdt and timestamp data from current condition
        widths = conditionData_trim2(:,11);         % col 11 = measured widths
        timestamps = conditionData_trim2(:,2);      % col 2  = timestamp in seconds
        isDrop = conditionData_trim2(:,5);          % col 5  = isDrop, 1 marks a birth event
        curveFinder = conditionData_trim2(:,6);     % col 6  = curve finder (ID of curve in condition)
        
        % calculate mean timestep
        if isempty(widths)
            dWdt = [];
            dWdt_normalizedByWidth = [];
        else
            curveIDs = unique(curveFinder);
            firstFullCurve = curveIDs(2);
            firstFullCurve_timestamps = timestamps(curveFinder == firstFullCurve);
            dt = mean(diff(firstFullCurve_timestamps)); % timestep in seconds
            
            dW_raw = [NaN; diff(widths)];
            dWdt = dW_raw/dt * 3600;                       % units = um/hr
            dWdt_normalizedByWidth = dWdt./widths;
            
            dWdt(isDrop == 1) = NaN;
            dWdt_normalizedByWidth(isDrop == 1) = NaN;
            
        end
        
        % 7. calculate average and s.e.m. of stabilized data        
        mean_dWdt = nanmean(dWdt);
        count_dWdt = length(dWdt(~isnan(dWdt)));
        std_dWdt = nanstd(dWdt);
        sem_dWdt = std_dWdt./sqrt(count_dWdt);
        
        mean_dWdt_normalized = nanmean(dWdt_normalizedByWidth);
        count_dWdt_normalized = length(dWdt_normalizedByWidth(~isnan(dWdt_normalizedByWidth)));
        std_dWdt_normalized = nanstd(dWdt_normalizedByWidth);
        sem_dWdt_normalized = std_dWdt_normalized./sqrt(count_dWdt_normalized);
        
        
        % 8. accumulate data for storage / plotting        
        compiled_dWdt{c}.mean = mean_dWdt;
        compiled_dWdt{c}.std = std_dWdt;
        compiled_dWdt{c}.count = count_dWdt;
        compiled_dWdt{c}.sem = sem_dWdt;
        
        compiled_dWdt_normalized{c}.mean = mean_dWdt_normalized;
        compiled_dWdt_normalized{c}.std = std_dWdt_normalized;
        compiled_dWdt_normalized{c}.count = count_dWdt_normalized;
        compiled_dWdt_normalized{c}.sem = sem_dWdt_normalized;
        
        clear mean_dWdt std_dWdt count_dWdt sem_dWdt
        clear mean_dWdt_normalized std_dWdt_normalized count_dWdt_normalized sem_dWdt_normalized
    
    end
    
    % 10. store data from all conditions into measured data structure        
    dWdtData_fullOnly{index} = compiled_dWdt;
    dWdtData_fullOnly_normalized{index} = compiled_dWdt_normalized;
    
    clear compiled_dWdt compiled_dWdt_normalized
end


%% 11. Save new data into stored data structure
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
save('dWdtData_fullOnly.mat','dWdtData_fullOnly')
save('dWdtData_fullOnly_normalized.mat','dWdtData_fullOnly_normalized')

%% 12. plot average biovolume production rate over time
clc
clear

cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
load('dWdtData_fullOnly.mat')
load('dWdtData_fullOnly_normalized.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));
experimentCount = length(dataIndex);


% initialize summary stats for fitting
counter = 0;
summaryMeans = zeros(1,(experimentCount-1)*3 + 6);
summaryConcentrations = zeros(1,(experimentCount-1)*3 + 6);

% initialize colors
palette = {'FireBrick','Chocolate','ForestGreen','Amethyst','MidnightBlue'};
shapes = {'o','x','square','*'};

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
    
    
    % load timescale
    timescale = storedMetaData{index}.timescale;
    
    % isolate biomass prod data for current experiment
    experiment_dWdt = dWdtData_fullOnly{index};
    experiment_dWdt_norm = dWdtData_fullOnly_normalized{index};
    
    % isolate concentration data for current experiment
    concentration = storedMetaData{index}.concentrations;
    
    
    for c = 1:length(concentration)
        
        % if monod experiment
        if ischar(timescale)
            color = rgb(palette(5));
            xmark = shapes{1};

        % if fluc experiment
        elseif timescale == 30 && c == 1
            color = rgb(palette(1));
            xmark = shapes{1};
        elseif timescale == 300 && c == 1
            color = rgb(palette(2));
            xmark = shapes{2};
        elseif timescale == 900 && c == 1
            color = rgb(palette(3));
            xmark = shapes{3};
        elseif timescale == 3600 && c == 1
            color = rgb(palette(4));
            xmark = shapes{4};
        else
            color = rgb(palette(5));
            xmark = shapes{1};
        end
        
        
        % plot dV/dt data, labeled by stable vs fluc
        figure(1)
        errorbar(log(concentration(c)), experiment_dWdt{c}.mean, experiment_dWdt{c}.sem,'Color',color);
        hold on
        plot(log(concentration(c)), experiment_dWdt{c}.mean,'Marker',xmark,'MarkerSize',10,'Color',color)
        hold on
        ylabel('dW/dt (um/hr)')
        xlabel('log fold LB dilution')
        title(strcat('growth rate in width vs nutrient concentration'))
        
        % plot normalized dV/dt data, labeled by stable vs fluc
        figure(2)
        errorbar(log(concentration(c)), experiment_dWdt_norm{c}.mean, experiment_dWdt_norm{c}.sem,'Color',color);
        hold on
        plot(log(concentration(c)), experiment_dWdt_norm{c}.mean,'Marker',xmark,'MarkerSize',10,'Color',color)
        hold on
        ylabel('(dW/dt)/W (1/hr)')
        xlabel('log fold LB dilution')
        title(strcat('normalized growth rate in width vs nutrient concentration'))
        
    end
     
end