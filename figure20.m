%% figure 20


% Goal: Monod plot of growth rate vs nutrient concentration,
%       where growth rate is calculated in two ways:
%
%       A. mean of all instantaneous dV/dt in experiment
%       B. mean of all instantaneous dV/dt, normalized by instantaneous volume


% Strategy:
%
%       1. collect dV/dt and normalized dV/dt from all experiment data, after 3 hours
%       2. calculate mean, std, counts, and sem for each condition of each experiment
%       3. store stats into a structure
%       4. save two structures: one for dV/dt, one for normalized dV/dt
%       5. call data from structures for plotting


% Last edit: jen, 2018 May 15
% Commit: bump 3hr trim BELOW dV/dt calculation, only plots dV/dt not
%         normalized version by initial volume.

% OK let's go!


%% Build compiled data structures with all experiments

% 0. initialize experiment data
clear
clc

cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')

dataIndex = find(~cellfun(@isempty,storedMetaData));
dVdtData_newdVdt = cell(size(storedMetaData));
dVdtData_normalized_newdVdt = cell(size(storedMetaData));

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
        
        
        % 4. isolate condition data to those with full cell cycles
        curveIDs = conditionData(:,6);           % col 6 = curve ID
        conditionData_fullOnly = conditionData(curveIDs > 0,:);
        clear curveIDs
        
        % 5. isolate volume (Va), dVdt and timestamp data from current condition
        volumes = conditionData_trim2(:,12);        % col 12 = calculated va_vals (cubic um)
        timestamps = conditionData_trim2(:,2);      % col 2  = timestamp in seconds
        isDrop = conditionData_trim2(:,5);          % col 5  = isDrop, 1 marks a birth event 
        curveFinder = conditionData_trim2(:,6);     % col 6  = curve finder (ID of curve in condition)
        
        
        % 6. calculate mean timestep and dV/dt
        curveIDs = unique(curveFinder);
        firstFullCurve = curveIDs(2);
        if length(firstFullCurve) > 1
            firstFullCurve_timestamps = timestamps(curveFinder == firstFullCurve);
        else
            firstFullCurve = curveIDs(3);
            firstFullCurve_timestamps = timestamps(curveFinder == firstFullCurve);
        end
        dt = mean(diff(firstFullCurve_timestamps)); % timestep in seconds
        
        dV_raw = [NaN; diff(volumes)];
        dVdt = dV_raw/dt * 3600;                    % final units = cubic um/sec
        dVdt(isDrop == 1) = NaN;

        
        % 7. isolate data to stabilized regions of growth
        minTime = 3;  % hr
        maxTime = bubbletime(condition);
        
        times_trim1 = timestamps(timestamps >= minTime);
        conditionData_trim1 = conditionData_fullOnly(timestamps >= minTime,:);
        dVdt_trim1 = dVdt(timestamps >= minTime,:);
        
        if maxTime > 0
            conditionData_fullOnly = conditionData_trim1(times_trim1 <= maxTime,:);
            dVdt_trim2 = dVdt_trim1(times_trim1 <= maxTime,:);
        else
            conditionData_fullOnly = conditionData_trim1;
            dVdt_trim2 = dVdt_trim1;
        end
        clear times_trim1 dVdt_trim1 timestamps minTime maxTime
        
        
        % 8. calculate average and s.e.m. of stabilized data        
        mean_dVdt = nanmean(dVdt_trim2);
        count_dVdt = length(dVdt_trim2(~isnan(dVdt_trim2)));
        std_dVdt = nanstd(dVdt_trim2);
        sem_dVdt = std_dVdt./sqrt(count_dVdt);
        
        
        % 9. accumulate data for storage / plotting        
        compiled_dVdt{c}.mean = mean_dVdt;
        compiled_dVdt{c}.std = std_dVdt;
        compiled_dVdt{c}.count = count_dVdt;
        compiled_dVdt{c}.sem = sem_dVdt;
        
        clear minTime maxTime 
        clear mean_dVdt std_dVdt count_dVdt sem_dVdt
    
    end
    
    % 10. store data from all conditions into measured data structure        
    dVdtData_newdVdt_fullONLY{index} = compiled_dVdt;

    
    clear compiled_dVdt
end


%% 11. Save new data into stored data structure
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
save('dVdtData_newdVdt.mat','dVdtData_newdVdt')
save('dVdtData_normalized_newdVdt.mat','dVdtData_normalized_newdVdt')

%% 12. plot average biovolume production rate over time
clc
clear

cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
load('dVdtData_newdVdt.mat')
load('dVdtData_normalized_newdVdt.mat')
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
    experiment_dVdt_data = dVdtData_newdVdt{index};
    experiment_dVdt_norm = dVdtData_normalized_newdVdt{index};
    
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
        errorbar(log(concentration(c)), experiment_dVdt_data{c}.mean, experiment_dVdt_data{c}.sem,'Color',color);
        hold on
        plot(log(concentration(c)), experiment_dVdt_data{c}.mean,'Marker',xmark,'MarkerSize',10,'Color',color)
        hold on
        ylabel('dV/dt (cubic um/hr)')
        xlabel('log fold LB dilution')
        title(strcat('Population-averaged dV/dt vs log LB dilution'))
        
        % plot normalized dV/dt data, labeled by stable vs fluc
        figure(2)
        errorbar(log(concentration(c)), experiment_dVdt_norm{c}.mean, experiment_dVdt_norm{c}.sem,'Color',color);
        hold on
        plot(log(concentration(c)), experiment_dVdt_norm{c}.mean,'Marker',xmark,'MarkerSize',10,'Color',color)
        hold on
        ylabel('(dV/dt)/V (1/hr)')
        xlabel('log fold LB dilution')
        title(strcat('Population-averaged volume-normalized dV/dt vs log LB dilution'))
        
    end
     
end