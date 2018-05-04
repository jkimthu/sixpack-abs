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


% Last edit: jen, 2018 May 4
% Commit: plot monod using instantaneous dV/dt from all data after 3 hours

% OK let's go!


%% Build compiled data structures with all experiments

% 0. initialize experiment data
clear
clc

cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')

dataIndex = find(~cellfun(@isempty,storedMetaData));
dVdtData = cell(size(storedMetaData));
dVdtData_normalized = cell(size(storedMetaData));

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
        
        
        % 4. isolate data to stabilized regions of growth
        minTime = 3;  % hr
        maxTime = storedMetaData{index}.bubbletime(c);
        timestamps = conditionData(:,2)/3600; % time in seconds converted to hours
        
        times_trim1 = timestamps(timestamps >= minTime);
        conditionData_trim1 = conditionData(timestamps >= minTime,:);
        
        if maxTime > 0
            conditionData_trim2 = conditionData_trim1(times_trim1 <= maxTime,:);
        else
            conditionData_trim2 = conditionData_trim1;
        end

        
        % 5. isolate volume (Va), dVdt and timestamp data from current condition
        volumes = conditionData_trim2(:,12);        % col 12 = calculated va_vals (cubic um)
        timestamps = conditionData_trim2(:,2);      % col 2  = timestamp in seconds
        isDrop = conditionData_trim2(:,5);          % col 5  = isDrop, 1 marks a birth event 
        
        dV_raw = [NaN; diff(volumes)];
        dt = [NaN; diff(timestamps)];
        dt(dt == 0) = NaN;
        dVdt = dV_raw./dt * 3600;
        dVdt_normalizedByVol = dVdt./volumes;
        
        dVdt(isDrop == 1) = NaN;
        dVdt_normalizedByVol(isDrop == 1) = NaN;
        %clear conditionData
        
        
        % 6. calculate average and s.e.m. of stabilized data        
        mean_dVdt = nanmean(dVdt);
        count_dVdt = length(dVdt(~isnan(dVdt)));
        std_dVdt = nanstd(dVdt);
        sem_dVdt = std_dVdt./sqrt(count_dVdt);
        
        mean_dVdt_normalized = nanmean(dVdt_normalizedByVol);
        count_dVdt_normalized = length(dVdt_normalizedByVol(~isnan(dVdt_normalizedByVol)));
        std_dVdt_normalized = nanstd(dVdt_normalizedByVol);
        sem_dVdt_normalized = std_dVdt_normalized./sqrt(count_dVdt_normalized);
        
        
        % 9. accumulate data for storage / plotting        
        compiled_dVdt{c}.mean = mean_dVdt;
        compiled_dVdt{c}.std = std_dVdt;
        compiled_dVdt{c}.count = count_dVdt;
        compiled_dVdt{c}.sem = sem_dVdt;
        
        compiled_dVdt_normalized{c}.mean = mean_dVdt_normalized;
        compiled_dVdt_normalized{c}.std = std_dVdt_normalized;
        compiled_dVdt_normalized{c}.count = count_dVdt_normalized;
        compiled_dVdt_normalized{c}.sem = sem_dVdt_normalized;
        
        clear minTime maxTime 
        clear mean_dVdt std_dVdt count_dVdt sem_dVdt
    
    end
    
    % 10. store data from all conditions into measured data structure        
    dVdtData{index} = compiled_dVdt;
    dVdtData_normalized{index} = compiled_dVdt_normalized;
    
    clear compiled_dVdt
end


%% 11. Save new data into stored data structure
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
save('dVdtData.mat','dVdtData')
save('dVdtData_normalized.mat','dVdtData_normalized')

%% 12. plot average biovolume production rate over time
clc
clear

cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
load('dVdtData.mat')
load('dVdtData_normalized.mat')
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
    experiment_dVdt_data = dVdtData{index};
    experiment_dVdt_norm = dVdtData_normalized{index};
    
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