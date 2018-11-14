%% figure 21


% Goal: Monod plot of growth rate vs nutrient concentration


% Strategy:
%
%       1. collect instantaneous growth rates from all curves, after 3 hours
%          note: for lowest monod point, need to allow all curves
%       2. calculate mean, std, counts, and sem for each condition of each experiment
%       3. store stats into a structure
%       4. save a data structure of stats from all experiments
%       5. call data from structures for plotting


% Last edit: jen, 2018 November 14

% Commit: add figure2 which plots monod on linear x axis


% OK let's go!


%% Build compiled data structures with all experiments

% 0. initialize experiment data
clear
clc

cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')

% 0. define growth rate of interest
prompt = 'Enter specific growth rate definition as string (raw / norm / log2 / lognorm): ';
specificGrowthRate = input(prompt);

% 0. initialize cell array for data storage
growthRateData = cell(size(storedMetaData));


%%
% 1. for each experiment, move to folder and load data
exptArray = [2:4,5:7,9,10,11,12,13,14,15,17,18]; % list experiments by index

for e = 14:length(exptArray)     
    
    % identify experiment by date
    index = exptArray(e);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    expType = storedMetaData{index}.experimentType;
    
    % move directory to experiment data
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    
    
    % load data
    if strcmp(date,'2017-11-12') == 1
        filename = 'lb-fluc-2017-11-12-width1p4-jiggle-0p5.mat';
    elseif ischar(timescale) == 0
        filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
    elseif strcmp(date,'2017-09-26') == 1
        filename = 'lb-monod-2017-09-26-jiggle-c12-0p1-c3456-0p5-bigger1p8.mat';
    elseif strcmp(date, '2017-11-09') == 1
        filename = 'lb-control-2017-11-09-width1p4-jiggle-0p5.mat';
    end
    load(filename,'D5','T')
    
    
    % build experiment data matrix
    display(strcat('Experiment (', num2str(e),') of (', num2str(length(exptArray)),')'))
    xy_start = 1;
    xy_end = length(D5);
    exptData = buildDM(D5,T,xy_start,xy_end,index,expType);
    
    clear D5 T filename experimentFolder
   
    
    % 2. for each condition, calculate mean biovolume production rate per condition
    xys = storedMetaData{index}.xys;
    xy_dimensions = size(xys);
    totalConditions = xy_dimensions(1);
    
    
    for c = 1:totalConditions
        
        % 3. isolate all data from current condition
        conditionData = exptData(exptData(:,21) == c,:); % col 21 = condition vals
        
        % 4. trim data to full cell cycles ONLY
        %curveFinder = conditionData(:,5);        % col 5 = curveFinder, ID of full cell cycles
        %conditionData_fullOnly = conditionData(curveFinder > 0,:);
        %clear curveFinder
        
        
        % 5. isolate volume (Va), dVdt and timestamp data from current condition
        volumes = conditionData(:,11);        % col 11 = calculated va_vals (cubic um)
        timestamps_sec = conditionData(:,2);  % col 2  = timestamp in seconds
        isDrop = conditionData(:,4);          % col 4  = isDrop, 1 marks a birth event
        curveFinder = conditionData(:,5);     % col 5  = curve finder (ID of curve in condition)
        trackNum = conditionData(:,20);       % col 20 = track number (not ID from particle tracking)
        
        
        % 6. calculate growth rate
        growthRates = calculateGrowthRate(volumes,timestamps_sec,isDrop,curveFinder,trackNum);
        
        
        
        
        % 7. isolate data to stabilized regions of growth
        %    NOTE: errors (excessive negative growth rates) occur at trimming
        %          point if growth rate calculation occurs AFTER time trim.

        minTime = 3;  % hr
        maxTime = floor(storedMetaData{index}.bubbletime(c)); % limit analysis to whole integer # of periods
        timestamps_hr = conditionData(:,2)/3600; % time in seconds converted to hours
        
        % trim to minumum
        times_trim1 = timestamps_hr(timestamps_hr >= minTime);
        conditionData_trim1 = conditionData(timestamps_hr >= minTime,:);
        growthRates_trim1 = growthRates(timestamps_hr >= minTime,:);
        
        % trim to maximum
        if maxTime > 0
            conditionData_trim2 = conditionData_trim1(times_trim1 <= maxTime,:);
            growthRates_trim2 = growthRates_trim1(times_trim1 <= maxTime,:);
        else
            conditionData_trim2 = conditionData_trim1;
            growthRates_trim2 = growthRates_trim1;
        end
        clear times_trim1 timestamps_hr minTime maxTime
        clear growthRates conditionData_fullOnly conditionData
        
        
        
        % 8. isolate selected specific growth rate
        if strcmp(specificGrowthRate,'raw') == 1
            specificColumn = 1;         % for selecting appropriate column in growthRates
        elseif strcmp(specificGrowthRate,'norm') == 1
            specificColumn = 2;
        elseif strcmp(specificGrowthRate,'log2') == 1
            specificColumn = 3;
        elseif strcmp(specificGrowthRate,'lognorm') == 1
            specificColumn = 4;
        end
        
        growthRt = growthRates_trim2(:,specificColumn);

        
        
        % 9. calculate average and s.e.m. of stabilized data        
        mean_growthRate = nanmean(growthRt);
        count_growthRate = length(growthRt(~isnan(growthRt)));
        std_growthRate = nanstd(growthRt);
        sem_growthRate = std_growthRate./sqrt(count_growthRate);
        
        
        % 10. accumulate data for storage / plotting        
        compiled_growthRate{c}.mean = mean_growthRate;
        compiled_growthRate{c}.std = std_growthRate;
        compiled_growthRate{c}.count = count_growthRate;
        compiled_growthRate{c}.sem = sem_growthRate;
        
        clear mean_growthRate std_growthRate count_growthRate sem_growthRate
    
    end
    
    
    % 11. store data from all conditions into measured data structure        
    %growthRateData_fullOnly{index} = compiled_growthRate;
    growthRateData{index} = compiled_growthRate;
    clear compiled_growthRate 
    
end


%% 11. Save new data into stored data structure
cd('/Users/jen/Documents/StockerLab/Writing/manuscript 1/figure2')
%save(strcat('growthRateData_fullOnly_',specificGrowthRate,'.mat'),'growthRateData_fullOnly','specificGrowthRate')

save(strcat('growthRateData_',specificGrowthRate,'.mat'),'growthRateData','specificGrowthRate')

%% 12. plot growth rate over nutrient concentration
clc
clear

cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
cd('/Users/jen/Documents/StockerLab/Writing/manuscript 1/figure2')
%load('growthRateData_fullOnly_log2.mat')
load('growthRateData_log2.mat')

exptArray = [2:4,5:7,9,10,11,12,13,14,15,17,18]; % list experiments by index


%%

% initialize colors
palette = {'FireBrick','Chocolate','ForestGreen','Amethyst','MidnightBlue'};
%shapes = {'o','x','square','*'};
shapes = {'o','o','o','o'};

for e = 1:length(exptArray)
    
    % identify experiment by date
    index = exptArray(e);
    date = storedMetaData{index}.date;
    
    % exclude outlier from analysis
    if strcmp(date, '2017-10-31') == 1
        disp(strcat(date,': excluded from analysis'))
        continue
    end
    
    if isempty(growthRateData{index}) == 1
    %if isempty(growthRateData_fullOnly{index}) == 1    
        disp(strcat(date, ': empty, skipped!'))
        continue
    else
        disp(strcat(date, ': analyze!'))
    end
    
    % load timescale
    timescale = storedMetaData{index}.timescale;
    
    % isolate growth rate data for current experiment
    experiment_gr_data = growthRateData{index};
    %experiment_gr_data = growthRateData_fullOnly{index};
    
    % isolate concentration data for current experiment
    concentration = storedMetaData{index}.concentrations;
    
    
    % determine axis range for plotting
        if strcmp(specificGrowthRate,'raw') == 1
            xmin = -5;                  % lower limit for plotting x axis
            xmax = 25;                  % upper limit for plotting x axis
        elseif strcmp(specificGrowthRate,'norm') == 1
            xmin = -1;
            xmax = 5;
        elseif strcmp(specificGrowthRate,'log2') == 1
            ymin = 0;
            ymax = 4.1;
        elseif strcmp(specificGrowthRate,'lognorm') == 1
            xmin = -0.5;
            xmax = 1;
        end
        
        
    
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
        
        
        % plot growth rate data, labeled by stable vs fluc
        
        figure(1) % sem
        errorbar(log(concentration(c)), experiment_gr_data{c}.mean, experiment_gr_data{c}.sem,'Color',color);
        hold on
        plot(log(concentration(c)), experiment_gr_data{c}.mean,'Marker',xmark,'MarkerSize',10,'Color',color)
        hold on
        
        figure(2)
        errorbar(concentration(c), experiment_gr_data{c}.mean, experiment_gr_data{c}.sem,'Color',color);
        hold on
        plot(concentration(c), experiment_gr_data{c}.mean,'Marker',xmark,'MarkerSize',10,'Color',color)
        hold on
        
        
    end
    
    figure(1)
    ylabel(strcat('growth rate: (',specificGrowthRate',')'))
    xlabel('log fold LB dilution')
    title(strcat('Population-averaged growth rate vs log LB dilution'))
    axis([-10,1,ymin,ymax])
     
    figure(2)
    ylabel(strcat('growth rate: (',specificGrowthRate',')'))
    xlabel('fraction full LB')
    title(strcat('Population-averaged growth rate vs log LB dilution'))
    axis([0,1,ymin,ymax])
    
    
end

%% plot spread

% initialize colors
palette = {'FireBrick','Chocolate','ForestGreen','Amethyst','MidnightBlue'};
shapes = {'o','x','square','*'};
%shapes = {'o','o','o','o'};



gr_stable = [];
conc_stable = [];
condition_stable = [];

gr_30 =[];
gr_300 = [];
gr_900 = [];
gr_3600 = [];
gr_low = [];
gr_ave = [];
gr_high = [];
gr_notQuiteAve=[];

for e = 1:length(exptArray)
    
    % identify experiment by date
    index = exptArray(e);
    date = storedMetaData{index}.date;
    concentration = storedMetaData{index}.concentrations;
    experimentType = storedMetaData{index}.experimentType;
    
    if isempty(growthRateData{index}) == 1
        %if isempty(growthRateData_fullOnly{index}) == 1
        disp(strcat(date, ': empty, skipped!'))
        continue
    else
        disp(strcat(date, ': analyze!'))
    end
    
    % load timescale
    timescale = storedMetaData{index}.timescale;
    
    for c = 1:length(concentration)
        
        if strcmp(experimentType,'origFluc') == 1
            
            if c==1 %if fluc
                if timescale == 30
                    gr_30 = [gr_30; growthRateData{index}{c}.mean];
                elseif timescale == 300
                    gr_300 = [gr_300; growthRateData{index}{c}.mean];
                elseif timescale == 900
                    gr_900 = [gr_900; growthRateData{index}{c}.mean];
                elseif timescale == 3600
                    gr_3600 = [gr_3600; growthRateData{index}{c}.mean];
                end
            elseif c == 2
                gr_low = [gr_low; growthRateData{index}{c}.mean];
            elseif c == 3
                gr_ave = [gr_ave; growthRateData{index}{c}.mean];
            elseif c == 4
                gr_high = [gr_high; growthRateData{index}{c}.mean];
            end
            
        end
        
        if strcmp(date,'2017-09-26') == 1
            if c==1
                gr_stable{8} = growthRateData{index}{c}.mean;
            elseif c==2
                gr_stable{7} = growthRateData{index}{c}.mean;
            elseif c == 3
                gr_stable{6} = growthRateData{index}{c}.mean;
            elseif c == 4
                gr_notQuiteAve = growthRateData{index}{c}.mean;
            elseif c == 5
                gr_low = [gr_low; growthRateData{index}{c}.mean];
            elseif c == 6
                gr_stable{1} = growthRateData{index}{c}.mean;
            end
        end
        
        if strcmp(date,'2017-11-09') == 1
            gr_notQuiteAve = [gr_notQuiteAve; growthRateData{index}{c}.mean]
        end
        
        %         % isolate growth rate data and condition for current experiment
        %             gr_stable = [gr_stable; growthRateData{index}{c}.mean];
        %             gr_high = [gr_high; growthRateData{index}{c}.mean];
        %             condition_stable = [condition_stable;c];
        %     % isolate concentration data for current experiment
        %     conc_stable = [conc_stable; concentration];
        
        
    end
end
%%

gr_stable{2} = gr_low;
gr_stable{3} = gr_notQuiteAve;
gr_stable{4} = gr_ave;
gr_stable{5} = gr_high;

gr_fluc_30 = cell(1,8);
gr_fluc_30{4} = gr_30;

gr_fluc_300 = cell(1,8);
gr_fluc_300{4} = gr_300;

gr_fluc_900 = cell(1,8);
gr_fluc_900{4} = gr_900;

gr_fluc_3600 = cell(1,8);
gr_fluc_3600{4} = gr_3600;


plotSpread(gr_stable, ...
    'xNames', {'1/10000', '1/1000', '1/100','1/95','1/50','1/32','1/8','Full'}, ...
    'distributionMarkers', {'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o'});
hold on
plotSpread(gr_fluc_30, ...
    'xNames', {'1/10000', '1/1000', '1/100','1/95','1/50','1/32','1/8','Full'}, ...
    'distributionMarkers', {'o', 'o', 'o', '*', 'o', 'o', 'o', 'o'});
hold on
plotSpread(gr_fluc_300, ...
    'xNames', {'1/10000', '1/1000', '1/100','1/95','1/50','1/32','1/8','Full'}, ...
    'distributionMarkers', {'o', 'o', 'o', 'sq', 'o', 'o', 'o', 'o'});
hold on
plotSpread(gr_fluc_900, ...
    'xNames', {'1/10000', '1/1000', '1/100','1/95','1/50','1/32','1/8','Full'}, ...
    'distributionMarkers', {'o', 'o', 'o', 'd', 'o', 'o', 'o', 'o'});
hold on
plotSpread(gr_fluc_3600, ...
    'xNames', {'1/10000', '1/1000', '1/100','1/95','1/50','1/32','1/8','Full'}, ...
    'distributionMarkers', {'o', 'o', 'o', 'p', 'o', 'o', 'o', 'o'});

%%

title(strcat('Population-averaged growth rate vs LB dilution'))
ylabel(strcat('growth rate: (',specificGrowthRate',')'))
xlabel('LB dilution')