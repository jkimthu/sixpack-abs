% singlePeriod_parameter

%  Goal: this script bins and plots all growth parameters by time increment (period fraction)
%        output is a parameter v nutrient phase plot with four panels, one
%        for each quarter of the cell cycle. each condition for each experiment
%        is plotted individually, which color indicating data from a specific timescale or
%        stable low/high


%  Last edit: jen, 2018 April 11

%  commit: edit colors for stable environments such that low, ave and high
%          are easily distiguished


%  Strategy:
%
%  Part A:
%     0. initialize analysis parameters
%     0. initialize complete meta data
%     1. for all experiments in dataset
%              2. collect experiment date
%              3. if outlier, exclude from analysis
%              4. else, initialize experiment meta data
%              5. load measured data
%              6. for each condition...
%                       6. gather specified condition data
%                       7. isolate ND2 timestamps
%                       8. remove data not in stabilized region
%                       9. remove data that is not part of a full cell cycle
%                      10. isolate parameters and corrected timestamp (corrected for signal lag)
%                      11. exclude conditions with no data after bubble trimming
%                      12. bin cell cycle fractions into quarters
%                      13. for each quarter period, generate a subplot
%                               14. isolate data from current quarter
%                               15. bin volumes by 20th of period
%                               16. calculate average parameter and s.e.m. per timebin
%                      17. plot, see groups below for figure designations
%             18. repeat for all conditions
%    19. repeat for all experiments


% OK let's go!

%% Group A. size parameters vs. nutrient phase
clc
clear


% 0. initialize analysis parameters
binsPerPeriod = 20;

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')

dataIndex = find(~cellfun(@isempty,storedMetaData));
experimentCount = length(dataIndex);

%%
% 1. for all experiments in dataset
exptCounter = 0;
for e = 1:experimentCount
    
    
    % 2. collect experiment date
    index = dataIndex(e);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    
    
    % 3. if outlier, exclude from analysis
    if strcmp(date, '2017-10-31') == 1 || strcmp (timescale, 'monod') == 1
        disp(strcat(date,': excluded from analysis'))
        continue
    end
    disp(strcat(date, ': analyze!'))
    exptCounter = exptCounter + 1;
    datesForLegend{exptCounter} = date;
    
    
    % 4. else, initialize experiment meta data
    xys = storedMetaData{index}.xys;
    bubbletime = storedMetaData{index}.bubbletime;
    
    
    % 5. load measured data
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-window5-width1p4-1p7-jiggle-0p5.mat');
    load(filename,'D5','M','M_va','T');
    
    
    % 6. gather specified condition data
    for condition = 1:4
        xy_start = storedMetaData{index}.xys(condition,1);
        xy_end = storedMetaData{index}.xys(condition,end);
        conditionData = buildDM(D5, M, M_va, T, xy_start, xy_end,e);
        
        
        % 7. isolate ND2 timestamps
        timestamps = conditionData(:,2)/3600;       % col 2 = ND2 timestamps (convert from sec to hr)
        clear filename xy_start xy_end xys
        
        
        % 8. remove data not in stabilized region
        minTime = 3;  % hr converted to min
        conditionData_trim1 = conditionData(timestamps >= minTime,:);
        timestamps_trim1 = timestamps(timestamps >= minTime);
        
        if bubbletime(condition) == 0
            conditionData_trim2 = conditionData_trim1;
            timestamps_trim2 = timestamps_trim1;
        else
            maxTime = bubbletime(condition);
            conditionData_trim2 = conditionData_trim1(timestamps_trim1 <= maxTime,:);
            timestamps_trim2 = timestamps_trim1(timestamps_trim1 <= maxTime);
        end
        clear minTime maxTime timestamps timestamps_trim1 timestamps_trim2
        
        
        % 9. remove data that is not part of a full cell cycle
        ccFraction = conditionData_trim2(:,9);              % col 9 = cell cycle fraction
        ccFraction_fullOnly = ccFraction(~isnan(ccFraction));
        conditionData_trim3 = conditionData_trim2(~isnan(ccFraction),:);
        clear conditionData conditionData_trim1 conditionData_trim2
        
        
        % 10. isolate parameters and corrected timestamp (corrected for signal lag)
        
        %          a. length
        %          b. width
        %          c. volume
        %          d. surface area
        %          e. SA/V
        
        lengths = conditionData_trim3(:,3);          % col 3 = length
        widths = conditionData_trim3(:,11);          % col 11 = width
        volumes = conditionData_trim3(:,12);         % col 12 = volumes (va)
        surfaceAreas = conditionData_trim3(:,13);     % col 13 = surfaceArea
        SA2V = surfaceAreas./volumes;
        
        if strcmp(date, '2017-10-10') == 1
            correctedTime = conditionData_trim3(:,2)/3600;
        else
            correctedTime = conditionData_trim3(:,25)/3600; % col 25 = timestamps corrected for signal lag
        end
        %clear ccFraction
        
        
        % 11. exclude conditions with no data after bubble trimming
        if isempty(volumes) == 1
            disp(strcat(date,': excluded from analysis, no data'))
            continue
        else
            groupA_growthRate = [lengths widths volumes surfaceAreas SA2V];
            groupA_parameters = {'L','W','V','SA','SA/V'};
        end
        
        
        % 12. bin cell cycle fractions into quarters
        ccInQuarters = ceil(ccFraction_fullOnly*4);
        clear lengths widths volumes surfaceAreas SA2V
        
        
        % 13. generate a subplot for each quarter period
        for q = 1:4
            
            % 14. isolate data from current quarter
            currentQuarter_data = groupA_growthRate(ccInQuarters == q,:);
            currectQuarter_times = correctedTime(ccInQuarters == q);
            
            % 15. bin volumes by 20th of period
            timeInSeconds = currectQuarter_times*3600;
            timeInPeriods = timeInSeconds/timescale; % units = sec/sec
            timeInPeriods_floors = floor(timeInPeriods);
            timeInPeriodFraction = timeInPeriods - timeInPeriods_floors;
            assignedBin = timeInPeriodFraction * binsPerPeriod;
            assignedBin = ceil(assignedBin);
            
            
            % 16. calculate average parameter and s.e.m. per timebin
            for pp = 1:length(groupA_parameters)
                
                binnedData{1,pp} = accumarray(assignedBin,currentQuarter_data(:,pp),[],@(x) {x});
                parameter_mean{1,pp} = cellfun(@mean,binnedData{pp});
                parameter_count{1,pp} = cellfun(@length,binnedData{pp});
                parameter_std{1,pp} = cellfun(@std,binnedData{pp});
                parameter_sem{1,pp} = parameter_std{1,pp}./sqrt(parameter_count{1,pp});
                
            end
            
            if q == 1
                
                q1.parameter_means = parameter_mean;
                q1.parameter_counts = parameter_count;
                q1.parameter_stds = parameter_std;
                q1.parameter_sems = parameter_sem;
                parameter_stats(exptCounter).q1 = q1;
                
            elseif q ==2
                
                q2.parameter_means = parameter_mean;
                q2.parameter_counts = parameter_count;
                q2.parameter_stds = parameter_std;
                q2.parameter_sems = parameter_sem;
                parameter_stats(exptCounter).q2 = q2;
                
            elseif q ==3
                
                q3.parameter_means = parameter_mean;
                q3.parameter_counts = parameter_count;
                q3.parameter_stds = parameter_std;
                q3.parameter_sems = parameter_sem;
                parameter_stats(exptCounter).q3 = q3;
                
            elseif q ==4
                
                q4.parameter_means = parameter_mean;
                q4.parameter_counts = parameter_count;
                q4.parameter_stds = parameter_std;
                q4.parameter_sems = parameter_sem;
                parameter_stats(exptCounter).q4 = q4;
                
            end
            clear timeInSeconds timeInPeriods timeInPeriodFraction assignedBin
            clear parameter_mean parameter_count parameter_std parameter_sem
            
        end
        clear q q1 q2 q3 q4 binnedVolumes currentQuarter_data currectQuarter_times
        clear pp groupA_growthRate
        clear ccInQuarters ccFraction_fullOnly
        
        
        % 17. plot
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
        
        
        if condition == 2
            color = rgb('MidnightBlue');
        elseif condition == 3
            color = rgb('RoyalBlue');
        elseif condition == 4
            color = rgb('Navy');
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
        
        for pp = 1:length(groupA_parameters)
            
            if pp == 1
                ax = [0,21,1.5,10];
            elseif pp == 2
                ax = [0,21,0,4];
            elseif pp == 3
                ax = [0,21,1.5,12];
            elseif pp == 4
                ax = [0,21,6,20];
            elseif pp == 5
                ax = [0,21,0,6];
            end
            
            figure(pp)
            q=1;
            d = parameter_stats(exptCounter).q1.parameter_means{pp};
            err = parameter_stats(exptCounter).q1.parameter_sems{pp};
            subplot(1,4,q)
            errorbar(d,err,'Color',color,'Marker',shape)
            hold on
            grid on
            axis(ax)

            ylabel(groupA_parameters{pp})
            
            q=2;
            d = parameter_stats(exptCounter).q2.parameter_means{pp};
            err = parameter_stats(exptCounter).q2.parameter_sems{pp};
            subplot(1,4,q)
            errorbar(d,err,'Color',color,'Marker',shape)
            hold on
            grid on
            axis(ax)
            
            q=3;
            d = parameter_stats(exptCounter).q3.parameter_means{pp};
            err = parameter_stats(exptCounter).q3.parameter_sems{pp};
            subplot(1,4,q)
            errorbar(d,err,'Color',color,'Marker',shape)
            hold on
            grid on
            axis(ax)
            
            q=4;
            subplot(1,4,q)
            d = parameter_stats(exptCounter).q4.parameter_means{pp};
            err = parameter_stats(exptCounter).q4.parameter_sems{pp};
            errorbar(d,err,'Color',color,'Marker',shape)
            hold on
            grid on
            axis(ax)
            
            title(strcat('err = sem'))
            xlabel('period bin (1/20)')
            legend(datesForLegend)
            
        end
        % 18. repeat for all conditions
    end
    clear conditionData
    
    % 19. repeat for all experiments
end




%%

