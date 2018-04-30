%% figures 8 and 9


%  Goal: determine how well birth size and interdivision time are predicted
%        by population-averaged dV/dt, plotting:

%        (i) birth size vs dV/dt
%       (ii) inter-division time vs dV/dt


%  Strategy: 
%       0. initialize experiment data
%       1. for all cell cycles in all experiments,
%       2. determine birth size of each cycle, inter-division time and
%          mean dV/dt of each cycle
%       3. plot population-averaged birth size vs dV/dt
%       4. plot population-averaged interdiv time vs mean dV/dt




%  Last edit: jen, 2018 April 30
%  Commit: re-run after editing dvdt calculations to include 2017-10-10 in figure

%  OK let's go!


%% initialize

clc
clear

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')

dataIndex = find(~cellfun(@isempty,storedMetaData));
experimentCount = length(dataIndex);


% 1. for all experiments in dataset
exptCounter = 0;
for e = 1:experimentCount
    
    
    % 1. collect experiment meta data
    index = dataIndex(e);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    bubbletime = storedMetaData{index}.bubbletime;
    
    % exclude outliers from analysis (2017-10-31 and monod experiments)
    if strcmp(date, '2017-10-31') == 1 || strcmp (timescale, 'monod') == 1
        disp(strcat(date,': excluded from analysis'))
        continue
    end
    disp(strcat(date, ': analyze!'))
    exptCounter = exptCounter + 1;
    
    
    % 2. load measured data
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-window5-width1p4-1p7-jiggle-0p5.mat');
    load(filename,'D5','M','M_va','T');
    
    
    % 3. compile experiment data matrix
    xy_start = min(min(storedMetaData{index}.xys));
    xy_end = max(max(storedMetaData{index}.xys));
    exptData = buildDM(D5, M, M_va, T, xy_start, xy_end,e);
    %exptData = buildDM(D5, M, M_va, T, 1, 10,e);
    clear D5 M M_va T xy_start xy_end
    
    
    % 5. initialize colors for plotting
    palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
    shapes = {'o','x','square','*'};
    
    for condition = 1:length(bubbletime)
        
        
        % 5. isolate condition specific data
        conditionData = exptData(exptData(:,23) == condition,:);  % col 23 = cond vals
        
        
        % 6. trim data to full cell cycles ONLY
        curveFinder = conditionData(:,6);        % col 6 = curveFinder, ID of full cell cycles
        conditionData_fullOnly = conditionData(curveFinder > 0,:);
        clear curveFinder
        
        
        % 7. isolate interdivision time
        interdivTime = conditionData_fullOnly(:,8)/60;     % col 8 = curve duration (sec converted to min)
        
        
        % 8. trim data to include only full cell cycles longer than 10 min
        conditionData_trim = conditionData_fullOnly(interdivTime > 10,:);
        clear interdivTime
        
        
        % 9. isolate volume, interdiv time and dV/dt data
        interdivTime = conditionData_trim(:,8)/60;     % col 8  = curve duration (sec converted to min)
        volume = conditionData_trim(:,12);             % col 12  = volume (Va)
        
        dvdt_data = dvdt(conditionData_trim, timescale, date);
        dVdt = dvdt_data(:,1);
        
        
        % 10. identify unique interdivision times (unique cell cycles)
        unique_interdivs = unique(interdivTime);
        
        
        % 11. for each unique cell cycle...
        %           i. store birth and division size
        %          ii. remove NaN values from dV/dt vector
        growthRate = nan(length(unique_interdivs),1);
        V_birth = nan(length(unique_interdivs),1);
        
        for cc = 1:length(unique_interdivs)
            
            currentVolumes = volume(interdivTime == unique_interdivs(cc));
            currentGrowthRates = dVdt(interdivTime == unique_interdivs(cc));

            V_birth(cc,1) = currentVolumes(1);
            growthRate(cc,1) = nanmean(currentGrowthRates);
            
        end
        
        % 12. plot
        color = rgb(palette(condition));
        %color = rgb(palette(e-3));
        
        if condition == 1 && timescale == 300
            xmark = shapes{2};
        elseif condition == 1 && timescale == 900
            xmark = shapes{3};
        elseif condition == 1 && timescale == 3600
            xmark = shapes{4};
        else
            xmark = shapes{1};
        end
        
        V_birth_means = mean(V_birth);
        growthRate_means = mean(growthRate);
        interdiv_means = mean(unique_interdivs);
        
        V_birth_stds = std(V_birth);
        growthRate_stds = std(growthRate);
        interdiv_stds = std(unique_interdivs);
   
        
        % (i) birth size vs. dVdt
        figure(8)
        errorbar(growthRate_means,V_birth_means,V_birth_stds,'Color',color)
        hold on
        plot(growthRate_means,V_birth_means,'Marker',xmark,'Color',color)
        hold on
        ylabel('birth size (cubic um)')
        xlabel('dV/dt (cubic um/hr)')
        title('population averages from all experiments')
        axis([0 19 0 8])
        
        % (ii) inter-division time vs. dVdt
        figure(9)
        errorbar(growthRate_means,interdiv_means,interdiv_stds,'Color',color)
        hold on
        plot(growthRate_means,interdiv_means,'Marker',xmark,'Color',color)
        hold on
        ylabel('inter-division time (min)')
        xlabel('dV/dt (cubic um/hr)')
        title('population averages from all experiments')
        axis([0 19 0 80])

    
    end
      
    
end
