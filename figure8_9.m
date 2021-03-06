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




%  Last edit: jen, 2018 May 14
%  Commit: replace dvdt function with in-script calculation of dV/dt,
%          replace method of id-ing unique cell cycles (use curveID instead of unique interdiv time) 

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
        
        
        % 9. isolate volume, interdiv time and other datas
        interdivTime = conditionData_trim(:,8)/60;     % col 8  = curve duration (sec converted to min)
        volumes = conditionData_trim(:,12);             % col 12  = volume (Va)
        timestamps = conditionData_trim(:,2);      % col 2  = timestamp in seconds
        isDrop = conditionData_trim(:,5);          % col 5  = isDrop, 1 marks a birth event
        curveFinder = conditionData_trim(:,6);     % col 6  = curve finder (ID of curve in condition)
        
        
        % 10. calculate mean timestep and dVdt
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
        
        
        % 11. for each unique cell cycle...
        %           i. store birth and division size
        %          ii. remove NaN values from dV/dt vector
        growthRate = nan(length(curveIDs),1);
        V_birth = nan(length(curveIDs),1);
        interDivs = nan(length(curveIDs),1);
        
        for cc = 1:length(curveIDs)
            
            currentVolumes = volumes(curveFinder == curveIDs(cc));
            currentGrowthRates = dVdt(curveFinder == curveIDs(cc));
            currentInterdivs = interdivTime(curveFinder == curveIDs(cc));

            V_birth(cc,1) = currentVolumes(1);
            growthRate(cc,1) = nanmean(currentGrowthRates);
            interDivs(cc,1) = currentInterdivs(1);
            
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
        interdiv_means = mean(interDivs);
        
        V_birth_stds = std(V_birth);
        growthRate_stds = std(growthRate);
        interdiv_stds = std(interDivs);
   
        
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
