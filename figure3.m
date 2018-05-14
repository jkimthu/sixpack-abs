%% figure 3


%  Goal: how does dV/dt relate to birth size?

%        plot population-averaged dV/dt vs birth size.
%        population averages are calculated by accumulating vectors of
%        birth size and mean dV/dt per cell cycle, and simply taking the
%        mean of both.


%  Strategy: 
%
%       a) initialize experimental data
%       b) identify complete cell cycles within each condition 
%       c) compile birth volumes, inter-division times, and dV/dt
%       d) calculate mean and plot condition data



%  Last edit: jen, 2018 May 14
%  Commit: replace dvdt function with in-script calculation


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
    xys = storedMetaData{index}.xys;
    
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
    
    
    % 3. initialize colors for plotting
    palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
    shapes = {'o','x','square','*'};
    
    
    for condition = 1:length(palette)
        
        
        % 4. compile experiment data matrix
        xy_start = min(xys(condition,:));
        xy_end = max(xys(condition,:));
        conditionData = buildDM(D5, M, M_va, T, xy_start, xy_end,e);
        clear xy_start xy_end
        
       
        % 5. trim data to full cell cycles ONLY
        curveFinder = conditionData(:,6);        % col 6 = curveFinder, ID of full cell cycles
        conditionData_fullOnly = conditionData(curveFinder > 0,:);
        clear curveFinder
        
        
        % 6. isolate interdivision time
        interdivTime = conditionData_fullOnly(:,8)/60;     % col 8 = curve duration (sec converted to min)
        
        
        % 7. trim data to include only full cell cycles longer than 10 min
        conditionData_trim = conditionData_fullOnly(interdivTime > 10,:);
        clear interdivTime
        
        
        % 8. isolate volume, interdiv time and other datas
        interdivTime = conditionData_trim(:,8)/60;     % col 8  = curve duration (sec converted to min)
        volumes = conditionData_trim(:,12);             % col 12  = volume (Va)
        timestamps = conditionData_trim(:,2);      % col 2  = timestamp in seconds
        isDrop = conditionData_trim(:,5);          % col 5  = isDrop, 1 marks a birth event
        curveFinder = conditionData_trim(:,6);     % col 6  = curve finder (ID of curve in condition)
        
        
        % 9. calculate mean timestep and dVdt
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
        
     
        % 10. for each unique cell cycle, calculate C+D period (min)
        growthRate = nan(length(curveIDs),1);
        V_birth = nan(length(curveIDs),1);
        
        for cc = 1:length(curveIDs)
            
            currentVolumes = volumes(curveFinder == curveIDs(cc));
            currentGrowthRates = dVdt(curveFinder == curveIDs(cc));

            V_birth(cc,1) = currentVolumes(1);
            growthRate(cc,1) = nanmean(currentGrowthRates);
            
        end
        
        
        % 11. plot
        color = rgb(palette(condition));
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
        
        V_birth_stds = std(V_birth);
        growthRate_stds = std(growthRate);
   
        
        % (i) dV/dt vs. birth size
        figure(3)
        errorbar(V_birth_means,growthRate_means,growthRate_stds,'Color',color)
        hold on
        plot(V_birth_means,growthRate_means,'Marker',xmark,'Color',color)
        hold on
        xlabel('birth size (cubic um)')
        ylabel('dV/dt (cubic um/hr)')
        title('population averages from all experiments')
        axis([1 6.5 0 23])
        
        
    end
    clear D5 M M_va T
    
    
end
