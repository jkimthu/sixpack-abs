%% figures 10 and 11


%  Goal: determine how well birth size and interdivision time are predicted
%        by single-cell average dV/dt, plotting single experiment data:

%        (i) birth size vs dV/dt
%       (ii) inter-division time vs dV/dt


%  Strategy: 

%       0. initialize experiment data
%       1. analyze cell cycles born later than 3 hrs into experiment
%       2. store birth size of each cycle, inter-division time and
%          mean dV/dt of each cycle
%       3. plot mean birth size vs mean dV/dt
%       4. plot mean interdiv time vs dV/dt




%  Last edit: jen, 2018 May 17
%  Commit: update method of dV/dt calculation.

%  OK let's go!

%% initialize

clc
clear

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')

dataIndex = find(~cellfun(@isempty,storedMetaData));
experimentCount = length(dataIndex);
exptsToInclude = [6,7,10:13];
%%
for i = 1:length(exptsToInclude)
    
    % 2. collect experiment date
    e = exptsToInclude(i);
    index = dataIndex(e);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    bubbletime = storedMetaData{index}.bubbletime;
    
    
    % 2. load measured data
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-window5-width1p4-1p7-jiggle-0p5.mat');
    load(filename,'D5','M','M_va','T');
    
    
    % 3. compile experiment data matrix
    xy_start = min(min(storedMetaData{index}.xys));
    xy_end = max(max(storedMetaData{index}.xys));
    exptData = buildDM(D5, M, M_va, T, xy_start, xy_end,e);
    clear D5 M M_va T xy_start xy_end e
    
    
    % 4. initialize colors for plotting
    palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
    environments = {'fluc','low','ave','high'};
    
    for condition = 1:length(environments)
        
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
        
        
        % 9. isolate volume (Va) and timestamp data
        volumes = conditionData_trim(:,12);        % col 12 = calculated va_vals (cubic um)
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
        
        clear isDrop
        
        
        % 11. for each unique cell cycle...
        %           i. determine if cell cycle is starts after experimental hr 3
        %          ii. if not, continue
        %         iii. if yes, store birth and division size data
        %          iv. remove NaN values from dV/dt vector
        
        curveCounter = 0;
        rawTime = timestamps/3600;      % convert to hours
        interdivs = conditionData_trim(:,8)/60;     % col 8 = curve duration (sec converted to min)
        
        for cc = 1:length(curveIDs)
            
            currentTimes = rawTime(curveFinder == curveIDs(cc));
            if currentTimes(1) < 3
                %disp(strcat(num2str(currentTimes(1)),': toss!'))
                continue
            else
                %disp(strcat(num2str(currentTimes(1)),': analyze!'))
                curveCounter = curveCounter + 1;
            end
            
            currentVolumes = volumes(curveFinder == curveIDs(cc));
            currentGrowthRates = dVdt(curveFinder == curveIDs(cc));
            currentInterDivs = interdivs(curveFinder == curveIDs(cc));
            
            V_birth(curveCounter,1) = currentVolumes(1);
            growthRate(curveCounter,1) = nanmean(currentGrowthRates);
            interdivisionz(curveCounter,1) = currentInterDivs(1);
            
        end
        
        
        % 12. plot
        color = rgb(palette(condition));
        
        % (i) birth size vs. dVdt
        figure(i)
        subplot(2,2,condition)
        plot(growthRate,V_birth,'o','Color',color)
        hold on
        legend(environments(condition))
        ylabel('birth size (cubic um)')
        xlabel('dV/dt (cubic um/hr)')
        title(strcat(date,': (',num2str(curveCounter),') of (',num2str(cc),') total cycles'))
        axis([-10 50 0 20])
        
        % (ii) inter-division time vs. dVdt
        figure(i+10)
        subplot(2,2,condition)
        plot(growthRate,interdivisionz,'o','Color',color)
        hold on
        legend(environments(condition))
        ylabel('inter-division time (min)')
        xlabel('dV/dt (cubic um/hr)')
        title(strcat(date,': (',num2str(curveCounter),') of (',num2str(cc),') total cycles'))
        axis([-10 50 0 120])
        
        
        
    end
end



