%% whichGrowthModel


%  Goal: determine how well birth size and interdivision time are predicted
%        by single-cell average dV/dt, plotting 2018-02-01 data:

%        (i) birth size vs dV/dt
%       (ii) inter-division time vs dV/dt


%  Strategy: 
%       0. initialize experiment data
%       1. for all cell cycles in all experiments,
%       2. determine birth size of each cycle, inter-division time and
%          mean dV/dt of each cycle
%       3. plot mean birth size vs mean dV/dt
%       4. plot mean interdiv time vs dV/dt




%  Last edit: jen, 2018 April 30

%  Commit: 

%  OK let's go!


%% initialize

clc
clear

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')

dataIndex = find(~cellfun(@isempty,storedMetaData));
experimentCount = length(dataIndex);


% 1. collect 2018-02-01 meta data
e = 14;
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
    
    
    % 9. isolate volume, interdiv time, raw timestamp and dV/dt data
    interdivTime = conditionData_trim(:,8)/60;     % col 8  = curve duration (sec converted to min)
    volume = conditionData_trim(:,12);             % col 12  = volume (Va)
    rawTime = conditionData_trim(:,2)/3600;        % col 2  = timestamps from ND2 (sec converted to hr)
    
    dvdt_data = dvdt(conditionData_trim, timescale);
    dVdt = dvdt_data(:,1);
    
    
    % 10. identify unique interdivision times (unique cell cycles)
    unique_interdivs = unique(interdivTime);
    
    
    % 11. for each unique cell cycle...
    %           i. determine if cell cycle is starts after experimental hr 3
    %          ii. if not, continue
    %         iii. if yes, store birth and division size data 
    %          iv. remove NaN values from dV/dt vector
          
    curveCounter = 0;
    for cc = 1:length(unique_interdivs)
        
        currentTimes = rawTime(interdivTime == unique_interdivs(cc));
        if currentTimes(1) < 3
            %disp(strcat(num2str(currentTimes(1)),': toss!'))
            continue
        else
            %disp(strcat(num2str(currentTimes(1)),': analyze!'))
            curveCounter = curveCounter + 1;
        end
        
        currentVolumes = volume(interdivTime == unique_interdivs(cc));
        currentGrowthRates = dVdt(interdivTime == unique_interdivs(cc));
        
        V_birth(curveCounter,1) = currentVolumes(1);
        growthRate(curveCounter,1) = nanmean(currentGrowthRates);
        interdivisionz(curveCounter,1) = unique_interdivs(cc);
        
    end
    
    
    % 12. plot
    color = rgb(palette(condition));
    
    % (i) birth size vs. dVdt
    figure(10)
    subplot(2,2,condition)
    plot(growthRate,V_birth,'o','Color',color)
    hold on
    legend(environments(condition))
    ylabel('birth size (cubic um)')
    xlabel('dV/dt (cubic um/hr)')
    title(strcat(date,': (',num2str(curveCounter),') of (',num2str(cc),') total cycles'))
    axis([-10 50 0 20])
    
    % (ii) inter-division time vs. dVdt
    figure(11)
    subplot(2,2,condition)
    plot(growthRate,interdivisionz,'o','Color',color)
    hold on
    legend(environments(condition))
    ylabel('inter-division time (min)')
    xlabel('dV/dt (cubic um/hr)')
    title(strcat(date,': (',num2str(curveCounter),') of (',num2str(cc),') total cycles'))
    axis([-10 50 0 120])
    
  
    
    
end



