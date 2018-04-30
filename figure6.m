%% figure6


%  Goal: determine whether pre-steady-state cells in stable environments
%        exhibit division sizes that are twice are large as their birth
%        sizes.

%        apparent from figure 1, population-averaged division size vs birth
%        size, all conditions including fluctuating environments follow
%        V_div = 2 * V_birth. How expected is this from cells adapting to
%        new environments?



%  Strategy: 
%
%       0. initialize experiment data
%       1. for experiment 2018-02-01, find all full cell cycles in
%          low, ave and high conditions PRIOR TO STEADY-STATE
%               - from figure 7, dV/dt vs time for 2018-02-01 data, we see
%                 that all stable environments are not yet at steady-state
%                 between hours 0 - 3.
%
%       2. determine division size and birth size of each cycle
%       3. plot division size vs birth size



%  Last edit: jen, 2018 April 29

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


% 1. for all experiments in dataset
e = 14;  % 14 = index for experiment 2018-02-01


% 1. collect experiment meta data
index = dataIndex(e);
date = storedMetaData{index}.date;
timescale = storedMetaData{index}.timescale;


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

for condition = 1:length(palette)
    
    
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
    
    
    % 9. isolate volume, interdiv time, and raw timestamp
    interdivTime = conditionData_trim(:,8)/60;     % col 8  = curve duration (sec converted to min)
    volume = conditionData_trim(:,12);             % col 12 = volume (Va)
    rawTime = conditionData_trim(:,2)/3600;        % col 2  = timestamps from ND2 (sec converted to hr)
    
    
    % 10. identify unique interdivision times (unique cell cycles)
    unique_interdivs = unique(interdivTime);
    
    
    % 11. for each unique cell cycle...
    %           i. determine if cell cycle is within 1st 3 experimental hrs
    %          ii. if not, continue
    %         iii. if yes, store birth and division size data 
    
    curveCounter = 0;
    for cc = 1:length(unique_interdivs)
        
        currentTimes = rawTime(interdivTime == unique_interdivs(cc));
        if currentTimes(end) > 3.2
            %disp(strcat(num2str(currentTimes(end)),': toss!'))
            continue
        else
            %disp(strcat(num2str(currentTimes(end)),': analyze!'))
            curveCounter = curveCounter + 1;
        end
        
        currentVolumes = volume(interdivTime == unique_interdivs(cc));
        V_division(curveCounter,1) = currentVolumes(end);
        V_birth(curveCounter,1) = currentVolumes(1);
        
    end
    
    
    % 12. plot
    color = rgb(palette(condition));
    
    % prior to steady-state division size vs. birth size
    figure(6)
    subplot(2,2,condition)
    plot(V_birth,V_division,'o','Color',color)
    hold on
    legend(environments(condition))
    title(strcat(date,': (',num2str(curveCounter),') of (',num2str(cc),') total cycles'))
    xlabel('birth size (cubic um)')
    ylabel('division size (cubic um)')
    axis([0 25 0 25])
    
    
end


