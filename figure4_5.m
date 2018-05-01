%% figures 4 and 5


%  Goal: determine whether growth in low, ave, high and fluc environments
%        are characterizable as adder, sizer, or timer

%        plot (1) single cell division size vs birth size
%        plot (2) single cell inter-division time vs birth time


%  Strategy: 
%
%       a) initialize experimental data
%       b) identify complete cell cycles within each condition 
%       c) compile birth volumes, division volumes and inter-division times
%       d) plot condition data, for each cycle
%



%  Last edit: jen, 2018 May 1
%  Commit: adder or sizer? single-cell division volume vs birth volume,
%          and inter-div time vs. birth volume for 2018-02-01 data


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

e = 14;


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
environment = {'fluc','low','ave','high'};

for condition = 1:length(environment)
    
    
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
    
    
    % 9. isolate volume and interdiv time
    interdivTime = conditionData_trim(:,8)/60;     % col 8  = curve duration (sec converted to min)
    volume = conditionData_trim(:,12);             % col 12  = volume (Va)
    
    
    % 11. identify unique interdivision times (unique cell cycles)
    unique_interdivs = unique(interdivTime);
    
    
    % 12. for each unique cell cycle, calculate C+D period (min)
    V_division = nan(length(unique_interdivs),1);
    V_birth = nan(length(unique_interdivs),1);
    
    for cc = 1:length(unique_interdivs)
        
        currentVolumes = volume(interdivTime == unique_interdivs(cc));
        V_division(cc,1) = currentVolumes(end);
        V_birth(cc,1) = currentVolumes(1);
        
    end
    
    
    % 13. plot
    color = rgb(palette(condition));
    
    % (i) division size vs. birth size
    figure(4)
    subplot(2,2,condition)
    plot(V_birth,V_division,'o','Color',color)
    hold on
    legend(environment(condition))
    title(date)
    xlabel('birth size (cubic um)')
    ylabel('division size (cubic um)')
    axis([0 20 0 25])
    
    % (ii) inter-division time vs. birth size
    figure(5)
    subplot(2,2,condition)
    plot(V_birth,unique_interdivs,'o','Color',color)
    hold on
    legend(environment(condition))
    title(date)
    xlabel('birth size (cubic um)')
    ylabel('inter-division time (min)')
    axis([0 15 0 130])
    
 
end

