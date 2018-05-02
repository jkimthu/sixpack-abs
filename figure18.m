%% figure 18


%  Goal: determine whether growth in low, ave, high and fluc environments
%        are characterizable as adder, sizer, or timer

%        plot (1) single cell added volume vs birth volume
%        plot (2) single cell added volume vs time of birth
%        plot (3) single cell added volume vs time of division


%  Strategy: 
%
%       a) initialize experimental data
%       b) identify complete cell cycles within each condition 
%       c) compile birth volumes and added mass, using isDrop
%       d) calculate mean and plot condition data
%



%  Last edit: jen, 2018 May 2
%  Commit: added volume vs (A) birth vol, (B) time of birth, and (C) time
%  of division. single cell data from 2018-02-01


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
%bubbletime = storedMetaData{index}.bubbletime;


% 2. load measured data
experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
cd(experimentFolder)
filename = strcat('lb-fluc-',date,'-window5-width1p4-1p7-jiggle-0p5.mat');
load(filename,'D5','M','M_va','T');


% 3. compile experiment data matrix
xy_start = min(min(storedMetaData{index}.xys));
xy_end = max(max(storedMetaData{index}.xys));
exptData = buildDM(D5, M, M_va, T, xy_start, xy_end,e);
clear D5 M M_va T xy_start xy_end


% 4. initialize colors for plotting
palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
environment = {'fluc','low','ave','high'};


%%
for condition = 1:length(environment)
    
    % 6. isolate condition specific data
    conditionData = exptData(exptData(:,23) == condition,:);  % col 23 = cond vals
    
    
    % 7. trim data to full cell cycles ONLY
    ccFraction = conditionData(:,9);            % col 9 = ccFraction
    conditionData_fullOnly = conditionData(~isnan(ccFraction),:);
    
    
    % 8. isolate corrected time, added mass, volume, and birth event data (drop)
    curveFinder = conditionData_fullOnly(:,6);      % col 6   = curve Finder
    addedVol = conditionData_fullOnly(:,15);        % col 15  = added vol (calculated from Va)
    volume = conditionData_fullOnly(:,12);          % col 12  = volume (Va)
    isBirth = conditionData_fullOnly(:,5);          % col 5   = isDrop (1 = birth event, 0 = not)
    correctedTime = conditionData_fullOnly(:,25)/3600;   % col 25  = corrected time
    
    
    % 9. for each unique cell cycle, collect V_birth and time
    V_added = addedVol(isBirth == 1);
    V_birth = volume(isBirth == 1);
    
    unique_cycles = unique(curveFinder);
    for cc = 1:length(unique_cycles)
        
        currentTimes = correctedTime(curveFinder == unique_cycles(cc));
        timestamps_division(cc,1) = currentTimes(end);
        timestamps_birth(cc,1) = currentTimes(1);
        
    end

    % 10. plot
    color = rgb(palette(condition));
    
    % added volume vs. birth size
    figure(18)
    subplot(2,2,condition)
    plot(V_birth,V_added,'o','Color',color)
    hold on
%     plot(expected_Vb,expected_Vadded,'-','Color',expected_color);
%     hold on
    xlabel('birth volume (cubic um)')
    ylabel('added volume (cubic um)')
    title(strcat(date,': all (',num2str(sum(isBirth)),') cycles'))
    axis([-5 20 -5 20])
    legend(environment(condition))
    
    % added volume vs. birth time
    figure(19)
    subplot(2,2,condition)
    plot(timestamps_birth,V_added,'o','Color',color)
    hold on
    xlabel('time of birth (hr)')
    ylabel('added volume (cubic um)')
    title(strcat(date,': all (',num2str(sum(isBirth)),') cycles'))
    legend(environment(condition))
    
    % added volume vs. birth time
    figure(20)
    subplot(2,2,condition)
    plot(timestamps_division,V_added,'o','Color',color)
    hold on
    xlabel('time of division (hr)')
    ylabel('added volume (cubic um)')
    title(strcat(date,': all (',num2str(sum(isBirth)),') cycles'))
    legend(environment(condition))
    
    clear timestamps_birth timestamps_division currentTimes
       
end


