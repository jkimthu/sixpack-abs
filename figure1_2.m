%% figures 1 and 2


%  Goal: determine whether growth in low, ave, high and fluc environments
%        are characterizable as adder, sizer, or timer

%        plot (1) population-averaged division size vs birth size
%        plot (2) population-averaged inter-division time vs birth time


%  Strategy: 
%
%       a) initialize experimental data
%       b) identify complete cell cycles within each condition 
%       c) compile birth volumes, division volumes and inter-division times
%       d) calculate mean and plot condition data
%



%  Last edit: jen, 2018 May 1
%  Commit: adder or sizer? population-average division volume vs birth volume,
%          and inter-div time vs. birth volume


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
    datesForLegend{exptCounter} = date;
    
    
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
    V_div_means = mean(V_division);
    interdiv_means = mean(unique_interdivs);
    
    V_birth_stds = std(V_birth);
    V_div_stds = std(V_division);
    interdiv_stds = std(unique_interdivs);
    
    % (i) inter-division time vs. birth size
    figure(3)
    errorbar(V_birth_means,interdiv_means,interdiv_stds,'Color',color)
    hold on
    plot(V_birth_means,interdiv_means,'Marker',xmark,'Color',color)
    hold on
    xlabel('birth size (cubic um)')
    ylabel('inter-division time (min)')
    title('population averages from all experiments')
    axis([0 10 0 100])
    
    
    % (ii) division size vs. birth size
    figure(4)
    errorbar(V_birth_means,V_div_means,V_div_stds,'Color',color)
    hold on
    plot(V_birth_means,V_div_means,'Marker',xmark,'Color',color)
    hold on
    xlabel('birth size (cubic um)')
    ylabel('division size (cubic um)')
    title('population averages from all experiments')
    axis([0 8 0 17])

    end  
    
    
end
