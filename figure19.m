%% figure 19


%  Goal: determine whether growth in low, ave, high and fluc environments
%        are characterizable as adder, sizer, or timer

%        plot (1) population-level V_div/V_birth ratio vs time of division


%  Strategy: 
%
%       a) initialize experimental data
%       b) identify complete cell cycles within each condition 
%       c) compile birth volumes and added mass, using isDrop
%       d) calculate mean and plot condition data
%



%  Last edit: jen, 2018 May 3
%  Commit: population-averaged ratio between V_div and V_birth vs time
%  of division.


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
for e = 1:experimentCount
    
    % 1. collect experiment meta data
    index = dataIndex(e);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    
    % exclude outliers from analysis (2017-10-31 and monod experiments)
    if strcmp(date, '2017-10-31') == 1 || strcmp (timescale, 'monod') == 1
        disp(strcat(date,': excluded from analysis'))
        continue
    end
    disp(strcat(date, ': analyze!'))
    
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
    
    
    % 4. initialize parameters for plotting
    palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
    shapes = {'o','x','square','*'};
    binsPerHour = 2;
    
    
    for condition = 1:length(palette)
        
        % 5. isolate condition specific data
        conditionData = exptData(exptData(:,23) == condition,:);  % col 23 = cond vals
        
        
        % 6. trim data to full cell cycles ONLY
        ccFraction = conditionData(:,9);            % col 9 = ccFraction
        conditionData_fullOnly = conditionData(~isnan(ccFraction),:);
        
        
        % 7. isolate corrected time, added mass, volume, and birth event data (drop)
        curveFinder = conditionData_fullOnly(:,6);           % col 6   = curve Finder
        volume = conditionData_fullOnly(:,12);               % col 12  = volume (Va)
        timestamp = conditionData_fullOnly(:,2)/3600;        % col 2   = corrected time
        
        
        % 8. for each unique cell cycle, collect V_birth, V_division and time
        unique_cycles = unique(curveFinder);
        clear V_birth V_div timestamps_division currentTimes currentVolumes
        for cc = 1:length(unique_cycles)
            
            currentTimes = timestamp(curveFinder == unique_cycles(cc));
            timestamps_division(cc,1) = currentTimes(end);
            
            currentVolumes = volume(curveFinder == unique_cycles(cc));
            V_birth(cc,1) = currentVolumes(1);
            V_div(cc,1) = currentVolumes(end);
            
        end
        
        
        % 9. calculate V_division:V_birth ratio
        V_ratio = V_div./V_birth;
        
        
        % 10. bin volume ratios into time bins
        bins_div = ceil(timestamps_division*binsPerHour);
        
        binnedRatios_mean = accumarray(bins_div,V_ratio,[],@mean);
        binnedRatios_std = accumarray(bins_div,V_ratio,[],@std);
        binnedRatios_counts = accumarray(bins_div,V_ratio,[],@length);
        binnedRatios_sem = binnedRatios_std./sqrt(binnedRatios_counts);
        
        
        % 11. if not all experiment bins in time vector are populated,
        %     extend vector with zeros for plotting
        timeVector = 1:(10*binsPerHour);
        
            
        % 12. plot
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
        
        % division:birth volume ratio vs. time of division
        figure(19)
        errorbar(timeVector(1,1:length(binnedRatios_mean)),binnedRatios_mean,binnedRatios_sem,'Color',color)
        hold on
        plot(timeVector(1,1:length(binnedRatios_mean),1),binnedRatios_mean,'Marker',xmark,'Color',color)
        hold on
        xlabel('time of division (hr)')
        ylabel('division volume-to-birth volume ratio')
        title('population-averaged V_div:V_birth ratio, binned by time of division')
        axis([0 timeVector(end)+1 0 5])
        
        
    end
    
end

% plot expected line of V_ratio = 2
figure(19)
hold on
plot(timeVector,ones(length(timeVector))*2,'Color',rgb('SlateGray'))

        