%% figure 7
%
%  Goals: plot instantaneous dV/dt over time
%
%  
% 

% last updated: jen, 2018 May 8
% commit: plot binned and averaged instantaneous dV/dt vs time for each 15
%         min and 60 min period experiment

% OK let's go!

%% 
clear
clc

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')

dataIndex = find(~cellfun(@isempty,storedMetaData));

% 1. for all experiments in dataset
binsPerHour = 30;
exptCounter = 0; % experiment counter

for e = 8:11
    
    % 2. collect experiment date
    index = dataIndex(e);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    
    % exclude outliers from analysis (2017-10-31 and monod experiments)
    if strcmp(date, '2017-10-31') == 1 || strcmp (timescale, 'monod') == 1
        disp(strcat(date,': excluded from analysis'))
        continue
    end
    disp(strcat(date, ': analyze!'))
    exptCounter = exptCounter + 1;
    
    
    % 3. initialize experiment meta data
    xys = storedMetaData{index}.xys;
    bubbletime = storedMetaData{index}.bubbletime;
    
    
    % 4. load measured data
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-window5-width1p4-1p7-jiggle-0p5.mat');
    load(filename,'D5','M','M_va','T');
    
    
    % 5. build data matrix from specified condition
    for condition = 1:4 
    
        xy_start = storedMetaData{index}.xys(condition,1);
        xy_end = storedMetaData{index}.xys(condition,end);
        conditionData = buildDM(D5, M, M_va, T, xy_start, xy_end,e);
        
        
        % 6. isolate condition data to those with full cell cycles
        curveIDs = conditionData(:,6);           % col 6 = curve ID
        conditionData_fullOnly = conditionData(curveIDs > 0,:);
        clear curveFinder
        
        
        % 7. isolate data to stabilized regions of growth
        maxTime = bubbletime(condition);
        timestamps = conditionData_fullOnly(:,2)/3600; % time in seconds converted to hours
        
        if maxTime > 0
            conditionData_bubbleTrimmed = conditionData_fullOnly(timestamps <= maxTime,:);
        else
            conditionData_bubbleTrimmed = conditionData_fullOnly;
        end
        clear timestamps maxTime
        
        
        % 8. isolate volume (Va) and timestamp data
        volumes = conditionData_bubbleTrimmed(:,12);        % col 12 = calculated va_vals (cubic um)
        timestamps = conditionData_bubbleTrimmed(:,2);      % col 2  = timestamp in seconds
        %isDrop = conditionData_bubbleTrimmed(:,5);          % col 5  = isDrop, 1 marks a birth event
        curveFinder = conditionData_bubbleTrimmed(:,6);     % col 6  = curve finder (ID of curve in condition)
        
        
        % 9. calculate dV/dt
        dVdt = [];
        curveIDs = unique(curveFinder);
        for cc = 1:length(curveIDs)
            
            if length(curveFinder == curveIDs(cc)) == 1
                continue
            end
            
            % (i) isolate volumes for current curve
            currentCurve_indeces = find(curveFinder == cc);
            currentCurve_volumes = volumes(currentCurve_indeces);
            currentCurve_times = timestamps(currentCurve_indeces); % sec
          
            dt = [NaN; diff(currentCurve_times)]; % timestep in seconds
            dV_raw = [NaN; diff(currentCurve_volumes)];
            
            currentCurve_dVdt = dV_raw./dt*3600;            % final units = cubic um/hr
            
            dVdt = [dVdt; currentCurve_dVdt];
            
        end
        clear dV_raw
        
        
        % 10. bin dV/dt into time bins by assigning timestamps to bins
        timeInHours = timestamps/3600;
        bins = ceil(timeInHours*binsPerHour);  
        
        % remove nans from dvdt
        dVdt_noNaNs = dVdt(~isnan(dVdt),:);
        bins_noNaNs = bins(~isnan(dVdt),:);
        
        binned_dvdt = accumarray(bins_noNaNs,dVdt_noNaNs,[],@(x) {x});
        bin_means = cellfun(@mean,binned_dvdt);
        bin_stds = cellfun(@std,binned_dvdt);
        bin_counts = cellfun(@length,binned_dvdt);
        bin_sems = bin_stds./sqrt(bin_counts);
        
        binVector = linspace(1,binsPerHour*10,binsPerHour*10);
  
        
        % 11. plot
        palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
        shapes = {'o','*','square'};
        
        color = rgb(palette(condition));
        xmark = shapes{1};
        
        figure(e)
        errorbar(binVector(1:length(bin_means))/binsPerHour,bin_means,bin_sems,'Color',color)
        hold on
        plot(binVector(1:length(bin_means))/binsPerHour,bin_means,'Color',color,'Marker',xmark)
        hold on
        grid on
        axis([0,10.1,-5,25])
        xlabel('Time (hr)')
        ylabel('mean instantaneous dV/dt (cubic um/hr)')
        title(date)
        legend('fluc','1/1000 LB','ave', '1/50 LB')
        
    end
end


%%

