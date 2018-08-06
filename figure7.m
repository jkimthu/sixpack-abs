%% figure 7
%
%  Goals: plot growth rate over time
%
%         version 1: dV/dt vs time
%         version 2: (dV/dt)/V vs time
% 

% last updated: jen, 2018 August 06

% commit: plot growth rate single upshift experiment, 2018-06-15, tracked
%         using the lower width bound 1.4 for fluc xys 1-10
%         

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

%%
for e = 17%:length(dataIndex)
    
    % 2. collect experiment date
    index = dataIndex(e);
    date = storedMetaData{index}.date;
    expType = storedMetaData{index}.experimentType;
    
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
    %filename = strcat('lb-fluc-',date,'-window5-width1p7-jiggle-0p5.mat');
    load(filename,'D5','M','M_va','T');
    
    
    % 5. build data matrix from specified condition
    for condition = 1:length(bubbletime)
    
        xy_start = storedMetaData{index}.xys(condition,1);
        xy_end = storedMetaData{index}.xys(condition,end);
        conditionData = buildDM(D5, M, M_va, T, xy_start, xy_end,e,expType);
        
        
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
        isDrop = conditionData_bubbleTrimmed(:,5);          % col 5  = isDrop, 1 marks a birth event
        curveFinder = conditionData_bubbleTrimmed(:,6);     % col 6  = curve finder (ID of curve in condition)
        
        
        % 9. calculate raw dV/dt and dV/dt normalized by initial volume
        curveIDs = unique(curveFinder);
        firstFullCurve = curveIDs(2);
        if length(firstFullCurve) > 1
            firstFullCurve_timestamps = timestamps(curveFinder == firstFullCurve);
        else
            firstFullCurve = curveIDs(3);
            firstFullCurve_timestamps = timestamps(curveFinder == firstFullCurve);
        end
        dt = mean(diff(firstFullCurve_timestamps)); % timestep in seconds
        
        dV_raw_noNan = diff(volumes);
        dV_raw = [NaN; dV_raw_noNan];
        dV_norm = [NaN; dV_raw_noNan./volumes(1:end-1)];
        
        dVdt_raw = dV_raw/dt * 3600;            % final units = cubic um/hr
        dVdt_overV = dV_norm/dt * 3600;         % final units = 1/hr           
        
        dVdt_raw(isDrop == 1) = NaN;
        dVdt_overV(isDrop == 1) = NaN;
        
        
        
        % 10. bin dV/dt into time bins by assigning timestamps to bins
        timeInHours = timestamps/3600;
        bins = ceil(timeInHours*binsPerHour);  
        
        % remove nans from raw dvdt
        dVdt_noNaNs = dVdt_raw(~isnan(dVdt_raw),:);
        bins_noNaNs = bins(~isnan(dVdt_raw),:);
        
        binned_dvdt_raw = accumarray(bins_noNaNs,dVdt_noNaNs,[],@(x) {x});
        bin_means_raw = cellfun(@mean,binned_dvdt_raw);
        bin_stds_raw = cellfun(@std,binned_dvdt_raw);
        bin_counts_raw = cellfun(@length,binned_dvdt_raw);
        bin_sems_raw = bin_stds_raw./sqrt(bin_counts_raw);
        
        % remove nans from dvdt normalized by initial volume
        dVdt_noNaNs_norm = dVdt_overV(~isnan(dVdt_overV),:);
        bins_noNaNs_norm = bins(~isnan(dVdt_overV),:);
        
        binned_dvdt_norm = accumarray(bins_noNaNs_norm,dVdt_noNaNs_norm,[],@(x) {x});
        bin_means_norm = cellfun(@mean,binned_dvdt_norm);
        bin_stds_norm = cellfun(@std,binned_dvdt_norm);
        bin_counts_norm = cellfun(@length,binned_dvdt_norm);
        bin_sems_norm = bin_stds_norm./sqrt(bin_counts_norm);
        
        binVector = linspace(1,binsPerHour*10,binsPerHour*10);
  
        
        % 11. plot raw dV/dt and normalized dV/dt over time
        palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
        shapes = {'o','*','square'};
        
        color = rgb(palette(condition));
        xmark = shapes{1};
        
        figure(e)
        errorbar(binVector(1:length(bin_means_raw))/binsPerHour,bin_means_raw,bin_sems_raw,'Color',color)
        hold on
        plot(binVector(1:length(bin_means_raw))/binsPerHour,bin_means_raw,'Color',color,'Marker',xmark)
        hold on
        grid on
        axis([0,10.1,-5,25])
        xlabel('Time (hr)')
        ylabel('mean instantaneous dV/dt (cubic um/hr)')
        title(date)
        legend('fluc','1/1000 LB','ave', '1/50 LB')
        
        figure(e+10)
        errorbar(binVector(1:length(bin_means_norm))/binsPerHour,bin_means_norm,bin_sems_norm,'Color',color)
        hold on
        plot(binVector(1:length(bin_means_norm))/binsPerHour,bin_means_norm,'Color',color,'Marker',xmark)
        hold on
        grid on
        axis([0,10.1,-1,5])
        xlabel('Time (hr)')
        ylabel('mean dV/dt/V (1/hr)')
        title(date)
        legend('fluc','1/1000 LB','ave', '1/50 LB')
        
    end
end


%%

