%% horoscope_growthTransitions


% Goal: are growth phenotypes predicted by when a cell is born?
%       plot growth rate before and after nutrient transitions,
%       as a function of birth phase 
%
%
%       Mean and scatter plots:
%          
%           1. growth rate in 30 sec prior to nutrient upshift
%           2. growth rate in 30 sec following nutrient upshift



%  Last edit: Jen Nguyen, 2018 Mar 2

%  Commit: plot growth rate before and after nutrient transitions, as a function of birth phase 
%


%  Strategy:
%
%     0. initialize complete meta data
%     0. for experiments of interest...
%           1. collect experiment date and exclude outliers (2017-10-31 and monod experiments)
%           2. initialize binning parameters
%           3. load measured data for stable condition
%           4. for stable average... 
%                  5. isolate relevant data
%                  6. trim data to timepoints after 3 hrs and before bubbles
%                  7. accumulate growth dadta for time window prior to (1) or following (2) upshift
%                           i. convert lag corrected timestamps into period fractions
%                          ii. generate highly resolved vector of period fraction
%                         iii. identify period fraction of upshift transition
%                                   - upshift is the transition between quarters 3 and 4
%                                   - downshift is the transition between quarters 1 and 2
%                          iv. identify period fraction values that bound time window
%                           v. trim data to include only points within window boundaries
%                  8. remove data associated with non-positive growth rates
%                  9. remove data not associated with a full curve
%                 10. bin growth data by nutrient phase of birth
%                           i. for each track and curve pair in isolated growth rate data,
%                              identify time of birth
%                          ii. isolate corresponding track data from dm trimmed only to stabilized region
%                         iii. identify and save corrected birth timestamp of paired curve
%                          iv. assign birth times to period fraction bin
%                           v. ...at last! bin growth data by period fraction
%                 11. plot mean and s.e.m. of growth rates immediately before transition over nutrient period        
%          12. repeat analysis for fluctuating environment, plotting fluc data over stable


% OK! Lez go!


%% ONE. growth rate immediately prior to nutrient upshift


% 0. initialize data

clc
clear

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));


% for target experiments...
for e = 6:14
    
    % 1. collect experiment date
    index = dataIndex(e);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    bubbleTime = storedMetaData{index}.bubbletime;
    
    % exclude outliers from analysis (2017-10-31 and monod experiments)
    if strcmp(date, '2017-10-31') == 1 || strcmp (timescale, 'monod') == 1
       disp(strcat(date,': excluded from analysis'))
       continue
    end
    disp(strcat(date, ': analyze!'))

    
    % 2. initialize binning parameters
    binsPerPeriod = 12;
    
    
    % 3. load measured data in stable average condition
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-window5-width1p4-1p7-jiggle-0p5.mat');
    load(filename,'D','D5','M','M_va','T');
    
    % 4. for stable average and then fluctuating environment...
    environment = [3;1];
    
    for i = 1:length(environment)
        
        condition = environment(i);
        xy_start = storedMetaData{index}.xys(condition,1);
        xy_end = storedMetaData{index}.xys(condition,end);
        conditionData = buildDM(D5, M, M_va, T, xy_start, xy_end, e);
        
        % 5. isolate relevant data
        isDrops = conditionData(:,5);           % col 5   =  isDrop; 0 during curve, 1 at birth event
        correctedTimes = conditionData(:,30);   % col 30  =  true times after lag correction
        bvpr = conditionData(:,29);             % col 29  =  biovol prod rate
        mu_va = conditionData(:,17);            % col 17  =  mu, calculated from volume approx'd as a cylinder with spherical capss
        curveID = conditionData(:,6);           % col 6   =  curve ID per track
        trackNum = conditionData(:,27);         % col 27  =  track number, not ID from particle tracking 
        
        
        % 6i. trim data to timepoints after 3 hrs 
        data = [isDrops correctedTimes bvpr mu_va curveID trackNum];
        data_trim1 = data(correctedTimes/3600 >= 3,:);
         
        % 6ii. trim data to timepoints before appearance of bubbles
        if bubbleTime(condition) > 0
            data_trim2 = data_trim1(data_trim1(:,2)/3600 <= bubbleTime(condition),:);
        else
            data_trim2 = data_trim1;        
        end
        
        
        % 7. ACCUMULATE GROWTH DATA for time window prior to upshift
        timeWindow = 30; % sec
        window_asPeriodFraction = timeWindow/timescale;
        
        % 7i. convert lag corrected timestamps into period fractions
        correctedTimes_trimmed = data_trim2(:,2); % still in sec
        timeInPeriods = correctedTimes_trimmed/timescale; % unit = sec/sec
        timeInPeriodFraction = timeInPeriods - floor(timeInPeriods);
        
        % 7ii. generate highly resolved vector of period fraction
        resolved_periodFractions = linspace(0,1,10000)';
        timeInQuarters = ceil(resolved_periodFractions * 4);
        
        % 7iii. identify period fraction of upshift transition
        %       upshift is the transition between quarter 3 and 4.
        q4_indeces = find(timeInQuarters == 4);
        
        % 7iv. identify period fraction values that bound time window
        boundingWindow_end = resolved_periodFractions(min(q4_indeces)-1);
        boundingWindow_start = boundingWindow_end - window_asPeriodFraction;
        
        % 7v. trim data to include only points within window boundaries
        timeInPeriodFraction_beforeUpshift = timeInPeriodFraction(timeInPeriodFraction <= boundingWindow_end,:);
        data_beforeUpshift = data_trim2(timeInPeriodFraction <= boundingWindow_end,:);
        
        timeInPeriodFraction_inWindow = timeInPeriodFraction_beforeUpshift(timeInPeriodFraction_beforeUpshift >= boundingWindow_start,:);
        data_inWindow = data_beforeUpshift(timeInPeriodFraction_beforeUpshift >= boundingWindow_start,:);
        
        
        % 8. remove data associated with non-positive growth rates
        mus_inWindow = data_inWindow(:,4);
        data_inWindow_positive = data_inWindow(mus_inWindow > 0,:);
        clear mus_inWindow
        
        
        % 9. remove data not associated with a full curve
        curves_positive = data_inWindow_positive(:,5);
        data_inWindow_fullCurves = data_inWindow_positive(curves_positive > 0,:);
        clear curves_positive
        
        
        % 10.  BIN GROWTH DATA by nutrient phase of birth
        % 10. (i) for each track and curve pair in isolated growth rate data, identify time of birth
        tracks_inWindow = data_inWindow_fullCurves(:,6);
        curve_birthTimes = nan(length(tracks_inWindow),1);
        for m = 1:length(tracks_inWindow)
            
            % 10. (ii) isolate corresponding track data from dm trimmed only to stabilized region
            currentTrack = data_inWindow_fullCurves(m,6);
            trackData = data_trim2(data_trim2(:,6) == currentTrack,:);
            
            % 10. (iii) identify and save corrected birth timestamp of paired curve
            pairedCurve = data_inWindow_fullCurves(m,5);
            curveData = trackData(trackData(:,5) == pairedCurve,:);
            curve_birthTimes(m) = curveData(1,2); % sec            
            
        end
        
        % 10. (iv) assign birth times to period fraction bin
        curve_birthTimeinPeriods = curve_birthTimes/timescale; % sec/sec
        curve_birthTimeinPeriodFraction = curve_birthTimeinPeriods - floor(curve_birthTimeinPeriods);        
        assignedBin = ceil(curve_birthTimeinPeriodFraction * binsPerPeriod);
        
        
        % 11. (v) ...at last! bin growth data by period fraction
        bvpr_final = data_inWindow_fullCurves(:,3);
        mus_final = data_inWindow_fullCurves(:,4);
        
        bvpr_binnedByPeriodFraction = accumarray(assignedBin, bvpr_final, [], @(x) {x});
        mus_binnedByPeriodFraction = accumarray(assignedBin, mus_final, [], @(x) {x});
        
        % if no data for all cells, add empty cells to fill out period
        if length(bvpr_binnedByPeriodFraction) < binsPerPeriod
            
            bvpr_binnedByPeriodFraction{binsPerPeriod,1} = [];
            mus_binnedByPeriodFraction{binsPerPeriod,1} = [];
            
        end
        
        bvpr_binned_mean = cellfun(@mean,bvpr_binnedByPeriodFraction);
        bvpr_binned_std = cellfun(@std,bvpr_binnedByPeriodFraction);
        bvpr_binned_count = cellfun(@length,bvpr_binnedByPeriodFraction);
        bvpr_binned_sem = bvpr_binned_std./bvpr_binned_count;
        
        mus_binned_mean = cellfun(@mean,mus_binnedByPeriodFraction);
        mus_binned_std = cellfun(@std,mus_binnedByPeriodFraction);
        mus_binned_count = cellfun(@length,mus_binnedByPeriodFraction);
        mus_binned_sem = mus_binned_std./mus_binned_count;
    
        
        
        % 11. plot mean and s.e.m. of growth rates immediately before transition over nutrient period 
        % 11. (i) convert bin # to absolute time (sec)
        timePerBin = timescale/binsPerPeriod;  % in sec
        binPeriod = linspace(1, binsPerPeriod, binsPerPeriod);
        timePeriod = timePerBin*binPeriod';

        % 11. (ii) repeat quarter period on both sides and plot over period fraction
        quarterOne = linspace(1,binsPerPeriod/4,binsPerPeriod/4);
        quarterFour = linspace(binsPerPeriod*3/4+1,binsPerPeriod,binsPerPeriod/4);
        quarterZero = linspace((binsPerPeriod/4-1),0,binsPerPeriod/4)*-1;
        
        % mean
        bvprSignal_quarterOne = bvpr_binned_mean(quarterOne,1);
        bvprSignal_quarterFour = bvpr_binned_mean(quarterFour,1);
        musSignal_quarterOne = mus_binned_mean(quarterOne,1);
        musSignal_quarterFour = mus_binned_mean(quarterFour,1);
        
        % sem
        bvprSEM_quarterOne = bvpr_binned_sem(quarterOne,1);
        bvprSEM_quarterFour = bvpr_binned_sem(quarterFour,1);
        musSEM_quarterOne = mus_binned_sem(quarterOne,1);
        musSEM_quarterFour = mus_binned_sem(quarterFour,1);

        time_quarterOne = timePeriod(quarterOne,1) + timescale;
        time_quarterZero = quarterZero*timePerBin;
        
        mean_bvpr_stitched = [bvprSignal_quarterFour; bvpr_binned_mean; bvprSignal_quarterOne];
        mean_mus_stitched = [musSignal_quarterFour; mus_binned_mean; musSignal_quarterOne];
        
        sem_bvpr_stitched = [bvprSEM_quarterFour; bvpr_binned_sem; bvprSEM_quarterOne];
        sem_mus_stitched = [musSEM_quarterFour; mus_binned_sem; musSEM_quarterOne];
        
        time_stitched = [time_quarterZero'; timePeriod; time_quarterOne];
  
        % bvpr
        figure(1)
        if condition == 3
            subplot(1,2,1)
            stem(time_stitched, mean_bvpr_stitched,'Color',[0.25 0.25 0.9])
            hold on
            axis([min(time_stitched) max(time_stitched) 0 20])
            title(strcat('stable: ',date))
        elseif condition == 1
            subplot(1,2,2)
            stem(time_stitched, mean_bvpr_stitched,'Color',[0 0.7 0.7])
            axis([min(time_stitched) max(time_stitched) 0 20])
            title(strcat('fluc: ',num2str(timescale)))
        end
        xlabel('time (sec)')
        ylabel('mean bvpr (min)')
        grid on
        
         % mus
        figure(2)
        if condition == 3
            subplot(1,2,1)
            stem(time_stitched, mean_mus_stitched,'Color',[0.25 0.25 0.9])
            hold on
            axis([min(time_stitched) max(time_stitched) 0 3])
            title(strcat('stable: ',date))
        elseif condition == 1
            subplot(1,2,2)
            stem(time_stitched, mean_mus_stitched,'Color',[0 0.7 0.7])
            axis([min(time_stitched) max(time_stitched) 0 3])
            title(strcat('fluc: ',num2str(timescale)))
        end
        xlabel('time (sec)')
        ylabel('mean mus (min)')
        grid on
        
        
        
        % 11. plot mean of growth rates with standard error
        % bvpr
        figure(3)
        if condition == 3
            subplot(1,2,1)
            errorbar(time_stitched, mean_bvpr_stitched,sem_bvpr_stitched,'Color',[0.25 0.25 0.9],'Marker','o')
            hold on
            axis([min(time_stitched) max(time_stitched) 0 20])
            title(strcat('stable: ',date))
        elseif condition == 1
            subplot(1,2,2)
            errorbar(time_stitched, mean_bvpr_stitched, sem_bvpr_stitched, 'Color',[0 0.7 0.7],'Marker','o')
            axis([min(time_stitched) max(time_stitched) 0 20])
            title(strcat('fluc: ',num2str(timescale)))
        end
        xlabel('time (sec)')
        ylabel('mean bvpr (min)')
        grid on
        
        % mus
        figure(4)
        if condition == 3
            subplot(1,2,1)
            errorbar(time_stitched,mean_mus_stitched, sem_mus_stitched,'Color',[0.25 0.25 0.9],'Marker','o')
            hold on
            axis([min(time_stitched) max(time_stitched) 0 3])
            title(strcat('stable: ',date))
        elseif condition == 1
            subplot(1,2,2)
            errorbar(time_stitched, mean_mus_stitched, sem_mus_stitched,'Color',[0 0.7 0.7],'Marker','o')
            axis([min(time_stitched) max(time_stitched) 0 3])
            title(strcat('fluc: ',num2str(timescale)))
        end
        xlabel('time (sec)')
        ylabel('mean mus with std error (1/hr)')
        grid on
        
       
        % 12. repeat analysis for fluctuating environment, plotting fluc data over stable
    end
    clear i environment 
    
    cd('/Users/jen/Documents/StockerLab/Data_analysis/horoscope')
    
    Fig1 = figure(1);
    saveas(Fig1,strcat('horoscope-growthTransitions-mean-bvpr-',num2str(timescale),'-',date),'epsc')
    close(Fig1)
    
    Fig2 = figure(2);
    saveas(Fig2,strcat('horoscope-growthTransitions-mean-mus-',num2str(timescale),'-',date),'epsc')
    close(Fig2)
    
    Fig3 = figure(3);
    saveas(Fig3,strcat('horoscope-growthTransitions-mean&sem-bvpr-',num2str(timescale),'-',date),'epsc')
    close(Fig3)
    
    Fig4 = figure(4);
    saveas(Fig4,strcat('horoscope-growthTransitions-mean&sem-mus-',num2str(timescale),'-',date),'epsc')
    close(Fig4)
    
end

%% TWO. growth rate immediately after nutrient upshift

% 0. initialize data

clc
clear

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));
                            

% for target experiments...
for e = 2:4
    
    % 1. collect experiment date
    index = dataIndex(e);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    bubbleTime = storedMetaData{index}.bubbletime;
    
    % exclude outliers from analysis (2017-10-31 and monod experiments)
    if strcmp(date, '2017-10-31') == 1 || strcmp (timescale, 'monod') == 1
       disp(strcat(date,': excluded from analysis'))
       continue
    end
    disp(strcat(date, ': analyze!'))

    
    % 2. initialize binning parameters
    binsPerPeriod = 12;
    
    
    % 3. load measured data in stable average condition
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-window5-width1p4-1p7-jiggle-0p5.mat');
    load(filename,'D','D5','M','M_va','T');
    
    % 4. for stable average and then fluctuating environment...
    environment = [3;1];
    
    for i = 1:length(environment)
        
        condition = environment(i);
        xy_start = storedMetaData{index}.xys(condition,1);
        xy_end = storedMetaData{index}.xys(condition,end);
        conditionData = buildDM(D5, M, M_va, T, xy_start, xy_end, e);
        
        % 5. isolate relevant data
        isDrops = conditionData(:,5);           % col 5   =  isDrop; 0 during curve, 1 at birth event
        correctedTimes = conditionData(:,30);   % col 30  =  true times after lag correction
        bvpr = conditionData(:,29);             % col 29  =  biovol prod rate
        mu_va = conditionData(:,17);            % col 17  =  mu, calculated from volume approx'd as a cylinder with spherical capss
        curveID = conditionData(:,6);           % col 6   =  curve ID per track
        trackNum = conditionData(:,27);         % col 27  =  track number, not ID from particle tracking 
        
%         % 6. calculate nScores for condition
%         [~, nScore] = nutrientScore(timescale,conditionData);
%        
        
        % 7i. trim data to timepoints after 3 hrs 
        data = [isDrops correctedTimes bvpr mu_va curveID trackNum];
        data_trim1 = data(correctedTimes/3600 >= 3,:);
         
        % 7ii. trim data to timepoints before appearance of bubbles
        if bubbleTime(condition) > 0
            data_trim2 = data_trim1(data_trim1(:,2)/3600 <= bubbleTime(condition),:);
        else
            data_trim2 = data_trim1;        
        end
        
        
        % 8. ACCUMULATE GROWTH DATA for time window prior to upshift
        timeWindow = 30; % sec
        window_asPeriodFraction = timeWindow/timescale;
        
        % 8i. convert lag corrected timestamps into period fractions
        correctedTimes_trimmed = data_trim2(:,2); % still in sec
        timeInPeriods = correctedTimes_trimmed/timescale; % unit = sec/sec
        timeInPeriodFraction = timeInPeriods - floor(timeInPeriods);
        
        % 8ii. generate highly resolved vector of period fraction
        resolved_periodFractions = linspace(0,1,10000)';
        timeInQuarters = ceil(resolved_periodFractions * 4);
        
        % 8iii. identify period fraction of upshift transition
        %       upshift is the transition between quarter 1 and 2.
        q2_indeces = find(timeInQuarters == 2);
        
        % 8iv. identify period fraction values that bound time window
        boundingWindow_start = resolved_periodFractions(min(q2_indeces));
        boundingWindow_end = boundingWindow_start + window_asPeriodFraction;
        
        % 8v. trim data to include only points within window boundaries
        timeInPeriodFraction_beforeUpshift = timeInPeriodFraction(timeInPeriodFraction <= boundingWindow_end,:);
        data_afterUpshift = data_trim2(timeInPeriodFraction <= boundingWindow_end,:);
        
        timeInPeriodFraction_inWindow = timeInPeriodFraction_beforeUpshift(timeInPeriodFraction_beforeUpshift >= boundingWindow_start,:);
        data_inWindow = data_afterUpshift(timeInPeriodFraction_beforeUpshift >= boundingWindow_start,:);
        
        
        % 9. remove data associated with non-positive growth rates
        mus_inWindow = data_inWindow(:,4);
        data_inWindow_positive = data_inWindow(mus_inWindow > 0,:);
        clear mus_inWindow
        
        
        % 10. remove data not associated with a full curve
        curves_positive = data_inWindow_positive(:,5);
        data_inWindow_fullCurves = data_inWindow_positive(curves_positive > 0,:);
        clear curves_positive
        
        
        % 11.  BIN GROWTH DATA by nutrient phase of birth
        % 11. (i) for each track and curve pair in isolated growth rate data, identify time of birth
        tracks_inWindow = data_inWindow_fullCurves(:,6);
        curve_birthTimes = nan(length(tracks_inWindow),1);
        for m = 1:length(tracks_inWindow)
            
            % 11. (ii) isolate corresponding track data from dm trimmed only to stabilized region
            currentTrack = data_inWindow_fullCurves(m,6);
            trackData = data_trim2(data_trim2(:,6) == currentTrack,:);
            
            % 11. (iii) identify and save corrected birth timestamp of paired curve
            pairedCurve = data_inWindow_fullCurves(m,5);
            curveData = trackData(trackData(:,5) == pairedCurve,:);
            curve_birthTimes(m) = curveData(1,2); % sec            
            
        end
        
        % 11. (iv) assign birth times to period fraction bin
        curve_birthTimeinPeriods = curve_birthTimes/timescale; % sec/sec
        curve_birthTimeinPeriodFraction = curve_birthTimeinPeriods - floor(curve_birthTimeinPeriods);        
        assignedBin = ceil(curve_birthTimeinPeriodFraction * binsPerPeriod);
        
        
        % 11. (v) ...at last! bin growth data by period fraction
        bvpr_final = data_inWindow_fullCurves(:,3);
        mus_final = data_inWindow_fullCurves(:,4);
        
        bvpr_binnedByPeriodFraction = accumarray(assignedBin, bvpr_final, [], @(x) {x});
        mus_binnedByPeriodFraction = accumarray(assignedBin, mus_final, [], @(x) {x});
        
        % if no data for all cells, add empty cells to fill out period
        if length(bvpr_binnedByPeriodFraction) < binsPerPeriod
            
            bvpr_binnedByPeriodFraction{binsPerPeriod,1} = [];
            mus_binnedByPeriodFraction{binsPerPeriod,1} = [];
            
        end
        
        bvpr_binned_mean = cellfun(@mean,bvpr_binnedByPeriodFraction);
        bvpr_binned_std = cellfun(@std,bvpr_binnedByPeriodFraction);
        bvpr_binned_count = cellfun(@length,bvpr_binnedByPeriodFraction);
        bvpr_binned_sem = bvpr_binned_std./bvpr_binned_count;
        
        mus_binned_mean = cellfun(@mean,mus_binnedByPeriodFraction);
        mus_binned_std = cellfun(@std,mus_binnedByPeriodFraction);
        mus_binned_count = cellfun(@length,mus_binnedByPeriodFraction);
        mus_binned_sem = mus_binned_std./mus_binned_count;
    
        
        
        % 12. plot mean and s.e.m. of growth rates immediately before transition over nutrient period 
        % 12. (i) convert bin # to absolute time (sec)
        timePerBin = timescale/binsPerPeriod;  % in sec
        binPeriod = linspace(1, binsPerPeriod, binsPerPeriod);
        timePeriod = timePerBin*binPeriod';

        % 12. (ii) repeat quarter period on both sides and plot over period fraction
        quarterOne = linspace(1,binsPerPeriod/4,binsPerPeriod/4);
        quarterFour = linspace(binsPerPeriod*3/4+1,binsPerPeriod,binsPerPeriod/4);
        quarterZero = linspace((binsPerPeriod/4-1),0,binsPerPeriod/4)*-1;
        
        % mean
        bvprSignal_quarterOne = bvpr_binned_mean(quarterOne,1);
        bvprSignal_quarterFour = bvpr_binned_mean(quarterFour,1);
        musSignal_quarterOne = mus_binned_mean(quarterOne,1);
        musSignal_quarterFour = mus_binned_mean(quarterFour,1);
        
        % sem
        bvprSEM_quarterOne = bvpr_binned_sem(quarterOne,1);
        bvprSEM_quarterFour = bvpr_binned_sem(quarterFour,1);
        musSEM_quarterOne = mus_binned_sem(quarterOne,1);
        musSEM_quarterFour = mus_binned_sem(quarterFour,1);

        time_quarterOne = timePeriod(quarterOne,1) + timescale;
        time_quarterZero = quarterZero*timePerBin;
        
        mean_bvpr_stitched = [bvprSignal_quarterFour; bvpr_binned_mean; bvprSignal_quarterOne];
        mean_mus_stitched = [musSignal_quarterFour; mus_binned_mean; musSignal_quarterOne];
        
        sem_bvpr_stitched = [bvprSEM_quarterFour; bvpr_binned_sem; bvprSEM_quarterOne];
        sem_mus_stitched = [musSEM_quarterFour; mus_binned_sem; musSEM_quarterOne];
        
        time_stitched = [time_quarterZero'; timePeriod; time_quarterOne];
  
        % bvpr
        figure(1)
        if condition == 3
            subplot(1,2,1)
            stem(time_stitched, mean_bvpr_stitched,'Color',[0.25 0.25 0.9])
            hold on
            axis([min(time_stitched) max(time_stitched) 0 20])
            title(strcat('stable: ',date))
        elseif condition == 1
            subplot(1,2,2)
            stem(time_stitched, mean_bvpr_stitched,'Color',[0 0.7 0.7])
            axis([min(time_stitched) max(time_stitched) 0 20])
            title(strcat('fluc: ',num2str(timescale)))
        end
        xlabel('time (sec)')
        ylabel('mean bvpr (min)')
        grid on
        
         % mus
        figure(2)
        if condition == 3
            subplot(1,2,1)
            stem(time_stitched, mean_mus_stitched,'Color',[0.25 0.25 0.9])
            hold on
            axis([min(time_stitched) max(time_stitched) 0 3])
            title(strcat('stable: ',date))
        elseif condition == 1
            subplot(1,2,2)
            stem(time_stitched, mean_mus_stitched,'Color',[0 0.7 0.7])
            axis([min(time_stitched) max(time_stitched) 0 3])
            title(strcat('fluc: ',num2str(timescale)))
        end
        xlabel('time (sec)')
        ylabel('mean mus (min)')
        grid on
        
        
        
        % 12. plot mean of growth rates with standard error
        % bvpr
        figure(3)
        if condition == 3
            subplot(1,2,1)
            errorbar(time_stitched, mean_bvpr_stitched,sem_bvpr_stitched,'Color',[0.25 0.25 0.9],'Marker','o')
            hold on
            axis([min(time_stitched) max(time_stitched) 0 20])
            title(strcat('stable: ',date))
        elseif condition == 1
            subplot(1,2,2)
            errorbar(time_stitched, mean_bvpr_stitched, sem_bvpr_stitched, 'Color',[0 0.7 0.7],'Marker','o')
            axis([min(time_stitched) max(time_stitched) 0 20])
            title(strcat('fluc: ',num2str(timescale)))
        end
        xlabel('time (sec)')
        ylabel('mean bvpr (min)')
        grid on
        
        % mus
        figure(4)
        if condition == 3
            subplot(1,2,1)
            errorbar(time_stitched,mean_mus_stitched, sem_mus_stitched,'Color',[0.25 0.25 0.9],'Marker','o')
            hold on
            axis([min(time_stitched) max(time_stitched) 0 3])
            title(strcat('stable: ',date))
        elseif condition == 1
            subplot(1,2,2)
            errorbar(time_stitched, mean_mus_stitched, sem_mus_stitched,'Color',[0 0.7 0.7],'Marker','o')
            axis([min(time_stitched) max(time_stitched) 0 3])
            title(strcat('fluc: ',num2str(timescale)))
        end
        xlabel('time (sec)')
        ylabel('mean mus with std error (1/hr)')
        grid on
        
       
        % 13. repeat analysis for fluctuating environment, plotting fluc data over stable
    end
    clear i environment 
    
    cd('/Users/jen/Documents/StockerLab/Data_analysis/horoscope')
    
    Fig1 = figure(1);
    saveas(Fig1,strcat('horoscope-growthTransitions_after-mean-bvpr-',num2str(timescale),'-',date),'epsc')
    close(Fig1)
    
    Fig2 = figure(2);
    saveas(Fig2,strcat('horoscope-growthTransitions_after-mean-mus-',num2str(timescale),'-',date),'epsc')
    close(Fig2)
    
    Fig3 = figure(3);
    saveas(Fig3,strcat('horoscope-growthTransitions_after-mean&sem-bvpr-',num2str(timescale),'-',date),'epsc')
    close(Fig3)
    
    Fig4 = figure(4);
    saveas(Fig4,strcat('horoscope-growthTransitions_after-mean&sem-mus-',num2str(timescale),'-',date),'epsc')
    close(Fig4)
    
end

