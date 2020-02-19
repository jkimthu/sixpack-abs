%%  added mass (per cell) over time


%  Goal: Searching responses in biomass accumulation
%        This script plots:
%
%        1a. mean added mass (per cell) since birth over time
%        1b. mean instantaneous added mass (per cell) over time
%        2. mean instantaneous added mass over a nutrient period


%  Last edit: Jen Nguyen, March 28th 2016


% The intended input for these scripts is the following data matrix,
% saved with the naming convention of:

% dmMMDD-cond.mat

%      where,
%              dm  =  dataMatrix                  (see matrixBuilder.m)
%              MM  =  month of experimental date
%              DD  =  day of experimental date
%       condition  =  experimental condition      (fluc or const)
%

% OK! Lez go!



% 0a. Load data matrices

dmDirectory = dir('dm*.mat'); % note: this assumes the only two data matrices are 'const' and 'fluc'
names = {dmDirectory.name}; % loaded alphabetically

for dm = 1:length(names)
    load(names{dm});                
    dataMatrices{dm} = dataMatrix;                                         % for entire condition
end                                                                        
clear dataMatrix dmDirectory dm;
clear names;


%%   O N E.
%    Added mass since birth over time
%    Instantaneous added mass over time



%    Strategy:
%
%     0. initialize experiment and analysis parameters
%     1. isolate data of interest
%     2. accumulate data points (of cell cycle stage) by time bin
%     3. calculate mean cell cycle stage per time bin
%     4. plot!
%     


% 0b. Initialize parameters

expHours = 10; %  duration of experiment in hours                          % 0.  initialize parameters
binFactor = 20; % time bins of 0.05 hr  
hrPerBin = 1/binFactor; 

% 0c. initialize time window parameters for analysis
firstHour = 3.7;                 % time at which to initate analysis
finalHour = 8.5;                % time at which to terminate analysis

for condition = 1:2
 
    % 1a.  isolate time and added mass data
    interestingData = dataMatrices{condition};  % condition: 1 = constant, 2 = fluctuating
    addedMass = interestingData(:,10);
    timeStamps = interestingData(:,2);
    
    %  trim off timepoints earlier than first
    addedMass = addedMass(timeStamps >= firstHour);
    lowTrimmed_timeStamps = timeStamps(timeStamps >= firstHour);
    
    %  trim off timepoints later than last
    addedMass = addedMass(lowTrimmed_timeStamps <= finalHour);
    finalTrimmed_timeStamps = lowTrimmed_timeStamps(lowTrimmed_timeStamps <= finalHour);
    
    % 1b.  calculate instantaneous added mass
    instaMass = diff(addedMass);
    instaMass = [0; instaMass];
   
    % 1c.  eliminate zeros (non-full track data and births) and negatives (divisions and noise) from data
    addedMass(addedMass <= 0) = NaN;
    instaMass(instaMass <= 0) = NaN;
    
    % 2.  accumulate data by associated time bin
    timeBins = ceil(finalTrimmed_timeStamps*binFactor);                                 
    binnedByTime_added = accumarray(timeBins,addedMass,[],@(x) {x});               
    binnedByTime_insta = accumarray(timeBins,instaMass,[],@(x) {x});
    
    % 3a.  calculate mean cell cycle stage per bin
    meanAdded = cell2mat( cellfun(@nanmean,binnedByTime_added,'UniformOutput',false) );      
    meanInsta = cell2mat( cellfun(@nanmean,binnedByTime_insta,'UniformOutput',false) );
    
    % 3b.  calculate std and error
    devAdded = cell2mat(cellfun(@nanstd,binnedByTime_added,'UniformOutput',false));
    nAdded = cell2mat(cellfun(@length,binnedByTime_added,'UniformOutput',false));
    errorAdded = devAdded./sqrt(nAdded);
    
    devInsta = cell2mat( cellfun(@nanstd,binnedByTime_insta,'UniformOutput',false) );
    nInsta = cell2mat( cellfun(@length,binnedByTime_insta,'UniformOutput',false) );
    errorInsta = devInsta./sqrt(nInsta);
    
    % 4a. create a time vector to convert bin # to absolute time
    timeVector = linspace(1, expHours/hrPerBin, expHours/hrPerBin);
    timeVector = hrPerBin*timeVector';                                       
    
    % 4b. before plotting, a little matrix manipulation to get around nans                                                            
    eventMask = find(~isnan(meanAdded));                                   % find indices of cell cycle data
    timeVector = timeVector(eventMask);                                    % trim all nans from time vector
    
    meanAdded = meanAdded(eventMask);                                      % trim all nans from mass vectors
    devAdded = devAdded(eventMask);
    errorAdded = errorAdded(eventMask);
    
    meanInsta = meanInsta(eventMask);                                      % same mask applies!
    devInsta = devInsta(eventMask);
    errorInsta = errorInsta(eventMask);
   
    % 4c. mean added (since birth) with standard deviation
    figure(1)
    if condition == 1
         plot(timeVector,meanAdded,'k')
         axis([firstHour,finalHour,0,3])
         hold on
         grid on
         %errorbar(timeVector,meanAdded, devAdded,'k')
    else
        plot(timeVector,meanAdded,'b')
        hold on
        %errorbar(timeVector,meanAdded,devAdded,'b')
    end

    % 4d. mean added (since birth) with standard error
    figure(2)
    if condition == 1
        plot(timeVector,meanAdded,'k')
        axis([firstHour,finalHour,0,3])
        hold on
        grid on
        errorbar(timeVector,meanAdded,errorAdded,'k')
    else
        plot(timeVector,meanAdded,'b')
        hold on
        errorbar(timeVector,meanAdded,errorAdded,'b')
    end
    
    
    % 4e. insta added with standard deviation
    figure(3)
    if condition == 1
         plot(timeVector,meanInsta,'k')
         axis([4,finalHour,0,.05])
         hold on
         grid on
         %errorbar(timeVector,meanInsta, devInsta,'k')
    else
        plot(timeVector,meanInsta,'b')
        hold on
        %errorbar(timeVector,meanInsta,devInsta,'b')
    end

    % 4d. insta added with standard error
    figure(4)
    if condition == 1
        plot(timeVector,meanInsta,'k')
        axis([4,finalHour,0,.05])
        hold on
        grid on
        errorbar(timeVector,meanInsta,errorInsta,'k')
    else
        plot(timeVector,meanInsta,'b')
        hold on
        errorbar(timeVector,meanInsta,errorInsta,'b')
    end
end


%%   T W O.
%    Instantaneous added mass over a nutrient period


%    Strategy:
%
%     0. initialize experiment and analysis parameters
%     1. isolate timepoint, curve i.d., and added mass since birth
%     2. use curve i.d. and added mass to find magnitude of incremental mass additions
%     3. accumulate data points (of cell cycle stage) by period fraction
%     4. calculate mean cell cycle stage per time bin
%     5. plot!
     


for condition = 1:2;                                  % for looping between conditions
    
    clearvars -except condition;
    
    % 0. initialize time binning parameters
    periodDuration = 1;                             % duration of nutrient period in hours
    binsPerHour = 20;                               % self-explanatory
    hrPerBin = 1/binsPerHour;                       % time bins of 0.005 hr
    
    % 0. initialize time vector for plotting
    binsPerPeriod = periodDuration/hrPerBin;
    periodTime = linspace(1, binsPerPeriod, binsPerPeriod);
    periodTime = hrPerBin*periodTime';
    periodFraction = periodTime/periodDuration;
    
    % 0. initialize time window parameters for analysis
    firstHour = 5;                                  % time at which to initate analysis
    finalHour = 10;                                 % time at which to terminate analysis

    % 0a. Load data matrices
    
    dmDirectory = dir('dm*.mat'); % note: this assumes the only two data matrices are 'const' and 'fluc'
    names = {dmDirectory.name}; % loaded alphabetically
    
    for dm = 1:length(names)
        load(names{dm});
        dataMatrices{dm} = dataMatrix;                                         % for entire condition
    end
    clear dataMatrix dmDirectory dm;
    clear names;

    interestingData = dataMatrices{condition};      % condition: 1 = constant, 2 = fluctuating
    
    % 1. isolate time and massAdded data
    timeStamps = interestingData(:,2);
    curveID = interestingData(:,6);
    massAdded = interestingData(:,10);              % mass added since birth
    
    %  trim off timepoints earlier than first
    curveID = curveID(timeStamps >= firstHour);
    massAdded = massAdded(timeStamps >= firstHour);
    lowTrimmed_timeStamps = timeStamps(timeStamps >= firstHour);
    
    %  trim off timepoints later than last
    curveID = curveID(lowTrimmed_timeStamps <= finalHour);
    massAdded = massAdded(lowTrimmed_timeStamps <= finalHour);
    finalTrimmed_timeStamps = lowTrimmed_timeStamps(lowTrimmed_timeStamps <= finalHour);
    
   
    
    % 2. calculate instantaneous added mass  
    %    strategy: add differences from individual curves to a vector of
    %    zeros, such that indexing is preserved
    
    incrementAdded = zeros(length(curveID),1);
    for id = 1:max(curveID)
        
        %  replace all values NOT part of the ID with zero
        currentCurves = curveID;
        currentCurves(currentCurves~=id)=0;     % keep values only for indices with current i.d.
        currentCurves = currentCurves/id;       % mask of ones and zeros
        currentAddedMSB = massAdded.*currentCurves;   % multiply so that all non i.d. indices become zero
        
        %  calcalute differences
        currentDiffs = diff(currentAddedMSB);
        currentDiffs = [0; currentDiffs];
        
        %  replace all negatives (end of curve, noise) with zero
        currentDiffs(currentDiffs<0)=0;
        
        %  add current differences to final accumulation of added increments
        incrementAdded = incrementAdded + currentDiffs;
    end
    
    %  eliminate non-complete curves from analysis
    incrementAdded(curveID == 0) = NaN;
    
    
    % 3. accumulate data by associated period fraction
    currentHours = floor(finalTrimmed_timeStamps);
    periodFractions = finalTrimmed_timeStamps-currentHours;
    fractionBins = ceil(periodFractions*binsPerPeriod);
    binnedByPeriodFraction = accumarray(fractionBins,incrementAdded,[],@(x) {x});
    
    % add length to binned vector, if needed
    if length(binnedByPeriodFraction) < binsPerPeriod
        binnedByPeriodFraction{binsPerPeriod} = [];
    end
    
    
    % 4. calculate mean cell cycle stage per time bin within current period
    meanIncrement = cell2mat( cellfun(@nanmean,binnedByPeriodFraction,'UniformOutput',false) );
    stdIncrement = cell2mat( cellfun(@nanstd,binnedByPeriodFraction,'UniformOutput',false) );
    nIncrement = cell2mat( cellfun(@length,binnedByPeriodFraction,'UniformOutput',false) );
    errorIncrement = stdIncrement./sqrt(nIncrement);
    
    
    % 5. plot mean for current period
    %    ...but before plotting, remove nans from cell cycle and time data
    eventMask = find(~isnan(meanIncrement));                                       % find indices of cell cycle data
    meanIncrement = meanIncrement(eventMask);                                      % trim all nans from ccStage vector
    stdIncrement = stdIncrement(eventMask);
    errorIncrement = errorIncrement(eventMask);
    periodFraction = periodFraction(eventMask);
    
    %  !
    %  cheating..... to fix within matrixBuilder
    meanIncrement(1) = meanIncrement(end);
    errorIncrement(1) = errorIncrement(end);
    meanIncrement(2) = meanIncrement(end);
    errorIncrement(2) = errorIncrement(end);
    
    figure(2)
    if condition == 1
        plot(periodFraction,meanIncrement,'k')   % constant
        axis([0,1.05,0,.2])
        grid on
        hold on
        errorbar(periodFraction,meanIncrement,errorIncrement,'k')
    else
        plot(periodFraction,meanIncrement,'b')   % fluctuating (blue)
        hold on
        errorbar(periodFraction,meanIncrement,errorIncrement,'b')
    end
    
end




%%   T H R E E.
%    Instantaneous added mass over a nutrient period (per period) 


%    Strategy:
%
%     0. initialize experiment and analysis parameters
%     1. isolate data of interest
%     2. accumulate data points (of cell cycle stage) by time bin
%     3. isolate current period of interest
%     4. calculate mean cell cycle stage over timsteps of current period
%     5. loops through all periods, to overlay each xy overlaid on a single plot!
%     

%  The data matrix this script arranges and plots:
%
%           columns: periods of interest, start with first full period after growth rate equilibrates 
%              rows: each timestep in period time (period fraction)



% 0. initialize time binning parameters
periodDuration = .25;                             % duration of nutrient period in hours                 
binsPerHour = 20;                               % self-explanatory 
hrPerBin = 1/binsPerHour;                       % time bins of 0.05 hr

% 0. initialize time vector for plotting
binsPerPeriod = periodDuration/hrPerBin;
periodTime = linspace(1, binsPerPeriod, binsPerPeriod);
periodTime = hrPerBin*periodTime';                                       

% 0. initialize looping parameters for analysis
firstHour = 5;                                  % time at which to initate analysis
finalHour = 10;                                 % time at which to terminate analysis
firstTimepoint = firstHour*binsPerHour + 1;     % calculate first timepoint (row number in binnedByTime)
numPeriods = (finalHour-firstHour)/periodDuration;


for condition = 1:2;                                  % for looping between conditions
    
    interestingData = dataMatrices{condition};      % condition: 1 = constant, 2 = fluctuating
    currentTimepoint = firstTimepoint;              % initialize first timepoint as current timepoint
    
    % 1a. isolate time and massAdded data
    timeStamps = interestingData(:,2);
    massAdded = interestingData(:,10);              % mass added since birth
    
    % 1b.  calculate instantaneous added mass
    incrementAdded = diff(massAdded);
    incrementAdded = [0; incrementAdded];
    
    % 1c.  eliminate zeros (non-full track data and births) and negatives (divisions and noise) from data
    massAdded(massAdded <= 0) = NaN;
    incrementAdded(incrementAdded <= 0) = NaN;
    
    % 2. accumulate data by associated time bin
    timeBins = ceil(timeStamps*binsPerHour);
    binnedByTime_increments = accumarray(timeBins,incrementAdded,[],@(x) {x});
    
    for period = 1:numPeriods
        
        % 3. establish current period in loop
        currentPeriod = currentTimepoint:(currentTimepoint + binsPerPeriod -1);
        
        if period < numPeriods
            currentIncrements = binnedByTime_increments(currentPeriod);
        else
            currentIncrements = binnedByTime_increments(currentTimepoint:end);
        end
        currentTimepoint = currentTimepoint + binsPerPeriod; % re-define currentTimepoint for next loop cycle
        
        % 4. calculate mean cell cycle stage per time bin within current period
        meanIncrement = cell2mat( cellfun(@nanmean,currentIncrements,'UniformOutput',false) );
        stdIncrement = cell2mat( cellfun(@nanstd,currentIncrements,'UniformOutput',false) );
        nIncrement = cell2mat( cellfun(@length,currentIncrements,'UniformOutput',false) );
        errorIncrement = stdIncrement./sqrt(nIncrement);
        
        % 5. plot mean for current period
        %    ...but before plotting, remove nans from cell cycle and time data
        eventMask = find(~isnan(meanIncrement));                                       % find indices of cell cycle data
        meanIncrement = meanIncrement(eventMask);                                      % trim all nans from ccStage vector
        stdIncrement = stdIncrement(eventMask);
        errorIncrement = errorIncrement(eventMask);
        
        fractionVector = periodTime;
        fractionVector = fractionVector(eventMask);                                          % trim all nans from time vector
        currentFraction = fractionVector./periodDuration;
        
        figure(1)
        if condition == 1
            plot(currentFraction,meanIncrement,'color',[0,0,0]+(1/period)*[1,1,1],'linewidth',1.01)   % constant
            axis([0,1,0,.1])
            grid on
            hold on
            errorbar(currentFraction,meanIncrement,errorIncrement,'k')
        else
            plot(currentFraction,meanIncrement,'color',[0.2,0,.2]+(1/period)*[0.7,0.7,0],'linewidth',1.01)   % fluctuating (blue)
            hold on
            errorbar(currentFraction,meanIncrement,errorIncrement,'b')
        end
        
    end
    clear period

end











