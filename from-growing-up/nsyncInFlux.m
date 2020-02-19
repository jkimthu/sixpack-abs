%% nsyncInFlux


% Goal: Searching for synchrony in growth data.

%       This script plots the population mean of cell cycle stage
%	    over a single fluctuating period    


%  Last edit: Jen Nguyen, February 19th 2016





% Growth phase is defined as a specific fraction of the growth curve, as
% calculated and assembled in matrixBuilder.
%
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

%%

% Initialize data.

dmDirectory = dir('dm*.mat'); % note: this assumes the only two data matrices are 'const' and 'fluc'
names = {dmDirectory.name}; % loaded alphabetically

for dm = 1:length(names)
    load(names{dm});                
    dataMatrices{dm} = dataMatrix;                                         % for entire condition
end                                                                        

clear dataMatrix dmDirectory dm;
clear names;

%

%  Line maps: average cell cycle stage per xy vs. time
%
%     -  recreating figure 3 of Mathis & Ackermann pre-print: line plots separate
%        experiments to illustrate reduced variation after pulsed shock



%  Strategy:
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
periodDuration = 0.25;                             % duration of nutrient period in hours                 
binsPerHour = 200;                              % self-explanatory 
hrPerBin = 1/binsPerHour;                       % time bins of 0.005 hr

% 0. initialize time vector for plotting
binsPerPeriod = periodDuration/hrPerBin;
periodTime = linspace(1, binsPerPeriod, binsPerPeriod);
periodTime = hrPerBin*periodTime';                                       

% 0. initialize looping parameters for analysis
firstHour = 5;                                  % time at which to initate analysis
finalHour = 10;                                 % time at which to terminate analysis
firstTimepoint = firstHour*binsPerHour + 1;     % calculate first timepoint (row number in binnedByTime)
numPeriods = (finalHour-firstHour)/periodDuration;

%
for condition = 1:2;                                  % for looping between conditions
    
    interestingData = dataMatrices{condition};      % condition: 1 = constant, 2 = fluctuating
    currentTimepoint = firstTimepoint;              % initialize first timepoint as current timepoint
    
    % 1. isolate time and ccStage data
    timeStamps = interestingData(:,2);
    ccStage = interestingData(:,9);
    
    % 2. accumulate data by associated time bin
    timeBins = ceil(timeStamps*binsPerHour);
    binnedByTime = accumarray(timeBins,ccStage,[],@(x) {x});
    
    for period = 1:numPeriods
        
        % 3. establish current period in loop
        currentPeriod = currentTimepoint:(currentTimepoint + binsPerPeriod -1);
        if period < numPeriods
            currentStages = binnedByTime(currentPeriod);
        else
            currentStages = binnedByTime(currentTimepoint:end);
        end
        currentTimepoint = currentTimepoint + binsPerPeriod; % re-define currentTimepoint for next loop cycle
        
        % 4. calculate mean cell cycle stage per time bin within current period
        meanStage = cell2mat( cellfun(@nanmean,currentStages,'UniformOutput',false) );
        stdStage = cell2mat( cellfun(@nanstd,currentStages,'UniformOutput',false) );
        nStage = cell2mat( cellfun(@length,currentStages,'UniformOutput',false) );
        errorStage = stdStage./sqrt(nStage);
        
        % 5. plot mean for current period
        %    ...but before plotting, remove nans from cell cycle and time data
        eventMask = find(~isnan(meanStage));                                       % find indices of cell cycle data
        meanStage = meanStage(eventMask);                                          % trim all nans from ccStage vector
        currentTime = periodTime;
        currentTime = currentTime(eventMask);                                          % trim all nans from time vector
        currentFraction = currentTime./periodDuration;
        
        figure(1)
        if condition == 1
            plot(currentFraction,meanStage,'color',[0,0,0]+(1/period)*[1,1,1],'linewidth',1.01)   % constant
            axis([0,1,0,1])
            grid on
            hold on
        else
            plot(currentFraction,meanStage,'color',[0.1,0,1]+(1/period)*[0.2,0.6,0],'linewidth',1.01)   % fluctuating (blue)
            hold on
        end
        
    end
    clear period

end

