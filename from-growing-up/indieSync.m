%%  indieSYNC


%  Goal: Searching for synchrony in growth data.
%        Like nSync.m, this script plots up multiple views of cell cycle stage.
%        Unlike nSync.m, indieSync works with isolated data matrices from
%        individual xy positions (in lieu of independent experiments)
%
%  Last edit: Jen Nguyen, February 17th 2016




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


%        col        1         2        3         4         5          6              7                8              9 
%       ------------------------------------------------------------------------------------------------------------------
%        row      Track#    Time     Lngth      Mu       drop?      curve#    timeSinceBirth    curveDuration    cc stage
%       ------------------------------------------------------------------------------------------------------------------
%         1         1         t        x         u         0*         1              0                3              0
%         2         1         t        x         u         0          1              1                3             .5
%         3         1         t        x         u         0          1              2                3              1
%         4         1         t        x         u         1          2              0                3              0
%         5         1         t        x         u         0          2              1                3             .5 
%         6         1         t        x         u         0          2              2                3              1
%         7         1         t        x         u         1          3              0                3              0
%         8         1         t        x         u         0          3              1                3             .5
%         9         1         t        x         u         0          3              2                3              1
%         10        1         t        x         u         1          4              0                3             nan


%       where,
%                row     =  row number, obvi
%                t       =  all timepoints associated with concatinated length trajectories
%                x       =  length values from concatentated length trajectories
%                mu      =  calculated growth rates from SlidingFits.m
%                drop?   =  finding where individual cell cycles start and end, a boolean
%                curve   =  an id number for each individual cell cycle
%                stage   =  time since birth / duration of entire cycle


% Strategy:
%
%       1.  for each curve, determine duration (time)
%       2.  for each time step, determine absolute time since birth
%       3.  for each data point in vector, record as fraction:
%                
%               ccStage = time since birth / total curve duration



% Considerations:
%
%       1. Does separation between phase-sorted subpopulations occur?
%       2. Vary number of fractions. Which leads to the best separation?
%       3. If there is separation, what explains it?



% OK! Lez go!

%%
%   Initialize data.

dmDirectory = dir('dm0818-xy*.mat');
names = {dmDirectory.name}; % loaded alphabetically

for dm = 1:length(names)
    load(names{dm});                                             
    dataMatrices{dm} = indivDM; % for individual positions
end                                                                        

clear dataMatrix dmDirectory dm;
clear names;


%%


%      O N E.  Cell cycle stage over experiment



%  Line maps: average cell cycle stage per xy vs. time
%
%     -  recreating figure 3 of Mathis & Ackermann pre-print: line plots separate
%        experiments to illustrate reduced variation after pulsed shock



%  Strategy:
%
%     0. initialize experiment and analysis parameters
%     1. isolate data of interest
%     2. accumulate data points (of cell cycle stage) by time bin
%     3. calculate mean cell cycle stage per time bin
%     4. overlay each xy overlaid on a single plot!
%     


                                 
expHours = 10; %  duration of experiment in hours                          % 0.  initialize parameters
binFactor = 200; % time bins of 0.005 hr  
hrPerBin = 1/binFactor; 

for xy = 1:12
%  
    interestingData = dataMatrices{xy};  % condition: 1 = constant, 2 = fluctuating
    timeStamps = interestingData(:,2);
    ccStage = interestingData(:,9);                                        % 1.  isolate time and ccStage vectors
    

    timeBins = ceil(timeStamps*binFactor);                                 % 2.  accumulate data by associated time bin
    binnedByTime = accumarray(timeBins,ccStage,[],@(x) {x});               
    meanStage = cellfun(@nanmean,binnedByTime,'UniformOutput',false);      % 3.  calculate mean cell cycle stage per bin
    meanStage = cell2mat(meanStage);   
    
    indieTime = linspace(1, expHours/hrPerBin, expHours/hrPerBin);         % as exact timestamps vary between xy positions,
    indieTime = hrPerBin*indieTime';                                       % create a time vector to convert bin # to absolute time
    
    % before plotting, a little matrix manipulation to get around nans
                                                                           
    eventMask = find(~isnan(meanStage)); % find indices of cell cycle data
    meanStage = meanStage(eventMask); % trim all nans from ccStage vector
    indieTime = indieTime(eventMask); % trim all nans from time vector
    
    
    figure(1)
    %plot(indieTime,meanStage,'color',[0,0,0]+(xy)*[.015,.015,.015])          % constant
    plot(indieTime,meanStage,'color',[0,0,0]+(xy)*[.01,.04,.05])          % fluctuating
    %plot(indieTime,meanStage)
    hold on
    
end

%%




%%

% B R A I N S T O R M


% 1. Cell cycle fraction over time
%
%           i. line (average)
%          ii. scatter
%         iii. heatmap
%          iv. percent dividing

% 2. Cycle durations over time 
%
%
%








