%% rightTimeLine


% Goal: Searching for any influence of nutrient environment on the cell cycle.

%       Like nsyncNFlux.m, this script follows cell cycle behavoir over a single nutrient period.
%       Key difference: this script plots the average cell cycle duration of
%                       cells born within a given increment of the nutrient period.


%  Last edit: Jen Nguyen, February 26th 2016




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



%  Strategy:
%
%     0. initialize experimental and analytical parameters
%     1. isolate data of interest
%     2. identify rows where drop = 1 AND timeSinceBirth = 0 AND curveDuration > 0, together these mark unique tracks      
%     3. find curveDuration and timestamp at those rows
%     4. trim all undesired timestamps (i.e. keep only those in final hrs)
%     5. assign all timestamps a period fraction
%     6. bin and average curveDurations by period fraction
%     7. plot !


% OK! Lez go!

%%

% Initialize data.

dmDirectory = dir('dm*.mat'); % note: this assumes the only two data matrices are 'const' and 'fluc'
names = {dmDirectory.name}; % loaded alphabetically

for dm = 1:length(names)
    load(names{dm});                
    dataMatrices{dm} = dataMatrix;
end                                                                        

clear dataMatrix dmDirectory dm;
clear names;


% 0. initialize experimental parameters
expHours = 10; %  duration of experiment in hours                      

% 0. initialize time binning parameters
periodDuration = 1; %0.25;                          % duration of nutrient period in hours                 
binsPerHour = 200;                                  % bPH of 200 = time bins of 0.005 hours (18 sec)
hrPerBin = 1/binsPerHour;                           % bPH of 40  = time bins of 0.025 hours (1.5 min)

% 0. initialize time vector for plotting
binsPerPeriod = periodDuration/hrPerBin;
periodTime = linspace(1, binsPerPeriod, binsPerPeriod);
periodTime = hrPerBin*periodTime';                                       

% 0. initialize looping parameters for analysis
firstHour = 5;                                      % time at which to initate analysis
finalHour = 10;                                     % time at which to terminate analysis
%firstTimepoint = firstHour*binsPerHour + 1;         % calculate first timepoint (row number in binnedByTime)
numPeriods = (finalHour-firstHour)/periodDuration;  % number of periods of interest
totalPeriods = finalHour/periodDuration;            % total periods in experiment


for condition = 1:2;                                  % for looping between conditions
 
    % 1. isolate data of interest
    interestingData = dataMatrices{condition};      % condition: 1 = constant, 2 = fluctuating
    %currentTimepoint = firstTimepoint;              % initialize first timepoint as current timepoint
    
    timeStamps = interestingData(:,2);
    drop = interestingData(:,5);                    % col #5 = drop boolean (1 = birth event)      
    timeSinceBirth = interestingData(:,7);          % col #7 = time since birth
    curveDurations = interestingData(:,8);          % col #8 = curve duration

    
    % 2. identify rows with unique tracks:  drop = 1  AND  timeSinceBirth = 0  AND  curveDuration > 0 
    durations = [];
    times = [];
    for r = 1:length(drop)
        if (drop(r) == 1) && (timeSinceBirth(r) == 0) && (curveDurations(r) > 0)
            % 3. find curveDuration and timestamp at those rows
            durations = [durations; curveDurations(r)];
            times = [times; timeStamps(r)];
        else
            continue
        end
    end
    clear r;
    
    % 4. trim off data points prior to firstHour
    keepTheseRows = find(times > firstHour);
    times = times(keepTheseRows);
    durations = durations(keepTheseRows);
   
    % 5a. assign all timestamps a period fraction
    timeBins = ceil(times*binsPerHour);
    
    % 6a. bin curveDurations by period fraction
    binnedByTime = accumarray(timeBins,durations,[],@(x) {x});
    binnedByTime{expHours/hrPerBin,1} = [];

    % 6b. average data points per time bin (and count, for normalization downstream)
    notNormalized = cell2mat( cellfun(@nanmean,binnedByTime,'UniformOutput',false) );
    countsPerTimeBin = cell2mat( cellfun(@length,binnedByTime,'UniformOutput',false) );
    
    % 6c. associate time bin with appropriate period fraction
    periodFraction = [];
    for p = 1:totalPeriods
        periodFraction = [periodFraction; periodTime];
    end
    clear p;
    fractionAssignments = floor(periodFraction.*binsPerHour);          % accumarray requires integers
    binnedByPeriod = accumarray(fractionAssignments,notNormalized,[],@(x) {x});
    countsByPeriod = accumarray(fractionAssignments,countsPerTimeBin,[],@(x) {x});

    % 6d. average accumulated averages, with weighting based on initial counts
    normalizedDuration = zeros(binsPerPeriod,1);
    for f = 1:binsPerPeriod
        if isempty(binnedByPeriod{f})
            continue
        else
            averages = binnedByPeriod{f};
            counts = sum(countsByPeriod{f});
            weights = countsByPeriod{f}/counts;
            normalized = averages.*weights;
            normalizedDuration(f) = nansum(normalized);
        end
    end
    clear f;
    
    
    % 7. LINE PLOT: mean over period fraction
    normalizedDuration(normalizedDuration==0) = NaN;                   % ...but before plotting, remove zeros and nans !
    eventMask = find(~isnan(normalizedDuration));                      % find indices of cell cycle data
    normalizedDuration = normalizedDuration(eventMask);                % trim all nans from ccStage vector
    currentTime = periodTime;
    currentTime = currentTime(eventMask);                              % trim all nans from time vector
    currentFraction = currentTime./periodDuration;
    
    figure(2)
    if condition == 1
        plot(currentFraction,normalizedDuration,'k')   % constant
        axis([0,1,0,5])
        grid on
        hold on
    else
        plot(currentFraction,normalizedDuration,'b')   % fluctuating (blue)
    end
    

end

