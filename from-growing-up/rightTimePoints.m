%% rightTimePoints


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
%     6. plot !


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
                    

% 0. initialize time binning parameters
periodDuration = 1; %0.25;                          % duration of nutrient period in hours                 
binsPerHour = 200;                                  % bPH of 200 = time bins of 0.005 hours (18 sec)
hrPerBin = 1/binsPerHour;                           % bPH of 40  = time bins of 0.025 hours (1.5 min)
                                  

% 0. initialize looping parameters for analysis
firstHour = 5;                                      % time at which to initate analysis


for condition = 1:2;                                  % for looping between conditions
 
    % 1. isolate data of interest
    interestingData = dataMatrices{condition};      % condition: 1 = constant, 2 = fluctuating
    
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

    % 4. trim off data points pre-equilibration
    keepTheseRows = find(times > firstHour);
    times = times(keepTheseRows);
    durations = durations(keepTheseRows);
    
    % 5. assign all timestamps a period fraction
    toSubtract = floor(times);  
    onlyFractions = times-toSubtract;
    
    % 6. plot
    figure(1)
    subplot(2,1,condition)
    if condition == 1
        scatter(onlyFractions,durations,'k')
        axis([0,1,0,5])
        grid on
    else
        scatter(onlyFractions,durations,'b')
        axis([0,1,0,5])
        grid on
    end
    
   

end

