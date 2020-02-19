% Incremental volume added vs. Stage of period

%  Last edit: Jen Nguyen, May 19th 2017


%  Goal: plot instantaneous V/time data per period fraction

    %   0.  Load data (measurements and meta data)
    %   1.  Isolate data from condition of interest
    %   2.  Initialize time data, determine time window of interest
    %   3.  Trim data to time windows of interest
    %   4.  Bin data by designated period fraction


% use data matrix to plot instantaneous / incremental increases in volume
% over a single period

%%
clear

% 0. Load data matrix
load('dm-t900-2016-10-20.mat');

% 0. Initialize meta data
load('meta.mat');
meta = meta_2016oct20;

% 0. Meta data for period fraction
periodLength = 900;                         % in seconds
binsPerPeriod = 20;

%%
for condition = 1%:2:3
    %%
    % 1. isolate data from condition of interest
    interestingData = find(dataMatrix(:,28) == condition);
    conditionData = dataMatrix(interestingData,:);
    clear interestingData;
    
    % 1. initialize time data
    %interestingTime = find(dataMatrix(:,2) == condition);
    timeTrack = conditionData(:,2); % in hours
    
    % 1. designate time window of analysis
    firstTimepoint = meta(condition,3); % in hours
    lastTimepoint = meta(condition,4);
    
    % 2. isolate incremental volume data
    incrementTrack = conditionData(:,27); % volume as approximated by capped cylinder

    % 3. remove data outside time window of interest
    timeTrack(timeTrack < firstTimepoint) = NaN;
    timeTrack(timeTrack > lastTimepoint) = NaN;
    
    timeFilter = find(~isnan(timeTrack));
    incrementTrack = incrementTrack(timeFilter);
    timeTrack = timeTrack(timeFilter);
    
    % 3. remove data with negative changes in volume
%     trimmedMu = incrementTrack;
%     trimmedMu(trimmedMu > 1) = NaN;
%     trimmedMu(trimmedMu <= 0) = NaN;
%     
%     % remove NaNs from data sets
%     nanFilter = find(~isnan(trimmedMu));
%     trimmedMu = trimmedMu(nanFilter);
%     trimmedTime = rightBin(nanFilter);

    % 4. create vector of bin assignments based on time
    timeInPeriods = (timeTrack*3600)/periodLength; % units = seconds/seconds
    timeInPeriods_floored = floor(timeInPeriods);
    timeAsPeriodFraction = timeInPeriods - timeInPeriods_floored;
    periodFractions = timeAsPeriodFraction * binsPerPeriod;
    periodFractions_binned = ceil(periodFractions);
      
    % 5. bin incremental volume data
    bin_mean = accumarray(periodFractions_binned,incrementTrack,[],@nanmean);
    bin_std = accumarray(periodFractions_binned,incrementTrack,[],@nanstd);
    
    % 5. for s.e.m. per bin
    countsPerBin = zeros(binsPerPeriod,1);
    for j = 1:binsPerPeriod
        countsPerBin(j) = length(find(periodFractions_binned == j));
    end
    clear j
    
    bin_sem = bin_std./sqrt(countsPerBin);

%%
    % 6. plot
    figure()
    errorbar(bin_mean,bin_sem)
    hold on
    grid on
    axis([-0.2,20.2,-0.008,.012])
    xlabel('Time')
    ylabel('incremental added volume, VA')

end

