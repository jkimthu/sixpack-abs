%% cell cycle, mew!


% Goal: Is cell cycle stage a function of instantaneous growth rate?

%       This script accumulates cell cycle fraction data into instantaneous
%       growth rate bins, and plots a mean line with error. 
%
%       x axis: instananeous growth rate
%       y axis: mean cell cycle stage


%  Last edit: Jen Nguyen, March 14th (pi day!) 2016



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

%  Stragety:
%   
%     0. designate time window of analysis
%     1. isolate data of interest (ccStage and mu)
%     2. determine bin size for mu
%     3. accumulate ccStage based on mu bins
%     4. calculate mean, std, n, error
%     5. plot!



% 0. designate time window of analysis

firstTimepoint = 5; % in hours
lastTimepoint = 10;



for condition = 1:2          % 1 = constant, 2 = fluctuating
    
    interestingData = dataMatrices{condition};
    
    % 1. isolate mu and ccStage data
    Mu = interestingData(:,4);
    ccStage = interestingData(:,9);
    timeStamps = interestingData(:,2);
   
    % trim off timepoints earlier than first
    Mu = Mu(timeStamps >= firstTimepoint);
    ccStage = ccStage(timeStamps >= firstTimepoint);
    lowTrimmed_timeStamps = timeStamps(timeStamps >= firstTimepoint);
    
    % trim off timepoints later than last
    Mu = Mu(lowTrimmed_timeStamps <= lastTimepoint);
    ccStage = ccStage(lowTrimmed_timeStamps <= lastTimepoint);
    finalTrimmed_timeStamps = lowTrimmed_timeStamps(lowTrimmed_timeStamps <= lastTimepoint);
    
    % 2. determine bin size for mu
    
    % replace all values of Mu > 1 or Mu <= 0 with NaN
    trimmedMu = Mu;
    trimmedMu(trimmedMu > 1) = NaN;
    trimmedMu(trimmedMu <= 0) = NaN;
    %histogram(trimmedMu)
    
    % remove NaNs from data sets
    nanFilter = find(~isnan(trimmedMu));
    trimmedMu = trimmedMu(nanFilter);
    trimmedStages = ccStage(nanFilter);
    
    % create binning vector such that bin size = 0.05 1/hr
    muBins = ceil(trimmedMu*20);
    
    
    % 3. accumulate ccStage by binned growth rates
    binnedByMu = accumarray(muBins,trimmedStages,[],@(x) {x});
    
    
    % 4. calculate mean, std, n, and error
    meanStage = cell2mat( cellfun(@nanmean,binnedByMu,'UniformOutput',false) );
    stdStage = cell2mat( cellfun(@nanstd,binnedByMu,'UniformOutput',false) );
    nStage = cell2mat( cellfun(@length,binnedByMu,'UniformOutput',false) );
    errorStage = stdStage./sqrt(nStage);
    
    
    % 5. plot mean and std
    figure(1)
    if condition == 1
        plot(meanStage,'k')
        axis([0,21,0,1])
        hold on
        grid on
        errorbar(meanStage,stdStage,'k')
    else
        plot(meanStage,'b')
        hold on
        errorbar(meanStage,stdStage,'b')
    end
    
    figure(2)
        if condition == 1
        plot(meanStage,'k')
        axis([0,21,0,1])
        hold on
        grid on
        errorbar(meanStage,errorStage,'k')
    else
        plot(meanStage,'b')
        hold on
        errorbar(meanStage,errorStage,'b')
    end
end


