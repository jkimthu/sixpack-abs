%% mew-mew cycle


% Goal: Is the instantaneous growth rate a function of current cell cycle stage?

%       This script accumulates mu's into increments of cell cycle fraction,
%       and plots a mean line with error. 
%
%       x axis: mean cell cycle stage
%       y axis: mass added since birth


%  Last edit: Jen Nguyen, March 15th 2016



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

clear
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
%     2. determine bin size for cell cycle fraction
%     3. accumulate growth rate based on cc bins
%     4. calculate mean, std, n, error
%     5. plot!



% 0. designate time window of analysis

firstTimepoint = 5; % in hours
lastTimepoint = 10;



for condition = 1:2          % 1 = constant, 2 = fluctuating
    
    interestingData = dataMatrices{condition};
    
    % 1. isolate addedMass and ccStage data
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
    
    % replace all values of Mu <= 0 and Mu > 1 with NaN
    trimmedMu = Mu;
    trimmedMu(trimmedMu > 1) = NaN;
    trimmedMu(trimmedMu <= 0) = NaN;
    %histogram(trimmedMu)
    
    % remove NaNs from data sets
    nanFilter = find(~isnan(trimmedMu));
    trimmedMu = trimmedMu(nanFilter);
    trimmedStages = ccStage(nanFilter);
    
    % replace all values of trimmedStages = 0 with NaN
    trimmedStages(trimmedStages == 0) = NaN;
    
    % remove NaNs from data sets
    zeroFilter = find(~isnan(trimmedStages)); % lists all indeces that are NOT NaN
    trimmedMu = trimmedMu(zeroFilter);
    trimmedStages = trimmedStages(zeroFilter);
    
    % create binning vector such that bin size = 0.05 or 1/20th of cell cycle
    stageBins = ceil(trimmedStages*20);
    
    
    % 3. accumulate ccStage by binned growth rates
    binnedByStage = accumarray(stageBins,trimmedMu,[],@(x) {x});
    
    
    % 4. calculate mean, std, n, and error
    meanMu = cell2mat( cellfun(@nanmean,binnedByStage,'UniformOutput',false) );
    stdMu = cell2mat( cellfun(@nanstd,binnedByStage,'UniformOutput',false) );
    nMu = cell2mat( cellfun(@length,binnedByStage,'UniformOutput',false) );
    errorMass = stdMu./sqrt(nMu);
    
    
    % 5. plot mean and std
    figure(1)
    if condition == 1
        plot(meanMu,'k')
        axis([0,21,0.3,0.65])
        hold on
        grid on
        errorbar(meanMu,stdMu,'k')
    else
        plot(meanMu,'b')
        hold on
        errorbar(meanMu,stdMu,'b')
    end
    
    figure(2)
        if condition == 1
        plot(meanMu,'k')
        axis([0,21,0.3,0.65])
        hold on
        grid on
        errorbar(meanMu,errorMass,'k')
    else
        plot(meanMu,'b')
        hold on
        errorbar(meanMu,errorMass,'b')
    end
end


