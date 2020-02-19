%%  cell cycle fraction over time


%  Goal: Searching for synchrony in growth data.
%        This script plots up multiple views of cell cycle stage:
%        
%        1. mean and std over time
%        2. heatmap over time


%  Last edit: Jen Nguyen, April 6th 2017




% Growth phase is defined as a specific fraction of the growth curve, as
% calculated and assembled in matrixBuilder.

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
%   Initialize.

dmDirectory = dir('dm*.mat'); % note: this assumes the only two data matrices are 'const' and 'fluc'
names = {dmDirectory.name}; % loaded alphabetically

for dm = 1:length(names)
    load(names{dm});                
    dataMatrices{dm} = dataMatrix;                                         % for entire condition
end                                                                        

clear dataMatrix dmDirectory dm;
clear names;


%%   O N E.
%    Mean cell cycle stage over time 


%  Goal: does cell cycle stage oscillate over time?
%        less standard deviation around mean suggests greater synchrony


%  Strategy:
%
%     0. initialize experiment and analysis parameters
%     1. isolate data of interest
%     2. accumulate data points (of cell cycle stage) by time bin
%     3. calculate mean cell cycle stage per time bin
%     4. plot!
%     

expHours = 11; %  duration of experiment in hours                          % 0.  initialize parameters
binFactor = 20; % time bins of 0.05 hr  
hrPerBin = 1/binFactor; 
counter = 1;

for condition = 5
   
    % in fluctuation experiments:
    %       1 = constant
    %       2 = fluctuating
    
    % in poly-K challenge 0316:
    %       1 = no treatment, zero glucose
    %       2 = no treatment, 1 uM glucose
    %       3 = no treatment, 1 uM glucose with gentle inoculation
    %       4 = poly-lysine,  zero glucose
    %       5 = poly-lysine,  1 uM glucose
    %       6 = poly-lysine,  1 uM glucose
    
    
    % 1.  isolate time and ccStage vectors
    interestingData = dataMatrices{condition};  % condition: 1 = constant, 2 = fluctuating
    timeStamps = interestingData(:,2);
    ccStage = interestingData(:,9);                                        
    
    % 2.  accumulate data by associated time bin
    timeBins = ceil(timeStamps*binFactor);                                 
    binnedByTime = accumarray(timeBins,ccStage,[],@(x) {x});               
    
    % 3a.  calculate mean cell cycle stage per bin
    meanStage = cellfun(@nanmean,binnedByTime,'UniformOutput',false);      
    meanStage = cell2mat(meanStage);
    
    % 3b.  calculate std and error
    devStage = cell2mat(cellfun(@nanstd,binnedByTime,'UniformOutput',false));
    nStage = cell2mat(cellfun(@length,binnedByTime,'UniformOutput',false));
    errorStage = devStage./sqrt(nStage);
    
    % 4a. create a time vector to convert bin # to absolute time
    timeVector = linspace(1, expHours/hrPerBin, expHours/hrPerBin);
    timeVector = hrPerBin*timeVector';                                       
    
    % 4b. before plotting, a little matrix manipulation to get around nans                                                            
    eventMask = find(~isnan(meanStage));                                   % find indices of cell cycle data
    meanStage = meanStage(eventMask);                                      % trim all nans from ccStage vector
    devStage = devStage(eventMask);
    errorStage = errorStage(eventMask);
    timeVector = timeVector(eventMask);                                    % trim all nans from time vector
    
    % 4c. mean with standard deviation
    figure(1)
    if counter == 1
         plot(timeVector,meanStage,'k')
         axis([0,11,0,1])
         hold on
         grid on
         errorbar(timeVector,meanStage, devStage,'k')
    else
        plot(timeVector,meanStage,'b')
        hold on
        errorbar(timeVector,meanStage,devStage,'b')
    end
    
    
    % 4d. mean with standard error
    figure(2)
    if counter == 1
        plot(timeVector,meanStage,'k')
        axis([0,11,0,1])
        hold on
        grid on
        errorbar(timeVector,meanStage,errorStage,'k')
    else
        plot(timeVector,meanStage,'b')
        hold on
        errorbar(timeVector,meanStage,errorStage,'b')
    end
    counter = counter + 1;
    
    % 4c. mean with standard error for all, without color designation
    figure(3)
    plot(timeVector,meanStage)
    axis([0,11,0,1])
    hold on
    grid on
    errorbar(timeVector,meanStage,errorStage)
    
end
legend('condition 1', 'condition 2', 'condition 3', 'condition 4', 'condition 5', 'condition 6');

%%   T W O.
%    Heatmap of cell cycle stage over time


%  Goal: display fraction of population in each cell cycle stage vs. time
%        striations suggest synchrony, whereas uniformity indicates
%        heterogeneity


%  Strategy:
%
%     0. isolate data of interest
%     1. define bin sizes: time and cell cycle stage
%     2. accumulate data points (of cell cycle stage) by time bin
%     3. spread time binned data into vertical cell cycle stage bins
%     4. count number of points in each bin
%     5. with counts, generate normalized plot!


for condition = 1:2
    
    interestingData = dataMatrices{condition};  % condition: 1 = constant, 2 = fluctuating
    timestamps = interestingData(:,2);
    ccStage = interestingData(:,9);                                            % 0.  isolate time and ccStage vectors
    
    timeBins = ceil(timestamps*200);  % time bins of 0.005 hr                        % 1a. define bin size (time)
    binnedByTime = accumarray(timeBins,ccStage,[],@(x) {x});                   % 2.  accumulate data by associated time bin
    
    
    
    
    % A. Generate grid of absolute counts
    
    dataGrid = zeros(10,2000);                                                 % (rows) ccStage: 0 - 1, with 0.1 incr.
    % (columns) time: 0 - 10, with 0.005 incr
    for i = 1:length(binnedByTime)
        if isempty(binnedByTime{i})
            continue
        else
            stageBins = ceil(binnedByTime{i}/.1);                              % 1b. define bin size (cell cycle stage)
            stageBins(stageBins==0) = 1;                                       %     manually include birth into 1st bin
            
            currentTimeStep = binnedByTime{i};                                 % 3.  accumulate ccStage into bins
            binnedByStage = accumarray(stageBins(~isnan(stageBins)),currentTimeStep(~isnan(currentTimeStep)),[],@(x){x});
            
            if isempty(binnedByStage)                                          % 4.  some timepoints are empty vectors
                continue                                                       %     due to non-counting of NaNs.
            else                                                               %     by-pass these empty vectors!
                counts = cellfun(@length,binnedByStage,'UniformOutput',false);
                counts = cell2mat(counts);
            end
            
            dataGrid(1:length(counts),i) = counts;
        end
    end
    clear i currentTimeStep binnedByStage counts;
    
    
    
    % B. Normalized grid
    
    countsPerTimePoint = sum(dataGrid);                                        % 5.  find total counts per timepoint
    normalizedByCount = zeros(10,2000);                                        %     divide each bin count by total
    
    for ii = 1:length(dataGrid)
        if countsPerTimePoint(ii) > 0
            normalizedByCount(:,ii) = dataGrid(:,ii)./countsPerTimePoint(ii);
        else
            continue
        end
    end
    
    
    
    % C. Plot!!
    
    figure(1)
    subplot(2,1,condition)
    imagesc(normalizedByCount,[0.02 .1])
    axis xy
    axis([0,2000,1,10])
    colorbar
    hold on
    
end










