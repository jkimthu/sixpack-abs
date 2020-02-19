%%  cell cycle fraction over time


%  Goal: Searching for synchrony in growth data.
%        This script plots up multiple views of cell cycle stage.
%
%  Last edit: Jen Nguyen, March 27th 2016




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


%%


%      O N E.  Cell cycle stage over experiment



%  Heatmap: fraction of population in each cell cycle stage vs. time
%
%     -  figure 3 of Mathis & Ackermann pre-print: line plots separate
%        experiments to illustrate reduced variation after pulsed shock
%     -  here, let's plot all data points per timestep


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










