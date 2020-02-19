%% isSubpop?

% Goal: this script aims to find subpopulations in our growth experiments

%       the final figure (four panels) plots the distribution unique tracks
%       across bins of doubling rate.
%
%       the four panels sample four distinct time periods on population
%       level observations of doubling rate, and compare distributions
%       between two experimental conditions.
%

% Strategy:
%           0. initialize data and conditions of interest
%           0. initialize time windows of interest
%           1. define time windows of interest
%           2. for each time window...
%                   3. convert timestamps into bins by window
%                   4. isolate data from current window of interest
%                   5. for each condition of interest...
%                           6. collect mus and stats
%                           7. bin mus into mu bins
%                           8. plot pdf onto subplot for time window
%                   9.  repeat for all conditions in window
%          10. repeat for all windows
%

%               


% last edit: jen, 2017 Oct 13

% OK LEZ GO!

%%

% 0. initialize data matrix
clear
load('lb-fluc-2017-10-10-window5-width1p4v1p7-jiggle-0p5-bigger1p8.mat');
dataMatrix = buildDM(D5,M,T);

% 0. designate which conditions from experiment to be analyzed
conditions = [1,2,3]; % of 4 total experiments in 2017-10-10

% 1. define time windows of interest                     
windowLength = 30;      % in minutes
window = [2,4,8,13,20]; % 20 half hr windows total in 10 hr experiment

%%
% 2. for each time window...
for w = 1:length(window);

    % 3. convert timestamps into bins by window
    timestamps = dataMatrix(:,2)/60; % in minutes
    timeBins = ceil(timestamps/windowLength);   % sorted by window length
    
    % 4. isolate data from current window of interest
    currentWindow = dataMatrix(timeBins == window(w),:);
    
    % 5. for each condition...
    for c = 1:3
        
        % 6. collect mus and stats
        col4 = currentWindow(currentWindow(:,35) == c,4);
        mus = col4(col4 > 0);
        
        meanMu{w}(c) = mean(mus);
        stdMu{w}(c) = std(mus);
        countMu{w}(c) = length(mus);
        semMu{w}(c) = std(mus)/sqrt(length(mus));
        
        % 7. bin mus into mu bins
        %bins = linspace(1,70,70);
        binnedMus = ceil(mus*10);
        tracksBinnedByMu = accumarray(binnedMus,mus,[],@(x) {x});
        binCounts = cellfun(@length,tracksBinnedByMu);
        
        % 8. plot pdf onto subplot for time window
        popFractionPerBin = binCounts/length(mus);
        
        figure(1)
        if c == 1
            subplot(1,length(window),w)
            barh(popFractionPerBin,0.4)
            hold on
            title(window(w))
            axis([0,.35,0,50])
            xlabel('pdf')
            ylabel('doubling rate bin')
            
        elseif c == 2
            subplot(1,length(window),w)
            barh(popFractionPerBin,0.4,'FaceColor',[0 0.7 0.7])
            hold on
            
        else
            subplot(1,length(window),w)
            barh(popFractionPerBin,0.4,'FaceColor',[1 0.6 0]) 
            %barh(popFractionPerBin,0.4,'FaceColor',[0.5 0 0.9])
            hold on
            
        end 
        legend('fluc','low','ave')
        
        if c == 1
            figure(2)
            subplot(1,length(window),w)
            barh(popFractionPerBin,0.4)
            hold on
            title(window(w))
            axis([0,.35,0,50])
            xlabel('pdf')
            ylabel('doubling rate bin')
            legend('fluc')
        end
        
        if c == 2
            figure(3)
            subplot(1,length(window),w)
            barh(popFractionPerBin,0.4,'FaceColor',[0 0.7 0.7])
            hold on
            title(window(w))
            axis([0,.35,0,50])
            xlabel('pdf')
            ylabel('doubling rate bin')
            legend('low')
        end
        
        if c == 3
            figure(4)
            subplot(1,length(window),w)
            barh(popFractionPerBin,0.4,'FaceColor',[1 0.6 0])
            hold on
            title(window(w))
            axis([0,.35,0,50])
            xlabel('pdf')
            ylabel('doubling rate bin')
            legend('ave')
        end
        
        % 9.  repeat for all conditions in window
    end
    
    % 10.  repeat for all conditions in window 
end

%save('2017-09-26-isSubpop-muStats.mat','countMu','meanMu','semMu','stdMu')

%%
c = 1;
    
    
    % 2. isolate condition of interest
    for condition = 1:2:3
        
        interestingCondition = find(currentWindow(:,28) == condition);
        currentConditionData = currentWindow(interestingCondition,:);
        
        
        % 3. collect mu data and track IDs
        currentMus = currentConditionData(:,18);
        currentTracks = currentConditionData(:,1);
        
        
        % 4. bin mu data into 0.1 sized bins
        binnedMus = ceil(currentMus*10);
        
        % 4. remove negative growth rates
        binnedMus(binnedMus < 0) = NaN;
        nanFilter = find(~isnan(binnedMus)); %find values that are NOT nan
        binnedMus = binnedMus(nanFilter);
        currentTracks = currentTracks(nanFilter);
        
        
        % 4. accumulate track IDs according to mu bins
        binnedMus = binnedMus+1;
        tracksBinnedByMu = accumarray(binnedMus,currentTracks,[],@(x) {x});
        
        
        % 5. count unique tracks per bin
        for t = 1:length(tracksBinnedByMu)
            uniqueTracks = unique(tracksBinnedByMu{t});
            trackCountsPerBin(t,1) = length(uniqueTracks);
            clear uniqueTracks;
        end
        
        
        % 6. normalize by total counts for population fraction
        totalCounts = sum(trackCountsPerBin)
        popFractionPerBin = trackCountsPerBin/totalCounts;
        
        % 7. plot bars of population fraction
        if condition > c
            subplot(1,4,i)
            barh(popFractionPerBin,0.25,'FaceColor',[0 0.7 0.7])
            title(window(i))
            axis([0,.5,0,25])

        else
            subplot(1,4,i)
            barh(popFractionPerBin,0.5)
            title(window(i))
            axis([0,.5,0,25])
            
            hold on
        end
        xlabel('pdf')
        ylabel('doubling rate bin')
        legend('KS 1 uM','KS 100 uM','KS 10 mM')
        
    end
end