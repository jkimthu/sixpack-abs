%% tweakSize.m


% Goal: tweak measured lengths by some constant. can some tweak size
%       reduce sharp peaks in calculated mu?

%       assuming divisions half the true length of cells, solve for the
%       tweak constant using measured data of length drops during division


% Strategy: 

%    0. initialize trimmed track data
%    0. initialize division parameters (drop threshold)
%    1. for current movie, identify number of tracks
%           2. for each track, isolate track data
%                 3. find all division events
%                 4. for each division, solve for tweak size (tw)
%                        i. find all L1's (lengths pre-division) 
%                       ii. find all L2's (lengths at birth)
%                      iii. compute all tw's
%           5. repeat for all tracks, storing tw
%    6. find average tw in entire movie
%    7. plot average and spread of tw
%    8. subtrack average from track data


% Computing tw:

%    assuming L1 = 2 * L2...
%
%    L1 - tw / L2 - tw = 2
%    L1 - tw = 2*L2 - 2tw
%    L1 - 2*L2 = - tw


% last update: 2017 Jul 11
% 
% OK lez go!

%%
% 0. initialize 
clear
clc
experiment = '2017-06-12';

% 0. open folder for experiment of interest
newFolder = strcat('/Users/jen/Documents/StockerLab/Data/',experiment);
cd(newFolder);


% 0. initialize trimmed track data
load('letstry-2017-06-12-revisedTrimmer-xy1-xy52-noLinker.mat','D7','T');
numMovies = length(D7);


% 0. initialize division parameters
dropThreshold = -0.75;


%  1. for each movie, identify the number of tracks
n = 52;
numTracks = length(D7{n});


%  2. per track, solve for all tw and store
tw = [];

for track = 1:numTracks
    
    % 2. isolate length, frame and time data
    trackLengths = D7{n}(track).MajAx;
    trackFrames = D7{n}(track).Frame;
    trackTimes = T{n}(trackFrames)/3600;
    
    % 3. find all division events (all changes in size that are more negative than threshold)
    sizeChange = diff(trackLengths);
    dropTrack = find(sizeChange <= dropThreshold);
    
    clear sizeChange
    
    
    % 4. for each division, solve for tweak size (tw)
    %for i = 1:length(dropTrack)
    
    % i. find all L1's (lengths pre-division)
    L1 = trackLengths(dropTrack);
    
    % ii. find all L2's (lengths at birth)
    L2 = trackLengths(dropTrack+1);
    
    % iii. compute all tw's
    tw = [tw; L1 - 2*L2];
    

    % 5. repeat for all tracks, storing tw
end
%%
%  6. find average tw in entire movie
aveTweak = mean(tw);
stdTweak = std(tw);

% 7. plot average and spread of tw
figure(1)
histogram(tw)
title(aveTweak)


%% sliding fits (smoothed) for track 1
track = 1;

% 2. isolate length, frame and time data
trackLengths = D7{n}(track).MajAx;
trackFrames = D7{n}(track).Frame;
trackTimes = T{n}(trackFrames)/3600;

% 3. find all division events (all changes in size that are more negative than threshold)
sizeChange = diff(trackLengths);
dropTrack = find(sizeChange <= dropThreshold);

%  4. build an array that identifies curve #

% 0. initalize array, curveNum
curveNum = zeros(length(trackFrames),1);

% % i. identify all changes in size > threshold
% sizeChange = diff(trackLengths);
% dropTrack = find(sizeChange <= dropThreshold);

% ii. starting with zero, list curve # for each frame
currentCurve = 0;
nextCurve = 1;
if ~isempty(dropTrack)
    for i = 1:length(curveNum)
        if nextCurve <= length(dropTrack)
            curveNum(i) = currentCurve;
            
            if i == dropTrack(nextCurve)
                currentCurve = currentCurve+1;
                nextCurve = nextCurve+1;
            end
            
        else
            curveNum(i) = currentCurve;
        end
    end
end
clear currentCurve nextCurve i sizeChange;

%%
% 5. initialize windows for current track
    
windowSize = 5;

% i. define row numbers for first window
firstWindow = 1:windowSize; 

% ii. calculate number of windows in current track
numWindows = length(trackLengths) - windowSize +1;
 
% iii. initialize current window to begin calculations
currentWindow = firstWindow; % will slide by +1 for each iternation


%  6. per window
for w = 1:numWindows
    
    % 6. isolate current window's time, length, and curve #
    wLength = trackLengths(currentWindow) + aveTweak;
    wCurves = curveNum(currentWindow);
    wTime = trackTimes(currentWindow);
    
    % 7. if curve # changes, adjust window Lengths by one of two means:
    isDrop = diff(wCurves);
    if sum(isDrop) ~= 0
        
        % i. find point of drop
        dropPoint = find(isDrop ~= 0);
        
        %         % if drop is in latter half of curve...
        %         if dropPoint > 2
        
        % ii. double length values after change
        multiplier = NaN(windowSize,1);
        minCurve = min(wCurves);
        maxCurve = max(wCurves);
        multiplier(wCurves == minCurve) = 1;
        multiplier(wCurves == maxCurve) = 2;
        
        wLength_adjusted = wLength.*multiplier;
        
        %         else
        %
        %             % if drop is in latter half of curve...
        %             % ii. double length values after change
        %             multiplier = NaN(windowSize,1);
        %             minCurve = min(wCurves);
        %             maxCurve = max(wCurves);
        %             multiplier(wCurves == minCurve) = 0.5;
        %             multiplier(wCurves == maxCurve) = 1;
        %
        %             wLength_adjusted = wLength.*multiplier;
        %
        %         end
        
        % iii. use effective length to calculate mu (see below for comments)
        ln_length = log(wLength_adjusted);
        fitLine = polyfit(wTime,ln_length,1);
        mu = fitLine(1)/log(2);         % divide by log(2), as eqn raises 2 by mu*t
        
        
    %  8. if no change in curve #, use effective length to calculate mu
    else
        % i. ln(effective length) vs time
        ln_length = log(wLength);
        
        % ii. fit linear slope to ln(eL) vs time
        fitLine = polyfit(wTime,ln_length,1); % gives: fitLine = [0.6931, -0.0000]
                                        % where:   slope = fitLine(1)
                                        %          y-int = fitLine(2)
        
        % iii. mu = slope / ln(2)
        mu = fitLine(1)/log(2);         % divide by log(2), as eqn raises 2 by mu*t
                                 
    end
    
    % 9. save mu and y-intercept
    slidingData(w,1) = mu;
    slidingData(w,2) = fitLine(2); % log(initial length), y-int
    
    % 10. repeat for all windows
    currentWindow = currentWindow + 1;
    
    clear wLength wCurves mu fitLine ln_length dropPoint isDrop maxCurve minCurve;
    clear wLength_adjusted multiplier experiment newFolder wTime;
end

%                11. repeat for all tracks
%       12. repeat for all movies

%% checks
% plot mu over time (like length) 
                      
figure(2)
plot(trackFrames, trackLengths,'o')
hold on
plot(trackFrames, trackLengths,'r')
grid on
xlim([0 202])
title(track);


hold on
plot(trackFrames(3:end-2),slidingData(:,1)*2,'Color',[1 0.5 0],'Marker','o'); 
hold on
plot(trackFrames(3:end-2),slidingData(:,1)*2,'Color',[0.5 0 0.5]); 
hold on
plot(trackFrames(3:end-2),slidingData(:,1),'Color','k'); 
xlim([0 202])
title(track);
%             


