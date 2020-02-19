%% slidingFits_helpfulNeighbors_unsmoothed

% Goals: 
%       1. a method of calculating doubling rate (mu), in which we use the
%          points in a neighboring curve to add more points for fitting.
%
%       2. direct comparison with slidingFits_helpfulNeighbors.m with
%          window size = 5, which smooths across a few points and produces
%          oscillations in mu that coincide heavily with division cycles
%
% 
% Strategy:
%
%       0. initialize trimmed track data
%       0. initialize window parameters (number of frames)
%       0. initialize division parameters (drop threshold)
%       1. for each movie, identify the number of tracks
%                2. per track, isolate length, width, frame and time data
%                3. calculate volume data
%                        4. build an array with length(track) that identifies curve #
%                                i. identify all changes in size > threshold (-0.75 um)
%                               ii. starting with zero, list curve # for each frame
%                        5. initialize windows for current track
%                        6. per window, isolate effective length and curve #
%                        7. if curve # changes,
%                                i. double length values after change
%                               ii. use effective length to calculate mu
%                        8. if no change in curve #, use effective length to calculate mu
%                                i. ln(effective length) vs time
%                               ii. fit linear slope to ln(eL) vs time
%                              iii. mu = slope / ln(2)
%                        9. save mu and y-intercept
%                       10. repeat for all windows
%                11. repeat for all tracks
%       12. repeat for all movies


% last update: jen, 2017 Jul 11

% OK lez go!

%%
% mu calculations from two points

% point one
% L1 = A * 2^(mu*T1)

% point two
% L2 = A * 2^(mu*T2)

% assuming no change in mu between T1 and T2
% L2/L1 = 2^(mu * (T2-T1))
% or,
% L2/L1 = 2^(mu * deltaT)

% ln(L2/L1) = mu * deltaT * ln(2)

% mu = ln(L2/L1) / (deltaT * ln(2))


ln_ratioL = log(wLength(2)/wLength(1));
deltaT = wTime(2) - wTime(1);
                 
mu = ln_ratioL/(log(2)*deltaT);

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


% 0. initialize window parameters
windowSize = 2;


% 0. initialize division parameters
dropThreshold = -0.75;


%  1. for each movie, identify the number of tracks
n = 52;
numTracks = length(D7{n});



%  2. per track, isolate length, width, frame and time data
track = 1;
trackLengths = D7{n}(track).MajAx;
trackWidths = D7{n}(track).MinAx;
trackFrames = D7{n}(track).Frame;
trackTimes = T{n}(trackFrames)/3600;
%trackID = D7{n}(track).TrackID;


% 3. calculate volume data

% i. approximate volume as a cylinder = pi * r_squared * height
r = trackWidths/2;
r_squared = r.^2; 
v_cylinder = pi * trackLengths .* r_squared;               


% ii. approximate volume as an ellipsoid = 4/3 * pi * a * b * c
% for us, b = c = radius of circular cross-section
a = trackLengths/2;
v_ellipse = 4/3 * pi * a .* r_squared;         


% iii. approximate volume as a cylinder with half-sphere caps
% cylinder = pi * r_squared * height
shortenedHeight = trackLengths - trackWidths;
vol_smallCylinder = pi * r_squared .* shortenedHeight;

% sphere = 4/3 * pi * r_cubed
r_cubed = r.^3;
vol_sphere = 4/3 * pi * r_cubed;

% v = cylinder + sphere
v_anupam = vol_smallCylinder + vol_sphere;

clear r r_squared r_cubed a shortenedHeight


%  4. build an array that identifies curve #

% 0. initalize array, curveNum
curveNum = zeros(length(trackFrames),1);

% i. identify all changes in size > threshold
sizeChange = diff(trackLengths);
dropTrack = find(sizeChange <= dropThreshold);

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

% i. define row numbers for first window
firstWindow = 1:windowSize; 

% ii. calculate number of windows in current track
numWindows = length(trackLengths) - windowSize +1;
 
% iii. initialize current window to begin calculations
currentWindow = firstWindow; % will slide by +1 for each iternation


%  6. per window
for w = 1:numWindows
    
    % 6. isolate current window's time, length, and curve #
    wLength = trackLengths(currentWindow);
    %wVolume = v_anupam(currentWindow);
    wCurves = curveNum(currentWindow);
    wTime = trackTimes(currentWindow);
    
    % 7. if curve # changes, adjust window Lengths by one of two means:
    isDrop = diff(wCurves);
    if sum(isDrop) ~= 0
        
        % i. find point of drop
        dropPoint = find(isDrop ~= 0);
        
        
        % ii. double length values after change
        multiplier = NaN(windowSize,1);
        minCurve = min(wCurves);
        maxCurve = max(wCurves);
        multiplier(wCurves == minCurve) = 1;
        multiplier(wCurves == maxCurve) = 2;
        
        wLength_adjusted = wLength.*multiplier;
        
        % calculate mu between two timepoints
        ln_ratioL = log(wLength_adjusted(2)/wLength_adjusted(1));
        deltaT = wTime(2) - wTime(1);
        
        mu = ln_ratioL/(log(2)*deltaT);

        
        
    %  8. if no change in curve #, use effective length to calculate mu
    else
        % calculate mu between two timepoints
        ln_ratioL = log(wLength(2)/wLength(1));
        deltaT = wTime(2) - wTime(1);
        
        mu = ln_ratioL/(log(2)*deltaT);
        
    end
    
    % 9. save mu and y-intercept
    slidingData(w,1) = mu;
    %slidingData(w,2) = fitLine(2); % log(initial length), y-int
    
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
plot(trackFrames(2:end),slidingData(:,1)*2,'Color',[1 0.5 0],'Marker','o'); 
hold on
plot(trackFrames(2:end),slidingData(:,1)*2,'Color',[0.5 0 0.5]); 
hold on
plot(trackFrames(2:end),slidingData(:,1),'Color','k'); 
xlim([0 202])
title(track);


%%
% re-create length track with y-int
firstWindow = 1:windowSize; 
currentWindow = firstWindow;

for i = 1:length(slidingData)
    
    wTime = trackTimes(i+2);
    %timestep = max(wTime) - min(wTime);
    
    mu = slidingData(i,1);
    yInt = slidingData(i,2);
    ln_length = yInt + mu*wTime *log(2);
    
    fit_Length(i,1) = exp(ln_length);
    
    currentWindow = currentWindow +1;
end


plot(trackFrames, trackLengths,'o')
hold on
plot(trackFrames, trackLengths,'r')
hold on
plot(trackFrames(3:end-2), fit_Length, 'ko')

xlim([0 202])
title(track);

%fitTrack(w,:) = exp( pFit(w,1)*timeTrack(Wtrack(w,:)) + pFit(w,2) );

% log(exp(1)) = 1
%%
%                        
%                10. repeat for all tracks
%       11. repeat for all movies



