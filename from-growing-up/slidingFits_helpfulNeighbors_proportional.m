%% slidingFits_helpfulNeighbors

% Goal: a method of calculating doubling rate (mu), in which we use the
% points in a neighboring curve to add more points for fitting.
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


% last update: jen, 2017 Jul 12

% OK lez go!

%%
% testing mu calculations
doublingRate = 1.0000;
t = [0,1,2,3,4];
data = 1 * 2.^(doublingRate*t);  % gives: data = [1,2,4,8,16]

ln_length = log(data);  % gives: ln_data = [0, 0.6931, 1.3863, 2.0794, 2.7726]

plot(t,data) % exponential
plot(t,ln_length) % linear

fitLine = polyfit(t,ln_length,1); % gives: fitLine = [0.6931, -0.0000]
                                % where: slope = fitLine(1)
                                %        y-int = fitLine(2)
                                
mu = fitLine(1)/log(2);         % gives: mu = 1
                                % woot! we found the doublingRate from the data!
                         
%%
% 0. initialize 
clear
clc
experiment = '2017-01-16';

% 0. open folder for experiment of interest
newFolder = strcat('/Users/jen/Documents/StockerLab/Data/',experiment,'  (t300)');
cd(newFolder);


% 0. initialize trimmed track data

load('t300_2017-01-16-revisedTrimmer-jiggle0p4.mat','D5','T');
numMovies = length(D5);


% 0. initialize window parameters
windowSize = 5;


% 0. initialize division parameters
dropFrac = -0.3;

%%
%  1. for each movie, identify the number of tracks


for n = 1:10%length(D5)
    
    numTracks = length(D5{n});
 
    %  2. per track, isolate length, width, frame and time data
    for track = 1: numTracks
        trackLengths = D5{n}(track).MajAx;
        trackFrames = D5{n}(track).Frame;
        trackTimes = T{n}(trackFrames)/3600;

        %  4. build an array that identifies curve #
        
        % 0. initalize array, curveNum
        curveNum = zeros(length(trackFrames),1);
        
        % i. identify all changes in size > proportional threshold
        sizeChange = diff(trackLengths);
        growthFrac = sizeChange./trackLengths(1:length(sizeChange));
        
        dropTrack = find(growthFrac <= dropFrac);
    
        
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
        
        % 11. save data and repeat for all tracks
        
        trackData = struct('mu',slidingData(:,1),'yInt',slidingData(:,2));
        M_prop{n}(track) = trackData;
        
        clear slidingData trackData;
        
        disp(['Track ', num2str(track), ' from xy ', num2str(n), ' complete!'])
        
    end
    
    %       12. repeat for all movies
end

%%
save('t300_2017-01-16-neighbors-prop-jiggle0p4.mat', 'D5', 'M_prop', 'T') %'D'


%% checks
% plot mu over time (like length) 

load('t300_2017-01-16-neighbors-prop-jiggle0p4.mat');


n=1;
track = 1;

trackLengths = D5{n}(track).MajAx;
trackFrames = D5{n}(track).Frame;
trackTimes = T{n}(trackFrames)/3600;
trackMus = M_prop{n}(track).mu;


figure(1)
plot(trackFrames, trackLengths,'o')
hold on
plot(trackFrames, trackLengths,'r')
grid on
xlim([0 202])
title(track);


%hold on
%plot(trackFrames(3:end-2),trackMus*4,'Color',[1 0.5 0],'Marker','o'); 
%hold on
%plot(trackFrames(3:end-2),trackMus*4,'Color',[0.5 0 0.5]); 
hold on
plot(trackFrames(3:end-2),trackMus,'bo','MarkerSize',10); 
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



