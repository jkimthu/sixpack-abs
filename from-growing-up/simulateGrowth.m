%% simulateGrowth.m

% Goal: plot a length trajectory over time, of known mu.
%       test how adding different proprotions of size to the "true" length
%       can change how we calculate mu.


% Strategy: 
%       0. initialize mu, size at birth, and time
%       0. initialize array of different added proportions
%       1. plot defined trajectories
%       2. for each trajectory,
%               3. find smoothed mu
%               4. find non-smoothed mu
%       5. determine average mu for all cases
%       6. plot 
%



% last update: 2017 Jul 11
% 
% OK lez go!

%%
% 0. initialize mu, size at birth, and time
clear

mu = .5; % 1/hr
birthSize = 1.5; % um
doublingSize = 2 * birthSize;

time = linspace(0,300,101)/60; % every 3 min for 10 hours (in hours)
timestep = time(2) - time(1);
tolerance = 0.01;

% 0. initialize array of different added proportions
addedProps = -25:25:25;

% 1. plot defined trajectories
lengthTrack = NaN(length(time),1);


% i. Lt = birthSize * 2^(mu*time)
%
Li = birthSize;
for i = 1:length(time)
    
    Lt = Li * (2^(mu*timestep));
    
    if abs(Lt - doublingSize) < tolerance
        
        % a drop! record.
        lengthTrack(i) = Lt;
        
        % re-define current size as start size for next iteration
        Li = birthSize;
        clear Lt
        
    else
        
        % record current size
        lengthTrack(i) = Lt;
        
        % re-define current size as start size for next iteration
        Li = Lt;
        clear Lt
        
    end
    
end

figure(1)
plot(time,lengthTrack)
%axis([0 5.1 2 9])

clear i Li timestep
    
    %%
    
for a = 1:length(addedProps)
    % 2. for each trajectory, initialize data
    
    % 0. initialize size gain (error)
    sizeGain = addedProps(a) / 100;
    errLengthTrack = lengthTrack + (birthSize*sizeGain);
    
    figure()
    plot(time,lengthTrack)
    hold on
    plot(time,errLengthTrack,'r')
    axis([0 5.1 -1 5])
    
    
    windowSize = 5;
    dropThreshold = -0.75;
    
    track = a;
    trackLengths = errLengthTrack;
    %trackWidths = D7{n}(track).MinAx;
    %trackFrames = D7{n}(track).Frame;
    trackTimes = time';
    
    
    % 3. find smoothed mu
    
    % 0. initalize array, curveNum
    curveNum = zeros(length(trackLengths),1);
    
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
    
    % determine ave and std of mu
    aveMu_smoothed(a,1) = mean(slidingData(:,1));
    stdMu_smoothed(a,1) = std(slidingData(:,1));
    
    %%
    figure()
    plot(trackTimes, trackLengths,'o')
    hold on
    plot(trackTimes, trackLengths,'r')
    grid on
    %xlim([0 202])
    title(track);
    
    
    hold on
    plot(trackTimes(3:end-2),slidingData(:,1),'Color','k');
    text(4.5, aveMu_smoothed(a) + 0.2, num2str(aveMu_smoothed(a)),'Color','k','FontSize',12);
    %xlim([0 202])
    title(track);
    axis([0 5.1 -1 5])
    
    
    
    %%
    %               4. find non-smoothed mu
    
    windowSize =2;
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
        slidingData_nonSmoothed(w,1) = mu;
        %slidingData(w,2) = fitLine(2); % log(initial length), y-int
        
        % 10. repeat for all windows
        currentWindow = currentWindow + 1;
        
        clear wLength wCurves mu fitLine ln_length dropPoint isDrop maxCurve minCurve;
        clear wLength_adjusted multiplier experiment newFolder wTime;
    end
    
    % 5. determine average mu for all cases
    aveMu_unsmoothed(a,1) = mean(slidingData_nonSmoothed);
    stdMu_unsmoothed(a,1) = std(slidingData_nonSmoothed);
    %%
    figure()
    plot(trackTimes, trackLengths,'o')
    hold on
    plot(trackTimes, trackLengths,'r')
    grid on
    title(track);
    
    hold on
    plot(trackTimes(2:end),slidingData_nonSmoothed,'Color',[0.5 0 0.5]);
    text(4.5, aveMu_unsmoothed(a) + 0.2, num2str(aveMu_unsmoothed(a)),'Color','k','FontSize',12);
    axis([0 5.1 -1 5])
    
    title(track);
    

    
end

% 6. plot





