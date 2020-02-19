%% Determination of growth rate by fitting an exponential to cell size trajectories
%


%  After quality control step is complete, D6 is a row of cells -- each
%  representing a different experimental movie / xy-position.
%           ex.  D6{n} = movie n from dataset D6
%
%  Each cell containss a column of structures, each the data for one track
%  from that movie.
%           ex.  D6{n}(m} = track m from movie n
%
%  Each structure contains a variety of metrics measured for that track,
%  including: X, Y, Area, Frame, and MajAx
%           ex.  D6{n}(m).MajAx = column of cell length values (double) for
%                                 each frame (row #) of track m
%


%  GOAL: measure the instantaneous growth rate for each track
%


%  APPROACH: fit sections of each track to the following exponential
%
%                   S(t) = s * 2^(mu*t)
%
%            where, S = final cell size
%                   s = initial cell size
%                   mu = growth rate over time t
%                   t = time window
%
%       Note: frame is NOT equal to t, but multiple of it. Instead, t
%       depends on how many points are used to calculate fit (measure mu)
%
%       Sections in this script are defined by a sliding window, which contains an odd
%       number of points. The calculated mu is then assigned to the
%       centermost timepoint.
%
%

%  OK! HERE WE GO!


% last update: Mar 29rd, 2017
%              working to combine analysis for length, area, and volume

%%

load('fluorctl_2017-01-24-autoTrimmed.mat');
%load('t30_2017-01-09-autoTrimmed.mat');
clear D2 D3 D4 D5;

%%

clear m n Ltrack Ttotal dT t_hr Ltime Ldiff L_Fit w Fit pFit Wdiff;
clear total_windows init_window Ttrack Ttrack_row Test Fenster Fenster_track Screen;
clear Fenster_trim hr_trim log_Fit SlidingData;


% 1. Initialize:
%
%       i.      Load length trajectory for track m
%       ii.     Determine length of time window
%       iii.    Establish whether window contains a division event
%       iv.     If it does, then remove dip from fit and double lengths of subsequent points
%


for n = 1:length(D7)
    
    for m = 1:length(D7{n})
      
        % Original length data (microns)
        lengthTrack = D7{n}(m).MajAx;                                      % loads current length trajectory
        lengthDiffs = diff(lengthTrack);                                   % used to find sharp drops
        
        % Time data (hours)
        timeTrack = T{n}/3600;
     
        % Set-up windows
        pointsInWindow = 13;                                                % sets number of frames in one window
        firstWindow = linspace(1,pointsInWindow,pointsInWindow);           % defines frame numbers for first window
        numWindows = length(lengthDiffs) - (pointsInWindow-1);             % total windows in track
        
        
        
        
        % Fitting directions for special cases, where a window contains a break in growth (division event)
        for w = 1:numWindows
            
            % Determine frames of analysis
            currentWindow = firstWindow + (w-1);                           % defines vector of frame numbers
            Wtrack(w,:) = currentWindow;
            Wdiff = lengthDiffs(currentWindow(1:pointsInWindow-1));        % incremental length differences in current window
            
            % Working around sharp dips in length
            dipFinder = find(Wdiff < -.75);  % returns 1 if all diffs are above threshold, 0 if a dip was found
            
            % When a window has a dip, remove dip from analysis
            % Method used depends on where dip is located:
            if isempty(dipFinder) == 0

                
%       A.      If dip is found between window points 4 and 5,
%               trim last data point to remove dip, only use first 4 for fit
                if dipFinder >= 7 %length(Wdiff)
                    trimmedWindow = currentWindow(1:dipFinder);
                    %Wtrack(w,:) = [trimmedWindow 0];
                    %disp(['Window ', num2str(w), '... Red!'])
                    
                    % covert length to log scale for linear fit
                    logLength = log(lengthTrack(trimmedWindow));
                    trimmedTime = timeTrack(trimmedWindow);
                    Fit = polyfit(trimmedTime,logLength,1);
                    pFit(w,:) = Fit;
                    log_Fit = polyval(Fit,trimmedTime);
                    
                    % return to linear scale and generate exponential fit
                    Slope = Fit(1);
                    Intercept = Fit(2);
                    hr = timeTrack(currentWindow);
                    L_Fit(w,:) = exp(Intercept)*exp(hr*Slope);
                    
                    clear Fit Slope Intercept log_L hr Fenster;
                end

                
%       B.    
                
%       C.      If dip is found between window points 4 and 5, or 3 and 4, or 2 and 3,
%               no data recorded.
%                 if dipFinder <= 6 
%                     pFit(w,:) = [0 0];
%                     L_Fit(w,:) = zeros(1,pointsInWindow);  
%                 end

                if dipFinder <= 6
                     trimmedWindow = currentWindow(dipFinder+1:pointsInWindow);
                    %Wtrack(w,:) = [trimmedWindow 0];
                    %disp(['Window ', num2str(w), '... Red!'])
                    
                    % covert length to log scale for linear fit
                    logLength = log(lengthTrack(trimmedWindow));
                    trimmedTime = timeTrack(trimmedWindow);
                    Fit = polyfit(trimmedTime,logLength,1);
                    pFit(w,:) = Fit;
                    log_Fit = polyval(Fit,trimmedTime);
                    
                    % return to linear scale and generate exponential fit
                    Slope = Fit(1);
                    Intercept = Fit(2);
                    hr = timeTrack(currentWindow);
                    L_Fit(w,:) = exp(Intercept)*exp(hr*Slope);
                end
                  
                
%           when there are no length breaks in window
            else
                %disp(['Window ', num2str(w), '... smooth sailing!'])
                
                logLength = log(lengthTrack(currentWindow));                              
                hr = timeTrack(currentWindow);
                Fit = polyfit(hr,logLength,1);                                
                pFit(w,:) = Fit;
                
                Slope = Fit(1);
                Intercept = Fit(2);
                L_Fit(w,:) = exp(Intercept)*exp(hr*Slope);
                
                clear Fit Slope Intercept Screen Fenster Wdiff;
                clear log_L hr Window;
                
            end
            
            clear Fenster Fenster_trim Wdiff Screen log_L hr hr_trim Slope Intercept log_Fit;
            
        end

        % saving data
        SlidingData = struct('Parameters',pFit,'Smoothed_length',L_Fit,'Windows',Wtrack);
        M6{n}(m) = SlidingData;
        
        clear SlidingData total_windows init_window Ltrack Ldiff total_windows;
        clear pFit L_Fit Fenster_track;
    
        disp(['Track ', num2str(m), ' from xy ', num2str(n), ' complete!'])
    
    end
    
end


save('fluorctl_2017-01-24-increasedWindow-Mus-length.mat', 'D6', 'M6', 'T')
%save('t30_2017-01-09-increasedWindow-Mus-length.mat', 'D6', 'M6', 'T') %'D'
clear Fenster_track L_Fit Ltime pFit t_hr;

%%

% Visualizing the fit!
for w = 1:numWindows
    fitTrack(w,:) = exp( pFit(w,1)*timeTrack(Wtrack(w,:)) + pFit(w,2) );
end
clear w

% totalTpt = length(timeTrack);
% calcTpt = length(fitTrack);

figure(1)
plot(timeTrack,lengthTrack)
hold on
plot(timeTrack(7:492), fitTrack(:,7),'r.');
hold on
plot(timeTrack(7:492), L_Fit(:,7), 'g.');
grid on;
axis([0,9,1.5,5])
xlabel('Time (hours)')
ylabel('Microns')
legend('Original Length', 'Smoothed Length');

% Vis method 2
for w = 1:numWindows
    polyTrack(w,:) = polyval( pFit(w,:),timeTrack(Wtrack(w,:)) );
end
clear w

figure(2)
plot(timeTrack,lengthTrack)
hold on
plot(timeTrack(7:492),exp(polyTrack(:,7)),'r.')
grid on
axis([0,9,1.5,5.5])
xlabel('Time (hours)')
ylabel('Microns')
legend('Original Length', 'Smoothed Length');
