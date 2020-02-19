%% slidingFits_helpfulNeighbors (using volume, as approx. by a capped cylinder)

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
%                        3. calculate volume data
%                        4. build an array with length(track) that identifies curve #
%                                i. identify all changes in size > proportional threshold (drop greater than 30% previous size)
%                               ii. starting with zero, list curve # for each frame
%                        5. initialize windows for current track
%                               i. define row numbers for first window
%                              ii. calculate number of windows in current track
%                             iii. initialize current window to begin calculations
%                        6. per window, isolate effective time, length, volume, and curve #
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


% last update: 2018 September 18
% commit: 2018-09-17 with constant width threshold



% OK lez go!

                         
%%
% 0. initialize 
clear
clc
experiment = '2018-09-17';

% 0. open folder for experiment of interest
newFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',experiment);%,'  (t300)');
cd(newFolder);


% 0. initialize trimmed track data

%load(strcat('lb-fluc-',experiment,'-width1p7-jiggle-0p5.mat'),'D','D5','T','rejectD');
load(strcat('lb-fluc-',experiment,'-width1p7-jiggle-0p5.mat'),'D','D5','T','rejectD');
numMovies = length(D5);


% 0. initialize window parameters
windowSize = 5;


% 0. initialize division parameters
dropFrac = -0.3;


%%
%  1. for each movie, identify the number of tracks


for n = 1:length(D5)
    
    numTracks = length(D5{n});
 
    %  2. per track, isolate length, width, frame and time data
    for track = 1:numTracks
        trackLengths = D5{n}(track).MajAx;
        trackWidths = D5{n}(track).MinAx;
        trackFrames = D5{n}(track).Frame;
        trackTimes = T{n}(trackFrames)/3600; % in sec converted to hr
        
        
        %  3. calculate volume, as approximated by a capped cylinder 
        vol_smallCylinder = (pi * (trackWidths/2).^2 .* (trackLengths - trackWidths) );
        vol_sphere = 4/3 * pi * (trackWidths/2).^3;
        trackVolumes = vol_smallCylinder + vol_sphere;  % approx. volume as cylinder with spherical caps
        

        %  4. build an array that identifies curve #
                %  0. initialize array (curve data)
                %  i. identify all changes in size > threshold (drop greater than 30% previous size)
                % ii. starting with zero, list curve # for each frame
        
        % 0. initalize array, curveNum
        curveNum = zeros(length(trackFrames),1);
        
        % i. identify all changes in size are MORE NEGATIVE than -0.3 (proportional threshold)
        sizeChange = diff(trackVolumes);
        growthFrac = sizeChange./trackVolumes(1:length(sizeChange));
        
        trackDrops = growthFrac <= dropFrac; % 0 = no drop, 1 = drop
        
        
        % ii. starting with zero, list curve # for each frame
        % from buildDM
        numberFullCurves = sum(trackDrops) - 1;                                % all curves start and end with a division, isDrop = 1                                      
        curveTrack = zeros(length(trackDrops),1);
        
        % find and number the full curves within a single track
        curveCounter = 0;                                                  
        for i = 1:length(trackDrops) 
            
            if trackDrops(i) == 0                           % 1. disregard incomplete first curve
                curveTrack(i,1) = curveCounter;             %    by starting count at 0   
            
            elseif (trackDrops(i) == 1)
                curveCounter = curveCounter + 1;            % 2. how to disregard final incomplete segment? 
                                                            %    stop when curveCount exceeds number of fullCurves
                if curveCounter <= numberFullCurves         
                    curveTrack(i,1) = curveCounter;
                else                                        
                    break                                   % all incomplete curves are left as 0               
                end
                
            end
            
        end
        trackCurves = [0; curveTrack]; % add a zero, to match V track, not diff track
        clear curveCounter i
        
        
        
        % 5. initialize windows for current track
        % i. define row numbers for first window
        firstWindow = 1:windowSize;
        
        % ii. calculate number of windows in current track
        numWindows = length(trackLengths) - windowSize +1;
        
        % iii. initialize current window to begin calculations
        currentWindow = firstWindow; % will slide by +1 for each iternation

        
        %  6. per window
        for w = 1:numWindows
            
            % 6. isolate current window's time, length, volume, and curve #
            wLength = trackLengths(currentWindow);
            wVolume = trackVolumes(currentWindow);
            wCurves = curveNum(currentWindow);
            wTime = trackTimes(currentWindow);
            
            % 7. if curve # changes, adjust window Lengths by one of two means:
            isDrop = diff(wCurves);
            if sum(isDrop) ~= 0
                
                % i. find point of drop
                dropPoint = find(isDrop ~= 0);
              
                % ii. double length or volume values after change
                multiplier = NaN(windowSize,1);
                minCurve = min(wCurves);
                maxCurve = max(wCurves);
                multiplier(wCurves == minCurve) = 1;
                multiplier(wCurves == maxCurve) = 2;
                
                wLength_adjusted = wLength.*multiplier;
                wVolume_adjusted = wVolume.*multiplier;
                
                % iii. use effective length or volume to calculate mu (see below for comments)
                ln_length = log(wLength_adjusted);
                fitLine_l = polyfit(wTime,ln_length,1);
                mu_l = fitLine_l(1)/log(2);
                
                ln_volume = log(wVolume_adjusted);
                fitLine_va = polyfit(wTime,ln_volume,1);
                mu_va = fitLine_va(1)/log(2);         % divide by log(2), as eqn raises 2 by mu*t
                
                
                %  8. if no change in curve #, use effective length to calculate mu
            else
                % i. ln(effective length or volume) vs time
                ln_length = log(wLength);
                ln_volume = log(wVolume);
                
                % ii. fit linear slope to ln(eL or eV) vs time
                fitLine_l = polyfit(wTime,ln_length,1);
                fitLine_va = polyfit(wTime,ln_volume,1);
                % where:   slope = fitLine(1)
                %          y-int = fitLine(2)
                
                % iii. mu = slope / ln(2)
                mu_l = fitLine_l(1)/log(2);  
                mu_va = fitLine_va(1)/log(2);      % divide by log(2), as eqn raises 2 by mu*t
                
            end
            
            % 9. save mus and y-intercept
            slidingData_l(w,1) = mu_l;
            slidingData_l(w,2) = fitLine_l(2); % log(initial length), y-int
            
            slidingData_va(w,1) = mu_va;
            slidingData_va(w,2) = fitLine_va(2); % log(initial volume), y-int
            
            % 10. repeat for all windows
            currentWindow = currentWindow + 1;
            
            clear wVolume wCurves mu fitLine ln_volume dropPoint isDrop maxCurve minCurve;
            clear wVolume_adjusted multiplier newFolder wTime;
        end
        
        % 11. save data and repeat for all tracks
        trackData_l = struct('mu',slidingData_l(:,1),'yInt',slidingData_l(:,2));
        M{n}(track) = trackData_l;
        
        trackData_va = struct('mu_va',slidingData_va(:,1),'yInt',slidingData_va(:,2));
        M_va{n}(track) = trackData_va;
        
        clear slidingData_l slidingData_va trackData_l trackData_va;
        
        disp(['Track ', num2str(track), ' from xy ', num2str(n), ' complete!'])
        
    end
    
    %       12. repeat for all movies
end


save(strcat('lb-fluc-',experiment,'-window5-width1p7-jiggle-0p5.mat'), 'D','D5','M', 'M_va', 'T','rejectD') %'D'


%% checks
% plot mu over time (like length) 
clear
load('mopsvnc-2017-05-26-revisedTrimmer-jiggle0p4.mat','D5','D','T','rejectD');


%%
% re-create volume track with y-int

% initialize calculated mu and y-intercept data
slidingData_va = M_va{n}(track);
mus = slidingData_va.mu_va;
yInts = slidingData_va.yInt;

% trim original data to length of mu data
volumes = trackVolumes(3:length(mus)+2);
times = trackTimes(3:length(mus)+2);

% reconstruct volume from calculations
ln_calc_volume = yInts + mus.*times *log(2);
fit_volume = exp(ln_calc_volume);

figure(2)
plot(times, volumes,'o')
hold on
plot(times, fit_volume,'ro')
hold on
plot(times, mus, 'k-')
legend('true volume','fit volume','mu')
xlim([0 10])
title(track);

%% adding additional xys to D, etc
% used to compile data from 2018-01-29 xy1 and xy2 into previously analyzed
% matrixes for xy3-40

clear
clc
load('lb-fluc-2018-01-29-xy1and2-window5.mat')

D_12 = D;
D5_12 = D5;
M_12 = M;
M_va_12 = M_va;
rejectD_12 = rejectD;
T_12 = T;

clear D D5 M M_va rejectD T

load('lb-fluc-2018-01-29-window5-width1p4-1p7-jiggle-0p5.mat')

D(1:2) = D_12;
D5(1:2) = D5_12;
M(1:2) = M_12;
M_va(1:2) = M_va_12;
rejectD(:,1) = rejectD_12(:,1);
rejectD(:,2) = rejectD_12(:,2);
T(1:2) = T_12(1:2);

%%
save('lb-fluc-2018-01-29-window5-width1p4-1p7-jiggle-0p5.mat', 'D','D5','M', 'M_va', 'T','rejectD')






