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


% last update: June 27, 2017
%              working to combine analysis for length, area, and volume

%%
clear
%load('fluorctl_2017-01-24-autoTrimmed.mat');
load('letstry-2017-06-12-autoTrimmed-scrambledOrder-editedJumpTrack.mat');
%clear D2 D3 D4 D5;
%D7=D6;
D7 = Scram6;
clear Scram2 Scram3 Scram4 Scram5;

%%

clear m n Ltrack Ttotal dT t_hr Ltime Ldiff L_Fit w Fit pFit Wdiff;
clear init_window Ttrack Ttrack_row Test Fenster Fenster_track Screen;
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
        
        
        % Original width data (microns)
        widthTrack = D7{n}(m).MinAx;
        %widthDiffs = diff(widthTrack);
        
        % Approximate volume
        vcTrack = pi * lengthTrack .* (widthTrack/2).^2;
        veTrack = 4/3 * pi * lengthTrack/2 .* (widthTrack/2).^2;
        vol_smallCylinder = (pi * (widthTrack/2).^2 .* (lengthTrack - widthTrack) );
        vol_sphere = 4/3 * pi * (widthTrack/2).^3;
        vaTrack = vol_smallCylinder + vol_sphere;
        
        % Time data (hours)
        timeTrack = T{n}/3600;
     
        % Set-up windows
        pointsInWindow = 13;                                               % sets number of frames in one window
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

                
%       A.      If dip is found after window point 7,
%               trim until last data point before dip
                if dipFinder >= 7 
                    trimmedWindow = currentWindow(1:dipFinder);
   
                    
                    % LENGTH
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
                    
                    
                    % VOLUME
                    % covert volume to log scale for linear fit
                    logVC = log(vcTrack(trimmedWindow));
                    logVE = log(veTrack(trimmedWindow));
                    logVA = log(vaTrack(trimmedWindow));
                    
                    fitVC = polyfit(trimmedTime,logVC,1);
                    fitVE = polyfit(trimmedTime,logVE,1);
                    fitVA = polyfit(trimmedTime,logVA,1);
                    
                    pFit_VC(w,:) = fitVC;
                    pFit_VE(w,:) = fitVE;
                    pFit_VA(w,:) = fitVA;
                    
                    log_FitVC = polyval(fitVC,trimmedTime);
                    log_FitVE = polyval(fitVE,trimmedTime);
                    log_FitVE = polyval(fitVA,trimmedTime);
                    
                 
                    % return to linear scale and generate exponential fit
                    slopeVC = fitVC(1);
                    interceptVC = fitVC(2);
                    vc_Fit(w,:) = exp(interceptVC)*exp(hr*slopeVC);
                    
                    slopeVE = fitVE(1);
                    interceptVE = fitVE(2);
                    ve_Fit(w,:) = exp(interceptVE)*exp(hr*slopeVE);
                    
                    slopeVA = fitVA(1);
                    interceptVA = fitVA(2);
                    va_Fit(w,:) = exp(interceptVA)*exp(hr*slopeVA);
                    
                    clear Fit Slope Intercept hr Fenster;
                    clear fitVC fitVA fitVE slopeVC slopeVE slopeVA interceptVC interceptVE interceptVA;
                end

                
%       B.      If dip is found before point 6,
%               trim points early points until dip is removed.
                if dipFinder <= 6
                     trimmedWindow = currentWindow(dipFinder+1:pointsInWindow);

                    % LENGTH 
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
                    
                    
                    % VOLUME
                    % covert volume to log scale for linear fit
                    logVC = log(vcTrack(trimmedWindow));
                    logVE = log(veTrack(trimmedWindow));
                    logVA = log(vaTrack(trimmedWindow));
                    
                    fitVC = polyfit(trimmedTime,logVC,1);
                    fitVE = polyfit(trimmedTime,logVE,1);
                    fitVA = polyfit(trimmedTime,logVA,1);
                    
                    pFit_VC(w,:) = fitVC;
                    pFit_VE(w,:) = fitVE;
                    pFit_VA(w,:) = fitVA;
                    
                    log_FitVC = polyval(fitVC,trimmedTime);
                    log_FitVE = polyval(fitVE,trimmedTime);
                    log_FitVA = polyval(fitVA,trimmedTime);
                    
                 
                    % return to linear scale and generate exponential fit
                    slopeVC = fitVC(1);
                    interceptVC = fitVC(2);
                    vc_Fit(w,:) = exp(interceptVC)*exp(hr*slopeVC);
                    
                    slopeVE = fitVE(1);
                    interceptVE = fitVE(2);
                    ve_Fit(w,:) = exp(interceptVE)*exp(hr*slopeVE);
                    
                    slopeVA = fitVA(1);
                    interceptVA = fitVA(2);
                    va_Fit(w,:) = exp(interceptVA)*exp(hr*slopeVA);
                    
                    clear Fit Slope Intercept hr Fenster;
                    clear fitVC fitVE fitVA slopeVC slopeVE slopeVA interceptVC interceptVE interceptVA;
                end
                  
                
%           when there are no length breaks in window
            else
                %disp(['Window ', num2str(w), '... smooth sailing!'])
                
                % LENGTH
                logLength = log(lengthTrack(currentWindow));                              
                hr = timeTrack(currentWindow);
                Fit = polyfit(hr,logLength,1);
                pFit(w,:) = Fit;
                
                Slope = Fit(1);
                Intercept = Fit(2);
                L_Fit(w,:) = exp(Intercept)*exp(hr*Slope);
                
                
                % VOLUME
                % covert volume to log scale for linear fit
                logVC = log(vcTrack(currentWindow));
                logVE = log(veTrack(currentWindow));
                logVA = log(vaTrack(currentWindow));
                
                fitVC = polyfit(hr,logVC,1);
                fitVE = polyfit(hr,logVE,1);
                fitVA = polyfit(hr,logVA,1);
                
                pFit_VC(w,:) = fitVC;
                pFit_VE(w,:) = fitVE;
                pFit_VA(w,:) = fitVA;
                
                %log_FitVC = polyval(fitVC,trimmedTime);
                %log_FitVE = polyval(fitVC,trimmedTime);
                
                
                % return to linear scale and generate exponential fit
                slopeVC = fitVC(1);
                interceptVC = fitVC(2);
                vc_Fit(w,:) = exp(interceptVC)*exp(hr*slopeVC);
                
                slopeVE = fitVE(1);
                interceptVE = fitVE(2);
                ve_Fit(w,:) = exp(interceptVE)*exp(hr*slopeVE);
                
                slopeVA = fitVA(1);
                interceptVA = fitVA(2);
                va_Fit(w,:) = exp(interceptVA)*exp(hr*slopeVA);
                
                
                clear fitVC fitVE fitVA slopeVC slopeVE slopeVA interceptVC interceptVE interceptVA;
                clear Fit Slope Intercept Fenster Wdiff hr Window;

            end
            
            clear Fenster Fenster_trim Wdiff hr hr_trim Slope Intercept log_Fit;
            clear log_FitVC log_FitVE log_FitVA logVC logVE logVA logLength;
             
            
        end

        % saving data
        SlidingData = struct('Parameters_L',pFit,'Smoothed_length',L_Fit,'Parameters_VC',pFit_VC,'Smoothed_vc',vc_Fit,'Parameters_VE',pFit_VE,'Smoothed_ve',ve_Fit,'Parameters_VA',pFit_VA,'Smoothed_va',va_Fit,'Windows',Wtrack);
        %M6{n}(m) = SlidingData;
        M_scrambled_edited{n}(m) = SlidingData;
        
        clear SlidingData init_window lengthDiffs lengthTrack widthTrack Wtrack;
        clear pFit L_Fit Fenster_track veTrack vcTrack vaTrack vc_Fit ve_Fit va_Fit pFit_VC pFit_VE pFit_VA;
    
        disp(['Track ', num2str(m), ' from xy ', num2str(n), ' complete!'])
    
    end
    
end


save('letstry-2017-06-12-increasedWindow-Mus-LVVV-scrambled-edited.mat', 'D_smash', 'M_scrambled_edited', 'T')
%save('monod-2016-05-25-increasedWindow-Mus-LVVV.mat', 'D6', 'M6', 'T') %'D'

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
