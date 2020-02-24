%%  matrixBuilder_untrimmed


%  Goal: Assembling growth data to facilitate plotting the bejeezy out of it!
%        This script define growth phase as a specific fraction of the growth curve
%        and tables it along with all other associated data into an awesome,
%        organized matrix
%
%  Last edit: Jen Nguyen, June 29 2017




% Envisioned data matrix:

%        row      Track#    Time     Lngth      Mu       drop?      curve#    timeSinceBirth    curveDuration    cc stage    lngthAdded    addedLngth    Width    V_cyl   V_elpse  V_anu  mu_vc  mu_ve  mu_va  Condition o
%         1         1         t        x         u         0*         1              0                3              1           0            z-x          wx                                                     1  
%         2         1         t        y         u         0          1              1                3              2          y-x           z-x          wy      v                                              1
%         3         1         t        z         u         0          1              2                3              3          z-x           z-x          wz      v                                              1
%         4         1         t        a         u         1          2              0                3              1           0            c-a          wa      v                                              1 
%         5         1         t        b         u         0          2              1                3              2          b-a           c-a          wb      v                                              1
%         6         1         t        c         u         0          2              2                3              3          c-a           c-a          wc      v                                              1
%         7         1         t        q         u         1          3              0                3              1           0            s-q          wq      v                                              1
%         8         1         t        r         u         0          3              1                3              2          r-q           s-q          wr      v                                              1
%         9         1         t        s         u         0          3              2                3              3          s-q           s-q          ws      v                                              1
%         10        1         t        j         u         1          4              0                0              1           0             0           wj      v                                              1  



%     where,
%                  row       =   row number, obvi
%        1.        track     =   identifies track 
%        2.        time      =   all timepoints associated with concatinated length trajectories
%        3.        majAx     =   length values from concatentated length trajectories
%        4.        mu        =   calculated growth rates from SlidingFits.m
%        5.        drop?     =   finding where individual cell cycles start and end, a boolean
%        6.        curve     =   an id number for each individual cell cycle
%        7.        tSince    =   time since birth
%        8.        Duration  =   full duration of cell cycle pertaining to current row
%        9.        stage     =   time since birth / duration of entire cycle
%       10.        lAdded    =   increments of added length since time of birth
%       11.        addedL    =   total length added during cell cycle pertaining to current row
%       12.        width     =   width values
%       13.        v_cyl     =   volume approximated as a cylinder
%       14.        v_elspe   =   volume approximated as an ellipse
%       15.        v_anu     =   volume approximated as a cylinder with half sphere caps
%       16.        mu_vc     =   rate of doubling vol as cylinder
%       17.        mu_ve     =   rate of doubling vol as ellipse
%       18.        mu_va     =   rate of doubling col as cylinder with half sphere caps
%       19.        addedVC   =   volume (cylindrical) added since birth
%       20.        addedVE   =   volume (ellipsoidal) added since birth
%       21.        addedVA   =   volume (cylindrical with caps) added since birth
%       22.        tot_addedVC  =  volume (cylindrical) added during cell cycle
%       23.        tot_addedVE  =  volume (ellipsoidal) added during cell cycle
%       24.        tot_addedVA  =  volume (capped cylinder) added during cell cycle
%       25.        addedVC_incr  =  instantaneous added volume (cylindrical)    
%       26.        addedVE_incr  =  instantaneous added volume (ellipsoidal)
%       27.        addedVA_incr  =  instantaneous added volume (capped cylinder)
%       28.        x         =   x position in original image
%       29.        y         =   y position in original image
%       30.        frame     =   frame # in original image sequence
%       31.        n         =   movie number in original ND file
%       32.        ecc       =   eccentricity of fitted ellipse
%       33.        angle     =   as in matrix D
%       34.        track(m)  =   not based on ID, but rather linked entity
%       35.        condition =   1 fluc; 2 low; 3 ave; 4 high, unless noted otherwise





% Strategy (for determining cell cycle stage):
%
%        1.  for each curve, determine duration (time)
%        2.  for each time step, determine absolute time since birth
%        3.  for each data point in vector, record as fraction:
%                
%               ccStage = time since birth / total curve duration



% OK! Lez go!

%%
%   Initialize.
clear
load('letstry-2017-06-12-autoTrimmed-scrambled-proportional.mat','Scram6','T','rejectD_scram');
D7 = Scram6;
%D7 = M6;

%clear D2 D3 D4 D5 D6 M6 rejectD;


%%
%   Part One.
%   Assemble the ultimate data matrix!
%


% Initialize data vectors for concatenation

condVals = [];

trackID = [];
trackM = [];

Time = [];

x_pos = [];
y_pos = [];
orig_frame = [];
movie_num = [];
eccentricity = [];
angle = [];

lengthVals = [];
widthVals = [];
vcVals = [];
veVals = [];
vaVals = [];

muVals = [];
mu_vcVals = [];
mu_veVals = [];
mu_vaVals = [];

isDrop = []; 
dropThreshold = -0.75;                                                     % consider greater negatives a division event

curveFinder = [];                                                        

timeSinceBirth = [];
lengthAdded_incremental_sinceBirth = [];
vcAdded_incremental_sinceBirth = [];
veAdded_incremental_sinceBirth = [];
vaAdded_incremental_sinceBirth = [];

allDurations = [];
allDeltas = [];
allTimestamps = [];

birthSizes = [];
birthTimes = [];

curveDurations = [];
addedLength = [];
addedVC = [];
addedVE = [];
addedVA = [];

addedLength_incremental = [];
addedVC_incremental = [];
addedVE_incremental = [];
addedVA_incremental = [];


%%
% Select xy positions for analysis / concatenation

for n = 1:length(D7)
    
    
    for m = 1:length(D7{n})                                                % use length of growth rate data as it is
                                                                           % slightly truncated from full length track due
                                                                           % to sliding fit
                                                                           
        %   TRACK #                                                        
        lengthCurrentTrack = length(D7{n}(m).TrackID);
        Track = D7{n}(m).TrackID;
        trackID = [trackID; Track];
        
        
        %   frame number in original image
        frameTrack = D7{n}(m).Frame;%(7:lengthCurrentTrack+6);
        orig_frame = [orig_frame; frameTrack];
        
        
        %   TIME
        %timeTrack = T(3:lengthCurrentTrack+2,n)/(60*60);                  % collect timestamp (hr)
        timeTrack = T{n}(frameTrack(1):lengthCurrentTrack+frameTrack(1)-1);%(7:lengthCurrentTrack+6)./(3600);                  % data format, if all ND2s were processed individually
        Time = [Time; timeTrack];                                          % concat=enate timestamp
        
        
        
        %   lengths
        lengthTrack = D7{n}(m).MajAx;%(7:lengthCurrentTrack+6);              % collect lengths (um)
        lengthVals = [lengthVals; lengthTrack];                            % concatenate lengths
        dLengths = diff(lengthTrack);
        dLengths = [0; dLengths];
        addedLength_incremental = [addedLength_incremental; dLengths];
        
        
        %   widths
        widthTrack = D7{n}(m).MinAx;%(7:lengthCurrentTrack+6);               % collect widths (um)
        widthVals = [widthVals; widthTrack];                               % concatenate widths
        
        
        %   x positions in original image
        xTrack = D7{n}(m).X;%(7:lengthCurrentTrack+6); 
        x_pos = [x_pos; xTrack];
        
        
        %   y positions in original image
        yTrack = D7{n}(m).Y;%(7:lengthCurrentTrack+6);
        y_pos = [y_pos; yTrack];
        
        
        %   movie number in original ND2
        movieTrack = ones(length(Track),1)*n;
        movie_num = [movie_num; movieTrack];
        
        
        %   eccentricity of ellipses used in particle tracking
        eccTrack = D7{n}(m).Ecc;%(7:lengthCurrentTrack+6);
        eccentricity = [eccentricity; eccTrack];
        
        
        %   angle of ellipses used in particle tracking
        angTrack = D7{n}(m).Ang;%(7:lengthCurrentTrack+6);
        angle = [angle; angTrack];
        
        %   ELONGATION RATE
        %muTrack = D7{n}(m).Parameters_L(:,1);                                % collect elongation rates (1/hr)
        %muVals = [muVals; muTrack];                                        % concatenate growth rates
        
        
        %   VOLUME
       % v_cylinder = pi * lengthTrack .* (widthTrack/2).^2;                % approx. volume as a cylinder
       % v_ellipse = 4/3 * pi * lengthTrack/2 .* (widthTrack/2).^2;         % approx. volume as an ellipse
       % vol_smallCylinder = (pi * (widthTrack/2).^2 .* (lengthTrack - widthTrack) );
       % vol_sphere = 4/3 * pi * (widthTrack/2).^3;
       % v_anupam = vol_smallCylinder + vol_sphere;                          % approx. volume as cylinder with spherical caps
        
%         
%         vcVals = [vcVals; v_cylinder];                                     % concatenate values
%         veVals = [veVals; v_ellipse];
%         vaVals = [vaVals; v_anupam];
%         
%         dVC = diff(v_cylinder);
%         dVC = [0; dVC];
%         addedVC_incremental = [addedVC_incremental; dVC];
%         
%         dVE = diff(v_ellipse);
%         dVE = [0; dVE];
%         addedVE_incremental = [addedVE_incremental; dVE];
%         
%         dVA = diff(v_anupam);
%         dVA = [0; dVA];
%         addedVA_incremental = [addedVA_incremental; dVA];
%         
%         
%         %   GROWTH RATE (VOLUME)
%         mu_vcTrack = D7{n}(m).Parameters_VC(:,1);                          % as approximated by a cylinder
%         mu_vcVals = [mu_vcVals; mu_vcTrack];                               % see slidingFits
%         
%         mu_veTrack = D7{n}(m).Parameters_VE(:,1);                          % as approximated by an ellipse
%         mu_veVals = [mu_veVals; mu_veTrack];                               % see slidingFits
%         
%         mu_vaTrack = D7{n}(m).Parameters_VE(:,1);                          % as approximated by anupam's suggestion
%         mu_vaVals = [mu_vaVals; mu_vaTrack];       
%         
%         %   DROP?
%         dropTrack = diff(lengthTrack);
%         toBool = dropTrack < dropThreshold;                                % converts different to a Boolean based on dropThreshold
%         toBool = [0; toBool];                                              % * add zero to front, to even track lengths
%         isDrop = [isDrop; toBool];
%         
%         
%         
%         %   CURVE FINDER                                                  
%         numberFullCurves = sum(toBool) - 1;                                      
%         curveTrack = zeros(length(toBool),1);
%         curveCounter = 0;                                                  % finds and labels full curves within a single track
%                                                                            % hint: full curves are bounded by ones
%         for i = 1:length(toBool) 
%             if toBool(i) == 0                                              % 1. disregard incomplete first curve
%                 curveTrack(i,1) = curveCounter;                            %    by starting count at 0   
%             elseif (toBool(i) == 1)
%                 curveCounter = curveCounter + 1;                           % 2. how to disregard final incomplete segment? 
%                 if curveCounter <= numberFullCurves                              %    stop when curveCount exceeds number of fullCurves
%                     curveTrack(i,1) = curveCounter;
%                 else                                                       % all incomplete curves are filled with 0
%                     break                                                  
%                 end
%             end
%         end
%         %   generate column that identifies each tpt to a curve in given track
%         %   i.e. count starts over with each new track
%         curveFinder = [curveFinder; curveTrack];
%         
%         
%         
%         %   TIME SINCE BIRTH
%         
%         % 1. find timepoints with division/birth events
%         % 2. calculate, for each timestep: 
%         %       a)  time since birth
%         %       b)  added mass since birth
%         %       c)  final timestep is also cycle duration & added size 
%         
%         
%         % Part C.
%         % generate a vector of timestamps where events occur
%         isolateEvents = timeTrack.*toBool;
%         eventTimes = isolateEvents(isolateEvents~=0);                      
%         
%         % Part C.
%         % preparing to collect duration and size at the end of each full cell cycle (curve)             
%         durationsPerTrack = zeros(numberFullCurves,1);
%         lengthPerTrack = zeros(numberFullCurves,1);
%         vcPerTrack = zeros(numberFullCurves,1);
%         vePerTrack = zeros(numberFullCurves,1);
%         
%         % Part A & B.
%         % preparing to collect incremental time and mass for all timesteps
%         tsbPerTrack = zeros(lengthCurrentTrack,1);
%         lsbPerTrack = zeros(lengthCurrentTrack,1);
%         vcsbPerTrack = zeros(lengthCurrentTrack,1);
%         vesbPerTrack = zeros(lengthCurrentTrack,1);
%         vasbPerTrack = zeros(lengthCurrentTrack,1);
%         
%         % per individual curve
%         %       - identify current birth event and division event (i.e. next birth)
%         %       - isolate timepoints in between events for calculations specific to that curve
%         %       - time since birth = isolated timepoints minus time of birth
%         %       - added mass since birth = length(at timepoints) minus length at birth 
%         %       - final added mass and total curve duration = values at division event
%        
%         for currentCurve = 1:numberFullCurves; 
%             
%             % identify events bounding each curve
%             currentBirthRow = find(timeTrack == eventTimes(currentCurve)); 
%             nextBirthRow = find(timeTrack == eventTimes(currentCurve+1));  
%             currentTimes = timeTrack(currentBirthRow:nextBirthRow-1);
%             
%             % incremental time
%             tsbPerCurve = currentTimes - timeTrack(currentBirthRow);
%             tsbPerTrack(currentBirthRow:nextBirthRow-1,1) = tsbPerCurve;
%                       
%             % incremental length
%             lsbPerCurve = lengthTrack(currentBirthRow:nextBirthRow-1) - lengthTrack(currentBirthRow);
%             lsbPerTrack(currentBirthRow:nextBirthRow-1,1) = lsbPerCurve;
%             
%             % incremental volume (cylindrical)
%             vcsbPerCurve = v_cylinder(currentBirthRow:nextBirthRow-1) - v_cylinder(currentBirthRow);
%             vcsbPerTrack(currentBirthRow:nextBirthRow-1,1) = vcsbPerCurve;
%             
%             % incremental volume (ellipsoidal)
%             vesbPerCurve = v_ellipse(currentBirthRow:nextBirthRow-1) - v_ellipse(currentBirthRow);
%             vesbPerTrack(currentBirthRow:nextBirthRow-1,1) = vesbPerCurve;
%             
%             % incremental volume (anupam's method)
%             vasbPerCurve = v_anupam(currentBirthRow:nextBirthRow-1) - v_anupam(currentBirthRow);
%             vasbPerTrack(currentBirthRow:nextBirthRow-1,1) = vasbPerCurve;
%             
%             % final duration and mass
%             durationsPerTrack(currentCurve) = tsbPerCurve(end);            % tsb = time since brith
%             lengthPerTrack(currentCurve) = lsbPerCurve(end);               % lsb = length added since birth
%             vcPerTrack(currentCurve) = vcsbPerCurve(end);                  % vcsb = volume added since birth
%             vePerTrack(currentCurve) = vesbPerCurve(end);
%             vaPerTrack(currentCurve) = vasbPerCurve(end);
%             
%         end
%         
%         
%         %   SPIN-OFF DATA COMPILATIONS (two groups):
%         %       
%         %       1. "all" group:
%         %               - all cell cycle durations
%         %               - all added masses since birth
%         %               - all corresponding timestamps per cell cycle (end)
%         %
%         %       2. "birth" group
%         %               - birth lengths
%         %               - birth timestamps
%         %
%         
%         % "ALL" group
%         timeSinceBirth = [timeSinceBirth; tsbPerTrack]; % compiled values of time passed
%         allDurations = [allDurations; durationsPerTrack]; % compiled final cell cycle durations
%         
%         lengthAdded_incremental_sinceBirth = [lengthAdded_incremental_sinceBirth; lsbPerTrack]; % compiled increments of added length
%         vcAdded_incremental_sinceBirth = [vcAdded_incremental_sinceBirth; vcsbPerTrack];
%         veAdded_incremental_sinceBirth = [veAdded_incremental_sinceBirth; vesbPerTrack];
%         vaAdded_incremental_sinceBirth = [vaAdded_incremental_sinceBirth; vasbPerTrack];
%         
%         allDeltas = [allDeltas; lengthPerTrack]; % compiled final added mass per cell cycle
%         
%         if length(eventTimes) > 1
%             allTimestamps = [allTimestamps; eventTimes(2:end)]; % compiled timestamps for FULL cell cycles
%         end
%         
%         % "BIRTH" group
%         birthLengths = lengthTrack.*toBool; % isolate lengths at birth
%         birthRows = find(birthLengths > 0); % when drop=1, then that is the start of new curve
%         sizeAtBirth = birthLengths(birthRows);
%         birthSizes = [birthSizes; sizeAtBirth];
%         
%         timeAtBirth = timeTrack(birthRows);
%         birthTimes = [birthTimes; timeAtBirth];
%         
% 
%         
%         %   CURVE DURATION & ADDED LENGTH
%         
%         %   for calculations of cell cycle fraction, etc, generate:
%         %           1.  a vector of total cell cycle duration
%         %           2.  a vector of final length added in that cell cycle
%         %   compile individual curve durations in single vector
%         perTrack_duration = zeros(lengthCurrentTrack,1);
%         perTrack_length = zeros(lengthCurrentTrack,1);
%         perTrack_vc = zeros(lengthCurrentTrack,1);
%         perTrack_ve = zeros(lengthCurrentTrack,1);
%         perTrack_va = zeros(lengthCurrentTrack,1);
%         
%         
%         % for all timepoints in current track:           
%         %       - if timepoint is part of a full curve, move on so value remains zero
%         %       - if timepoint is part of a full curve, record final curve duration and added size 
%         
%         for j = 1:length(curveTrack) 
%             if curveTrack(j) == 0 
%                 continue
%             else
%                 perTrack_duration(j,1) = durationsPerTrack(curveTrack(j));
%                 perTrack_length(j,1) = lengthPerTrack(curveTrack(j));
%                 perTrack_vc(j,1) = vcPerTrack(curveTrack(j));
%                 perTrack_ve(j,1) = vePerTrack(curveTrack(j));
%                 perTrack_va(j,1) = vaPerTrack(curveTrack(j));
%             end
%         end
%         curveDurations = [curveDurations; perTrack_duration]; % collect all durations for analytical ease (ccStage)
%         addedLength = [addedLength; perTrack_length];
%         addedVC = [addedVC; perTrack_vc];
%         addedVE = [addedVE; perTrack_ve];
%         addedVA = [addedVA; perTrack_va];
%         
%         
%         %    CELL CYCLE FRACTION
%         
%         %   cc fraction = time since birth / total curve duration
%         ccFraction = timeSinceBirth./curveDurations;                       % NaN =  no full cycle
%                                                                            % 0   =  start of full cycle
%                                                                            % 1   =  end of full cycle
%     
                                                                           
        %   CONDITION
        % assign condition based on xy number
        condition = ceil(n/10);
        
        % label each row with a condition #
        condTrack = ones(lengthCurrentTrack,1)*condition;
        condVals = [condVals; condTrack];
        
    end % for m
    
    disp(['Tracks (', num2str(m), ') complete from xy ', num2str(n), ', condition ', num2str(condition), '!'])
    
    
end % for n

%%

vcVals = NaN(length(angle),1);
veVals = NaN(length(angle),1);
vaVals = NaN(length(angle),1);

muVals = NaN(length(angle),1);
mu_vcVals = NaN(length(angle),1);
mu_veVals = NaN(length(angle),1);
mu_vaVals = NaN(length(angle),1);

isDrop = NaN(length(angle),1);  
curveFinder = NaN(length(angle),1);                                                      

timeSinceBirth = NaN(length(angle),1);
lengthAdded_incremental_sinceBirth = NaN(length(angle),1);
vcAdded_incremental_sinceBirth = NaN(length(angle),1);
veAdded_incremental_sinceBirth = NaN(length(angle),1);
vaAdded_incremental_sinceBirth = NaN(length(angle),1);

allDurations = NaN(length(angle),1);
allDeltas = NaN(length(angle),1);
allTimestamps = NaN(length(angle),1);

birthSizes = NaN(length(angle),1);
birthTimes = NaN(length(angle),1);

curveDurations = NaN(length(angle),1);
addedLength = NaN(length(angle),1);
addedVC = NaN(length(angle),1);
addedVE = NaN(length(angle),1);
addedVA = NaN(length(angle),1);

addedLength_incremental = NaN(length(angle),1);
addedVC_incremental = NaN(length(angle),1);
addedVE_incremental = NaN(length(angle),1);
addedVA_incremental = NaN(length(angle),1);

ccFraction = NaN(length(angle),1);


%%
% Compile data into single matrix
dataMatrix = [trackID Time lengthVals muVals isDrop curveFinder timeSinceBirth curveDurations ccFraction lengthAdded_incremental_sinceBirth addedLength widthVals vcVals veVals vaVals mu_vcVals mu_veVals mu_vaVals vcAdded_incremental_sinceBirth veAdded_incremental_sinceBirth vaAdded_incremental_sinceBirth addedVC addedVE addedVA addedVC_incremental addedVE_incremental addedVA_incremental x_pos y_pos orig_frame movie_num eccentricity angle trackM condVals];

save('dm-2017-06-12_scrambled_proportional.mat', 'dataMatrix');
