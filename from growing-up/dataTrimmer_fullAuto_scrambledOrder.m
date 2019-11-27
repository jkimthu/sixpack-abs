%% dataTrimmer_fullAuto_scrambledOrder


%  Goal: automatedly, remove or snip tracks that are unlikely to be well-tracked cells.
%         Removal or trimming is as specified by the selection criteria below.



%  Selection Criteria:

%   1. Track must have only one TrackID
%           - tracks are snipped, such that data prior to change in TrackID number remain

%   2. Tracks cannot oscillate too quickly between gains and losses

%   3. Tracks that do not increase by more than JumpFrac (30% of size at previous timepoint)
%           - tracks are snipped, such that data prior to jump remain

%   4. Tracks must be at least the length of fitting window
%           - note: currently, windowSize = 5 frames

%   5. Tracks must be of reasonable size, at least SizeStrainer (1.5um)




% last edit: July 2, 2017 during image processing overhaul during crazy
%            confusion over glucose data

% retired: 2018 March 29


%% initialize

% particle tracking data
clear
load('letstry-2017-06-12-dSmash.mat');
D = D_smash;

% reject data matrix
rejectD_scram = cell(5,length(D));

% criteria counter
criteria_counter = 0;


%% Criteria ONE: tracks cannot contain multiple TrackIDs

criteria_counter = criteria_counter + 1;

% Goal: it seems that any tracks with changes in trackID are ones that are poorly joined 
%       but at least the first trackID is useable. Let's keep these first ones
%       and reject data from subsequent IDs.

% 0. for each track in current movie
%       1. determine whether trackID contains number changes
%               2. if so, trim track such that only first trackID remains
%                  (all following are very likely error prone)
%                         i. isolate entire data struct of current track,
%                            in prep to clip all variables (MajAx, X, Y, etc.)
%                        ii. isolate data corresponding to first TrackID
%               3. replace data from original track (containing multiple IDs) with trimmed data
%               4. add remainder of track to temporary (movie-specific) rejects collection
%       5. if no changes, continue to next track
% 6. when all tracks finished, save accmulated rejects.
% 7. repeat for next movie

for n = 1:length(D)
    
    
    % 0. initialize
    data = D{n};
    currentRejects = [];
    reject_counter = 0;
    
    for m = 1:length(data)
        
        % 1. determine whether trackID contains number changes
        trackIDs = data(m).TrackID;
        isChange = diff(trackIDs);
        
        % if so,
        if sum(isChange) ~= 0
            
            % 2. trim track such that only first trackID remains
            reject_counter = reject_counter +1;
            disp(strcat('Track (', num2str(m),') from xy (', num2str(n),') has multiple IDs! Trimming...'))
            
            % i. isolate entire data struct of current track, in prep to clip all variables (MajAx, X, Y, etc.)
            originalTrack = data(m);
            
            % ii. isolate data corresponding to first TrackID
            originalIDs = originalTrack.TrackID;
            firstIDs = originalIDs == originalTrack.TrackID(1);
            firstTrack = structfun(@(M) M(firstIDs), originalTrack, 'Uniform', 0);
            
            
            % 3. replace data from original track (containing multiple IDs) with trimmed data
            data(m) = firstTrack;
            
            
            % 4. add remainder of track to rejects collection
            rejectIDs = originalIDs ~= originalTrack.TrackID(1);
            rejectTrack = structfun(@(M) M(rejectIDs), originalTrack, 'Uniform', 0);
            currentRejects{reject_counter} = rejectTrack;
            
            % 5. if no changes, continue to next track
        end
        
    end
    
    % 6. when all tracks finished, save trimmed data and accmulated rejects
    D2{n} = data;
    rejectD{criteria_counter,n} = currentRejects;
    
    clear currentRejects data rejectTrack rejectIDs originalIDs originalTrack
    clear firstTrack firstIDs
end


%% Criteria Two: tracks cannot oscillate too quickly between gains and losses

criteria_counter = criteria_counter + 1;

D3 = D2;
gainLossRatio = 0.85;


for n = 1:length(D)
    
    for i = 1:length(D3{n})
        
        % determine derivative of each timestep in length
        Signs = diff(D3{n}(i).MajAx);
        
        % minute res is so noisy. use average change for every 10 steps.
        sampleLength = floor( length(Signs)/10 ) * 10;
        
        Signs = mean(reshape(Signs(1:sampleLength),10,[]))';
        Signs(Signs<0) = 0;
        Signs(Signs>0) = 1;
        
        % determine ratio of negatives to positives
        trackRatio = sum(Signs)/length(Signs);
        
        % store ratio for track removal
        allRatios(i,1) = trackRatio;
        
        clear Signs trackRatio sampleLength;
    end
    
    
    % identify tracks with over 15% negatives
    swigglyIDs = find(allRatios < 0.85);
    
    % report!
    X = ['Removing ', num2str(length(swigglyIDs)), ' swiggly tracks from D3(', num2str(n), ')...'];
    disp(X)
    
    % so loop doesn't crash if nothing is too wiggly
    if isempty(swigglyIDs) == 1
        continue
    end
    
    % remove structures based on row # (in reverse order)
    counter = 0;
    for q = 1:length(swigglyIDs)
        r = length(swigglyIDs) - counter;                                
        D3{n}(swigglyIDs(r)) = [];
        swigglyTracks(r,1) = D{n}(swigglyIDs(r));   % recording to add into reject data matrix  
        counter = counter + 1;
    end
    
    % save tracks that are too swiggly into reject data matrix
    rejectD_scram{1,n} = swigglyTracks;
    
    clear allRatios allTracks bottomTracks gainLossRatio swigglyIDs swigglyTracks q r i counter X; 
end

clear n gainLossRatio;




%% Criteria Two: clip tracks to remove >30% jumps in cell size
%          
%              - if too positive (cells shouldn't double within three minutes)
%              - in these cases, what causes these large jumps? check!
%              - negatives are OK because cells have to divide!

Scram3 = D3;  
JumpFrac = 0.3;                                                            % JumpFrac = threshold parameter
                                                                           % tracks that increase by a cell size fraction greater than JumpFrac will be eliminated from final dataset
for n = 1:length(D);                                                       
    
    counter = 0;                                                           
    Scram3{n} = rmfield(Scram3{n},'Conv');                                         % removes the 'Conversion' field, as it is only one element. it interferes with downstream clipping.
    
    for i = 1:length(Scram3{n})                                            % error "undefined variable" if some cells have no tracks, []     
                                                                           
        % derivative of current length tracjetory                            
        Rates = diff(Scram3{n}(i).MajAx);
        Lengths = Scram3{n}(i).MajAx(1:length(Rates));
        
        % derivative of as a fraction of cell size in previous timestep
        growthFrac = Rates./Lengths;
                                                                           
        % list rows in which Rate exceeds threshold
        jumpPoints = find(growthFrac > JumpFrac);                        

        
        % for tracks containing jumps...
        if isempty(jumpPoints) == 0
            
            % record track row, such that jumpTrack lists all jumpy tracks
            %currentJumpTrack = i;
            counter = counter + 1;                             
            %jumpTrack(counter) = i;
            
            % loop through and trim all jumpy tracks
            % uh oh- does this continuous loop keep snipping tracks when
            % they don't need to be?
            %for z = 1:length(jumpTrack) 
               
                % isolate structure of target track z, in prep to clip all variables (MajAx, X, Y, etc.)
                target2Snip = Scram3{1,n}(i);                            
                
                % defines timepoint at point of earliest jump
                snipPoint = jumpPoints(1);  
                
                % clips structure at desired timepoint
                snippedTarget = structfun(@(M) M(1:snipPoint), target2Snip, 'Uniform', 0);   
                
                % redefines track in data set as trimmed structure
                Scram3{1,n}(i) = snippedTarget;
                
                % store original Target for reject data matrix
                tracks_clipJump(counter,1) = target2Snip;
            %end
        end
        
    end
    
    % report!
    X = ['Clipping ', num2str(counter), ' jumps from Scram3(', num2str(n), ')...'];
    disp(X)
    
    % save tracks that are too small into reject data matrix
    rejectD_scram{2,n} = tracks_clipJump;
    
    if n < length(D)
        clear Rates jumpTrack z;                                        % erases info from current series, so that tracks don't roll into next iteration
    end
    
end
    
clear tf JumpFrac Target counter jumpTrack Rates Tpt clipTarget i n z X;
clear Jumpers tracks_clipJump;


%% Criteria Three: total track length must be at least size of fitting window

Scram4 = Scram3;
windowSize = 5;                                                             % each timepoint = 1:05 mins;

for n = 1:length(D);

    % determine number of timepoints in each track i 
    for i = 1:length(Scram4{n})
        cellLength{i} = length(Scram4{n}(i).MajAx);
    end
    
    % find tracks that are shorter than ___ mins
    cellLength_dbl = cell2mat(cellLength);
    shortGlimpse = find(cellLength_dbl < windowSize);       
    
    % report!
    X = ['Removing ', num2str(length(shortGlimpse)), ' short tracks from Scram4(', num2str(n), ')...'];
    disp(X)
    
    % to that loop doesn't crash if nothing is too short
    if isempty(shortGlimpse) == 1
        continue
    end
    
    % remove structures based on row # (in reverse order)
    fingers = 0;
    for q = 1:length(shortGlimpse)
        r = length(shortGlimpse) - fingers;                                
        Scram4{n}(shortGlimpse(r)) = [];
        tracks_shortGlimpses(r,1) = D{n}(shortGlimpse(r));   % recording to add into reject data matrix  
        fingers = fingers + 1;
    end
    
    % save tracks that are too small into reject data matrix
    rejectD_scram{3,n} = tracks_shortGlimpses;
    
    clear  cellLength cellLength_dbl i shortGlimpse fingers q r X tracks_shortGlimpses;
    
end

 clear Shortest n; 

 %% Criteria Four: tracks must increase in size by > 30%


Scram5 = Scram4;
GoldenRatio = 1.3;                  

for n = 1:length(D);      
    
    if isempty(Scram5{n}) == 1
        continue
        
    else
        % determine difference percent change between min and max length for each track, i
        for i = 1:length(Scram5{n})
            
            maxLengths(i) = arrayfun(@(Q) max(Q.MajAx), Scram5{n}(i));
            minLengths(i) = arrayfun(@(Q) min(Q.MajAx), Scram5{n}(i));
            
        end
           lengthRatio = maxLengths./minLengths;
    end
    

    
    % find tracks that show little growth                      
    littleGrowth = find(lengthRatio < GoldenRatio); 
    
    % report!
    X = ['Removing ', num2str(length(littleGrowth)), ' non-doublers from Scram5(', num2str(n), ')...']; 
    disp(X)                                                         
    
    % so loop doesn't crash if nothing grows too little
    if isempty(littleGrowth) == 1
        continue
    end
    
    % remove structures based on row # (in reverse order)
    counter = 0;
    for j = 1:length(littleGrowth)
        
        k = length(littleGrowth) - counter;
        Scram5{n}(littleGrowth(k)) = [];
        tracks_littleGrowth(k,1) = D{n}(littleGrowth(k));   % recording to add into reject data matrix
        counter = counter + 1;
        
    end
    
    % save tracks that hardly increase into reject data matrix
    rejectD_scram{4,n} = tracks_littleGrowth;
    clear maxLengths minLengths littleGrowth lengthRatio i j k X counter ToCut tracks_littleGrowth;
    
end

clear GoldenRatio n;
 
%% Criteria Five: maximum particle size must be greater than 1.5um

Scram6 = Scram5;
SizeStrainer = 1.5;

for n = 1:length(Scram5);                           
    
    for i = 1:length(Scram6{n})
        lengthTrack{i} = max(Scram6{n}(i).MajAx);                          
    end          
    
    % finds tracks that don't exceed __ um
    lengthTrack_dbl = cell2mat(lengthTrack);  
    tooSmalls = find(lengthTrack_dbl < SizeStrainer);                          
    
    % report!
    X = ['Removing ', num2str(length(tooSmalls)), ' small particles from Scram6(', num2str(n), ')...'];
    disp(X)
    
    % so loop doesn't crash if nothing is too small
    if isempty(tooSmalls) == 1
        continue
    end
    
    % remove too-small structures based on row # (in reverse order)
    countSmalls = 0;
    for s = 1:length(tooSmalls)
        t = length(tooSmalls) - countSmalls;
        Scram6{n}(tooSmalls(t)) = [];
        tracks_tooSmalls(t,1) = D{n}(tooSmalls(t));      %  recording to add into reject data matrix
        countSmalls = countSmalls + 1;
    end
    
    % save tracks that are too small into reject data matrix
    rejectD_scram{6,n} = tracks_tooSmalls;
    clear lengthTrack lengthTrack_dbl i tooSmalls countSmalls s t X tracks_tooSmalls;

    
end 

clear SizeStrainer n;


 

%% Saving results

save('letstry-2017-06-12-autoTrimmed-scrambled-proportional.mat', 'D_smash', 'D2', 'D3', 'D4', 'D5', 'D6', 'rejectD', 'T')%, 'reader', 'ConversionFactor')


%% dealing with improper track linking
% 
% Goal: it seems that any tracks with changes in trackID are ones that are poorly joined 
%       but at least the first trackID is useable. Let's keep these first ones
%       and reject data from subsequent IDs.

% 0. for each track in current movie
%       1. determine whether trackID contains number changes
%               2. if so, trim track such that only first trackID remains
%                  (all following are very likely error prone)
%                         i. isolate entire data struct of current track,
%                            in prep to clip all variables (MajAx, X, Y, etc.)
%                        ii. isolate data corresponding to first TrackID
%               3. replace data from original track (containing multiple IDs) with trimmed data
%               4. add remainder of track to temporary (movie-specific) rejects collection
%       5. if no changes, continue to next track
% 6. when all tracks finished, save accmulated rejects.
% 7. repeat for next movie


load('letstry-2017-06-12-autoTrimmed-scrambled-proportional.mat');

for n = 1:length(Scram6)
    
    
    % 0. initialize
    data = Scram6{n};
    currentRejects = [];
    reject_counter = 0;
    
    for m = 1:length(data)
        
        % 1. determine whether trackID contains number changes
        trackIDs = data(m).TrackID;
        isChange = diff(trackIDs);
        
        % if so,
        if sum(isChange) ~= 0
            
            % 2. trim track such that only first trackID remains
            reject_counter = reject_counter +1;
            disp(strcat('Track (', num2str(m),') from xy (', num2str(n),') has multiple IDs! Trimming...'))
            
            % i. isolate entire data struct of current track, in prep to clip all variables (MajAx, X, Y, etc.)
            originalTrack = data(m);
            
            % ii. isolate data corresponding to first TrackID
            originalIDs = originalTrack.TrackID;
            firstIDs = originalIDs == originalTrack.TrackID(1);
            firstTrack = structfun(@(M) M(firstIDs), originalTrack, 'Uniform', 0);
            
            
            % 3. replace data from original track (containing multiple IDs) with trimmed data
            data(m) = firstTrack;
            
            
            % 4. add remainder of track to rejects collection
            rejectIDs = originalIDs ~= originalTrack.TrackID(1);
            rejectTrack = structfun(@(M) M(rejectIDs), originalTrack, 'Uniform', 0);
            currentRejects{reject_counter} = rejectTrack;
            
            % 5. if no changes, continue to next track
        end
        
    end
    
    % 6. when all tracks finished, save trimmed data and accmulated rejects
    Scram7{n} = data;
    rejectD_scram{7,n} = currentRejects;
    
    clear currentRejects data rejectTrack rejectIDs originalIDs originalTrack
    clear firstTrack firstIDs
end


%%
i=18;
plot(tracks2Add{i}.Frame,tracks2Add{i}.MajAx,'o')
title(tracks2Add{i}.TrackID(1));


%% visualizing samples of data set

% -- criteria five, check
% select sample of tracks to visualize
bottomTracks = find(allRatios < 0.85);
sampleTracks = find(allRatios < 1);

c = ismember(sampleTracks, bottomTracks);
sampleTracks(c) = [];

for st = 1:length(sampleTracks)
    
    % designate subplot position
    subplot(ceil(length(sampleTracks)/5), 5, st)
    
    % plot
    figure(n)
    plot(T{n}(D6{n}(sampleTracks(st)).Frame(1:end))/3600,(D6{n}(sampleTracks(st)).MajAx),'Linewidth',2)
   
    
    % label
    title(sampleTracks(st));
    
end

%%
% -- final pass, check

%
for n = 52 %:length(D6)
    %counter = counter +1;
    
    for i = 1:20%length(D6{n})
        
        % designate subplot position
        %subplot(ceil(length(D6{n})/5), 5, i)
        subplot(ceil(20/5), 5, i)
        
        % plot
        %figure(counter)
        plot(T{n}(D6{n}(i).Frame(1:end))/3600,(D6{n}(i).MajAx),'Linewidth',2)
        
        % label
        title(i);
        
    end
end

%% visualize any one single curve

%
plot(T{n}(D6{n}(i).Frame(1:end))/3600,(D6{n}(i).MajAx),'Linewidth',2)
        axis([0,10,0,15])



