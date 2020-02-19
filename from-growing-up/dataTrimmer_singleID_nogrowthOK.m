%% dataTrimmer


%  Goal: automatedly, remove or clip tracks that are unlikely to be well-tracked cells.
%         Removal or trimming is as specified by the selection criteria below.



%  Selection Criteria:

%   1. Track must have only one TrackID
%           - tracks are clipped, such that data prior to change in TrackID number remain

%   2. Tracks must be at least the length of fitting window
%           - note: currently, windowSize = 5 frames

%   3. Tracks cannot oscillate too quickly between gains and losses

%   4. Tracks that do not increase by more than JumpFrac (30% of size at previous timepoint)
%           - tracks are clipped, such that data prior to jump remain

%   5. Tracks must be at least the length of fitting window
%           - note: currently, windowSize = 5 frames
%           -       this criteria is repeated to follow another clipping step

%   6. Tracks must be of reasonable size, at least SizeStrainer (1.5um)




% last edit: July 3, 2017 during intense image processing overhaul with
%            confusing glucose data

% retired: 2018 March 29

%% initialize

% particle tracking data
clear
load('letstry-2017-06-12-dSmash.mat');
D = D_smash;

% reject data matrix
rejectD = cell(6,length(D));

% criteria counter
criteria_counter = 0;


%% Criteria ONE: tracks cannot contain multiple TrackIDs

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
% 6. when all tracks finished, save accumulated rejects
% 7. report!
% 8. repeat for next movie


% for easy reshuffling of criteria
criteria_counter = criteria_counter + 1;

for n = 1:length(D)
    
    
    % 0. initialize
    data = D{n};
    
    % 0. remove 'Conversion' field, as it is only one element and interferes with downstream clipping.
    data = rmfield(data,'Conv'); 
    
    
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
    
    % 7. report !
    disp(strcat('Clipping (', num2str(length(currentRejects)),') tracks with multiple IDs from xy (', num2str(n),') !'))
    
    %8. repeat for all movies
    clear currentRejects data rejectTrack rejectIDs originalIDs originalTrack
    clear firstTrack firstIDs reject_counter
end

clear isChange trackIDs m n;

%%
% double-check:
% how many TrackID(1)s are represented only with < 5 data points?
% see dynamicOutlines commit 2017-07-03 for example movie with n=52
% and movie: shortTracks-afterTrackID-clipping-5fps.avi
% justifying the removal of these very short tracks

% for n = 1:length(D2)
%     
%     currentMovie = D2{n};
%     
%     for m = 1:length(currentMovie)
%         currentTrack = currentMovie(m);
%         trackLengths(m,1) = length(currentTrack.X);
%     end
%     shorties = find(trackLengths < 5);
%     allShorts{n} = shorties;
% end


%% Criteria Two: tracks must be at least window size in length (5 frames)


% Goal: after clipping tracks by track ID, some are quite short and lead to problems in
%       subsequent steps, i.e. those that require rates of change produce
%       an error if tracks are only 1 frame long. Remove these.

% 0. initiaize new dataset before trimming
% 0. in current movie
%           1. for each track, determine number of timepoints
%           2. find tracks that are shorter than threshold number of frames
%           3. report!
%           4. if no sub-threshold tracks, continue to next movie
%           5. else, remove structures based on row # (in reverse order)
%           6. save sub-threshold tracks into reject data matrix
% 7. repeat for next movie



criteria_counter = criteria_counter + 1;


% 0. initialize new dataset before trimming
D3 = D2;
windowSize = 5;                                                             % each timepoint = 1:05 mins;

for n = 1:length(D);

    % 1. determine number of timepoints in each track m 
    for m = 1:length(D3{n})
        numFrames(m) = length(D3{n}(m).MajAx);
    end
    
    % 2. find tracks that are shorter than threshold number of frames
    subThreshold = find(numFrames < windowSize);       
    
    % 3. report!
    X = ['Removing ', num2str(length(subThreshold)), ' short tracks from D3(', num2str(n), ')...'];
    disp(X)
    
    % 4. to that loop doesn't crash if nothing is too short
    if isempty(subThreshold) == 1
        continue
    end
    
    % 5. remove structures based on row # (in reverse order)
    fingers = 0;
    for toRemove = 1:length(subThreshold)
        
        r = length(subThreshold) - fingers;                  % reverse order                  
        D3{n}(subThreshold(r)) = [];                         % deletes data
        tracks_shortGlimpses(r,1) = D{n}(subThreshold(r));   % store data for reject data matrix
        fingers = fingers + 1;
        
    end
    
    % 6. save sub-threshold tracks into reject data matrix
    rejectD{criteria_counter,n} = tracks_shortGlimpses;
    
    % 7. repeat for all movies
    clear  numFrames m subThreshold fingers toRemove r X tracks_shortGlimpses;
    
end

 clear windowSize n; 
 

%% Criteria Three: tracks cannot oscillate too quickly between gains and losses


% Goal: jiggly tracks correspond to non-growing particles or errors in tracking...remove!
%       method: are most negatives divisions?

% 0. initialize copy of dataset before trimming
% 0. define thresholds
% 0. for each movie
%        1. for each track, collect % of non-drop negatives per track length (in frames)
%                2. isolate length data from current track
%                3. find change in length between frames
%                4. convert change into binary, where positives = 0 and negatives = 1
%                5. find all drops (negaitves that exceed drop threshold)
%                6. find the ratio of non-drop negatives per track length
%                7. store ratio for subsequent removal
%       8. repeat for all tracks
%       9. determine which tracks fall under jiggle threshold
%      10. report!
%      11. remove jiggly tracks in reverse order
% 12. repeat for all movies


criteria_counter = criteria_counter + 1;

% 0. initiaize new dataset before trimming
D4 = D3;

% 0. define threshold change in length considered a division
dropThreshold = -0.75;

% 0. define threshold under which tracks are too jiggly
jiggleThreshold = -0.1;

for n = 1:length(D)
    
    % 1. for each track, collect % of non-drop negatives per track length (in frames)
    nonDropRatio = NaN(length(D{n}),1);
    
    for m = 1:length(D{n})
        
        % 2. isolate length data from current track
        lengthTrack = D{n}(m).MajAx;
        
        % 3. find change in length between frames
        diffTrack = diff(lengthTrack);
        
        % 4. convert change into binary, where positives = 0 and negatives = 1
        binaryTrack = logical(diffTrack < 0);
        
        % 5. find all drops (negatives that exceed drop threshold)
        dropTrack = diffTrack < dropThreshold;
        
        % 6. find the ratio of non-drop negatives per track length
        nonDropNegs = sum(dropTrack - binaryTrack);
        squiggleFactor = nonDropNegs/length(lengthTrack);
        
        % 7. store ratio for subsequent removal
        nonDropRatio(m) = squiggleFactor;
    
    % 8. repeat for all tracks
    end
    
    % 9. determine which tracks fall under jiggle threshold
    belowThreshold = find(nonDropRatio < -0.1);
    
    % 10. report!
    X = ['Removing ', num2str(length(belowThreshold)), ' jiggly tracks from D4(', num2str(n), ')...'];
    disp(X)
    
    % 11. remove jiggly structures based on row # (in reverse order)
    counter = 0;
    for toRemove = 1:length(belowThreshold)
        
        r = length(belowThreshold) - counter;                  % reverse order
        D4{n}(belowThreshold(r)) = [];                         % deletes data
        jigglers(r,1) = D{n}(belowThreshold(r));   % store data for reject data matrix
        counter = counter + 1;
        
    end

    % 13. save sub-threshold tracks into reject data matrix
    rejectD{criteria_counter,n} = jigglers;
    
    % 14. repeat for all movies
    clear nonDropRatio lengthTrack diffTrack dropTrack nonDropNegs squiggleFactor belowThreshold
  
end

clear n gainLossRatio;




%% Criteria Four: clip tracks to remove >30% jumps in cell size


% Goal: huge jumps in cell size correspond to.... thus, clip off!

%  0. initialize threshold value
%  0. isolate data for current movie
%        0. for each track, search for large positive increases: > %30 jumps in size
%               1. determine change in length between each timestep
%               2. express change of as a fraction of cell size in previous timestep
%               3. list rows in which the change exceeds size jump threshold
%               4. if the track contains jumps,
%                        i. isolate structure of target track z, in prep to clip all variables (MajAx, X, Y, etc.)
%                       ii. define timepoint of earliest jump
%                      iii. clip data structure at earliest jump
%                       iv. redefine track in data set as clipped structure
%                        v. store remainder of original Target for reject data matrix
%        5. when all tracks finished, save trimmed data and accmulated rejects
%        6. report number of clipped tracks
%  7. repeat for next movie



criteria_counter = criteria_counter + 1;

% 0. initialize
jumpFrac = 0.3;                                                            
                                                                           
for n = 1:length(D);                                                       
    
    % 0. isolate data for current movie
    data = D4{n};
    
    jump_counter = 0;
    for m = 1:length(data)   
                                                                           
        % 1. determine change in length between each timestep                           
        Rates = diff(data(m).MajAx);

        
        % 2. express change of as a fraction of the cell size in previous timestep
        Lengths = data(m).MajAx(1:length(Rates));
        growthFrac = Rates./Lengths;
                                                                           
        % 3. list rows in which the change exceeds size jump threshold
        jumpPoints = find(growthFrac > jumpFrac);                        

        
        % 4. if the track contains jumps...
        if isempty(jumpPoints) == 0
            
            % i. isolate structure of target track z, in prep to clip all variables (MajAx, X, Y, etc.)
            originalTrack = data(m);
            
            % ii. define timepoint of earliest jump
            clipPoint = jumpPoints(1);
            
            % iii. clip data structure at desired timepoint
            clippedTarget = structfun(@(M) M(1:clipPoint), originalTrack, 'Uniform', 0);
            
            % iv. redefine track in data set as clipped structure
            data(m) = clippedTarget;
            
            % v. store remainder of original Target for reject data matrix
            remainderTrack = structfun(@(M) M(clipPoint+1:end), originalTrack, 'Uniform', 0);
            
            jump_counter = jump_counter + 1;
            trackScraps(jump_counter,1) = remainderTrack;
            
        end
        
    end
    
    % 5. when all tracks finished, save trimmed data and accmulated rejects
    D5{n} = data;
    rejectD{criteria_counter,n} = trackScraps;
    
    % 6. report!
    X = ['Clipping ', num2str(jump_counter), ' jumps in D5(', num2str(n), ')...'];
    disp(X)
    
    
end
    
clear growthFrac Lengths jumpFrac jump_counter Rates clippedTarget clipPoint m n;
clear trackScraps remainderTrack originalTrack jumpPoints X data;




%% Criteria Five: repeat - tracks must be at least window size in length (5 frames)


% Goal: after clipping off large size jumps, some tracks can be quite short
%       Remove tracks shorter than windowSize, such that mu calculations in
%       slidingFits does not encounter complications

% 0. initiaize new dataset before trimming
% 0. in current movie
%           1. for each track, determine number of timepoints
%           2. find tracks that are shorter than threshold number of frames
%           3. report!
%           4. if no sub-threshold tracks, continue to next movie
%           5. else, remove structures based on row # (in reverse order)
%           6. save sub-threshold tracks into reject data matrix
% 7. repeat for next movie



criteria_counter = criteria_counter + 1;


% 0. initialize new dataset before trimming
D6 = D5;
windowSize = 5;                                                             % each timepoint = 1:05 mins;

for n = 1:length(D);

    % 1. determine number of timepoints in each track m 
    for m = 1:length(D6{n})
        numFrames(m) = length(D6{n}(m).MajAx);
    end
    
    % 2. find tracks that are shorter than threshold number of frames
    subThreshold = find(numFrames < windowSize);       
    
    % 3. report!
    X = ['Removing ', num2str(length(subThreshold)), ' short tracks from D6(', num2str(n), ')...'];
    disp(X)
    
    % 4. to that loop doesn't crash if nothing is too short
    if isempty(subThreshold) == 1
        continue
    end
    
    % 5. remove structures based on row # (in reverse order)
    fingers = 0;
    for toRemove = 1:length(subThreshold)
        
        r = length(subThreshold) - fingers;                  % reverse order                  
        D6{n}(subThreshold(r)) = [];                         % deletes data
        rejectTracks(r,1) = D{n}(subThreshold(r));   % store data for reject data matrix
        fingers = fingers + 1;
        
    end
    
    % 6. save sub-threshold tracks into reject data matrix
    rejectD{criteria_counter,n} = rejectTracks;
    
    % 7. repeat for all movies
    clear  numFrames m subThreshold fingers toRemove r X rejectTracks;
    
end

clear n; 


 
%% Criteria Six: maximum particle size must be greater than 1.5um

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
    rejectD{6,n} = tracks_tooSmalls;
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
    rejectD{7,n} = currentRejects;
    
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



