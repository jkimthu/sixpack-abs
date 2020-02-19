%% squiggleThreshold.m


% Goal: this script visualizes the tracks suriving and rejected by the
%       "squiggle threshold" - a threshold intended to filture out junk
%       tracks from particle tracking data.
%
%       the initial issue, which leads to this separate testing script
%       here, was that fast growing but otherwise perfect tracks (LB)
%       were getting discarded under the following criteria:

%       1. gainLossRatio = 0.85
%       2. ratio calculated from change in size, diff(length),
%          time-averaged across 10 frames (10 diffs)


%       a more recent version was based on the ratio of
%       negative-non-divisions to the length of the track. this worked well
%       in fast growing conditions, when most the negatives are drops. the
%       slower ones however, are more noisy, so these are still harshly
%       lost.


%       latest attempt: a very smoothed sliding fit that should never be
%       negative.




% last edit: jen, 2017 Jul 18

% OK LEZ GO!



%% initialize

% particle tracking data
clear
load('letstry-2017-06-12-dSmash.mat');
D = D_smash;

% reject data matrix
rejectD = cell(6,length(D));

%%

for n = 1 %:length(D6)
    %counter = counter +1;
    
    subplot_counter = 0;
    for i = 61:80%length(D6{n})
        
        subplot_counter = subplot_counter + 1;
        % designate subplot position
        %subplot(ceil(length(D6{n})/5), 5, i)
        subplot(ceil(20/5), 5, subplot_counter)
        
        % plot
        %figure(counter)
        plot(T{n}(D{n}(i).Frame(1:end))/3600,(D{n}(i).MajAx),'Linewidth',2)
        
        % label
        title(i);
        
    end
end


%%

% 0. initiaize new dataset before trimming
slowData = D{1};
fastData = D{52};

%%
i = 200;
figure(1)
plot(slowData(i).MajAx)
uni_slow = unique(slowData(i).TrackID)
slowIDs = slowData(i).TrackID;

figure(2)
plot(fastData(i).MajAx)
uni_fast = unique(fastData(i).TrackID)
fastIDs = fastData(i).TrackID;

%%
% new scheme: are most negatives drops?

% 1. for each track,
%       2. find change in length between frames
%       3. convert change into binary, where positives = 0 and negatives = 1
%       4. find all drops (negaitves that exceed drop threshold)
%       5. find the ratio of non-drop negatives per track length
%       6. store ratio for subsequent removal
%
for n = 1:length(D)
    
    nonDropRatio = NaN(length(D{n}),1);
    dropThreshold = -0.75;
    
    
    for m = 1:length(D{n})
        
        % 1. isolate length data from current track
        lengthTrack = D{n}(m).MajAx;
        
        % 2. find change in length between frames
        diffTrack = diff(lengthTrack);
        
        % 3. convert change into binary, where positives = 0 and negatives = 1
        binaryTrack = logical(diffTrack < 0);
        
        % 4. find all drops (negatives that exceed drop threshold)
        dropTrack = diffTrack < dropThreshold;
        
        % 5. find the ratio of non-drop negatives per track length
        nonDropNegs = sum(dropTrack - binaryTrack);
        squiggleFactor = nonDropNegs/length(lengthTrack);
        
        
        nonDropRatio(m) = squiggleFactor;
        
    end
    
    belowThreshold{n} = find(nonDropRatio < -0.1);
    total2Trim(n) = length(belowThreshold{n});
    
    clear nonDropRatio lengthTrack diffTrack dropTrack nonDropNegs squiggleFactor
end

%%
% find TrackIDs of removed tracks from both methods:

% method 1. SIGNS
data = rejectD{1,1};
signs_IDs_slow = [];

for i = 1:length(data)

    trackIDs = unique(data(i).TrackID);
    signs_IDs_slow = [signs_IDs_slow; trackIDs];

end

data = rejectD{1,52};
signs_IDs_fast = [];

for i = 1:length(data)

    trackIDs = unique(data(i).TrackID);
    signs_IDs_fast = [signs_IDs_fast; trackIDs];

end

%%
% method 2. NONDROPS
nonDrop_IDs_slow = [];
data = belowThreshold{1,1};

for i = 1:length(data)
    
    trackIDs = unique(D{1}(data(i)).TrackID);
    nonDrop_IDs_slow = [nonDrop_IDs_slow; trackIDs];
    
end

nonDrop_IDs_fast = [];
data = belowThreshold{1,52};

for i = 1:length(data)
    
    trackIDs = unique(D{52}(data(i)).TrackID);
    nonDrop_IDs_fast = [nonDrop_IDs_fast; trackIDs];
    
end
%%
% build vectors of IDs in only Signs, only NonDrops, and both

% find IDs rejected in both methods, without repetitions
overLap_slow = intersect(signs_IDs_slow, nonDrop_IDs_slow);
overLap_fast = intersect(signs_IDs_fast, nonDrop_IDs_fast);

% find all unique
unique_rejectIDs_slow = unique([signs_IDs_slow; nonDrop_IDs_slow]);
unique_rejectIDs_fast = unique([signs_IDs_fast; nonDrop_IDs_fast]);

% find IDs ONLY in Signs
rows2keep = find(~ismember(signs_IDs_slow, overLap_slow));
onlySigns_slow = signs_IDs_slow(rows2keep);

rows2keep = find(~ismember(signs_IDs_fast, overLap_fast));
onlySigns_fast = signs_IDs_fast(rows2keep);


% find IDs ONLY in nonDrops
rows2keep = find(~ismember(nonDrop_IDs_slow, overLap_slow));
onlyNonDrops_slow = nonDrop_IDs_slow(rows2keep);

rows2keep = find(~ismember(nonDrop_IDs_fast, overLap_fast));
onlyNonDrops_fast = nonDrop_IDs_fast(rows2keep);

%%


% 2017 jul 18
% super smooth sliding fits:

n=1;
data = D3{n};


% 0. initialize window parameters
windowSize = 5;

% 0. initialize division parameters
dropFrac = -0.3;



numTracks = length(data);

%  2. per track, isolate length, width, frame and time data
for track = 1:numTracks
    
    trackLengths = data(track).MajAx;
    trackFrames = data(track).Frame;
    trackTimes = T{n}(trackFrames)/3600;
    
    
    %  4. build an array that identifies curve #
    
    % 0. initalize array, curveNum
    curveNum = zeros(length(trackFrames),1);
    
    % i. identify all changes in size > threshold
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
    %numWindows = length(trackLengths)/windowSize;
    
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
            multiplier(wCurves == minCurve+1) = 2;
            multiplier(wCurves == minCurve+2) = 4;
            %multiplier(wCurves == minCurve+3) = 8;
            
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
    
    trackData = slidingData(:,1);
    aveMu(track) = mean(trackData);
    
    clear slidingData trackData;
    
    disp(['Track ', num2str(track), ' from xy ', num2str(n), ' complete!'])
    
end

clear currentWindow curveNum dropFrac dropTrack firstWindow growthFrac
clear numTracks numWindows track trackFrames trackLengths trackTimes w windowSize

%%

cutoff = mean(aveMu) - std(aveMu);
isSquiggly = find(aveMu < cutoff);


% 4. if the track contains jumps...
test = data;    

% 3. report!
X = ['Removing ', num2str(length(isSquiggly)), ' short tracks from test(', num2str(n), ')...'];
disp(X)

% 4. to that loop doesn't crash if nothing is too short
%if isempty(isSquiggly) == 1
%    continue
%end

% 5. remove structures based on row # (in reverse order)
squig_counter = 0;
for toRemove = 1:length(isSquiggly)
    
    r = length(isSquiggly) - squig_counter;                  % reverse order
    swigglies(r,1) = data(isSquiggly(r));   % store data for reject data matrix
    test(isSquiggly(r)) = [];                         % deletes data
    squig_counter = squig_counter + 1;
    
end

% 6. save sub-threshold tracks into reject data matrix
reject_test = swigglies;




















