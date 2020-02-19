% checkRejects.m

% Goal: Lots of tracks are getting lost, but are they all that terrible?
%       Visualize reject tracks, identify cause of loss.


% Strategy: 
%           0. Initialize reject data
%           1. Build reject data matrix from movie of interest
%           2. Isolate rejects from trim stage of interest:
%                   - stage 1: excluded from D2, maxL < 1.5 um
%                   - stage 2: excluded from D3, maxL/minL ratio > 1.3
%                   - stage 3: excluded from D4, snip tracks when diff > 0.3
%                   - stage 4: excluded from D5, trackLength < 20 frames
%                   - stage 5: excluded from D6, gainLossRatio < 0.85
%           3. Visualize tracks from stage of interest
%           4. Sort between "definitely trash", and "what happened here"?



% last update: jen, 2017 Jun 27


% OK lez go!!



%%

% 0. initialize reject data
clear
clc
experiment = '2017-06-12';

% open folder for experiment of interest
newFolder = strcat('/Users/jen/Documents/StockerLab/Data/',experiment);
cd(newFolder);

% load data
load('letstry-2017-06-12-autoTrimmed-scrambledOrder-editedJumpTrack.mat');



% 1. Build reject data matrix from movie of interest
n = 52; % movie (xy position) of interest


% build data matrix of rejected tracks
%rejects_currentMovie = rejectD(:,n);
%dm_currentMovie_rejects = buildDM(rejects_currentMovie,T);


% sort by trim stage number
%for stage = 1:length(rejects_currentMovie)
    
%    dm_currentRejects = dm_currentMovie_rejects(dm_currentMovie_rejects(:,31) == stage,:); %not in D2
%    lostTracks{stage} = unique(dm_currentRejects(:,1));

%end
%clear stage;

%%

% save images of all tracks snipped  during specific stage of dataTrimmer

n=52;
stage = 2; % jumpTrack stage is stage 3 in original, 2 in scrambled

% define suplot parameters
total_subplots = 10;
cols = 2;
rz = total_subplots/cols;

numTracks = length(rejectD_scram{stage,n});
img = ones(ceil(numTracks/total_subplots),2);

for row = 1:length(img)
    % for 20 subplots per figure
    img(row,1) = img(row,1)+total_subplots*(row-1);
    img(row,2) = total_subplots*row;
    
end

img(length(img),2) = numTracks;

%%
for row = 1:length(img)
    %cla
    counter = 0;
    for i = img(row,1):img(row,2)%length(D6{n})
       
        counter = counter + 1;
        filename = strcat('checkingRejects-xy52-scrambled-edited-m',num2str(img(row,1)),'-',num2str(img(row,2)),'.tif');
        
        % designate subplot position
        %subplot(ceil(length(D6{n})/5), 5, i)
        subplot(rz, cols, counter)
        
        % plot
        %figure(counter)
        plot(T{n}(rejectD_scram{stage,n}(i).Frame(1:end))/3600,(rejectD_scram{stage,n}(i).MajAx),'Linewidth',2)
        
        % label
        title(i);
        xlim([0 11])
        
    end
    saveas(gcf,filename)
end