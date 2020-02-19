%% visualizeTrimmer.m


% Goal: plot tracks of surviving and trimmed tracks
%


clear
clc
experiment = '2017-01-16';


% TRACKING DATA
% open folder for experiment of interest
newFolder = strcat('/Users/jen/Documents/StockerLab/Data/',experiment,'  (t300)');
cd(newFolder);

% FROM DATA TRIMMER
% particle tracking data
clear
load('t300_2017-01-16-revisedTrimmer-jiggle0p3.mat','D7','D','T','rejectD');

%%
n = 1;
criteria = 3;

%subset = rejectD(:,n);
subset = D5{n}(1:100,:);
%subset = D3{n}(test,:);

numTracks = length(subset);%{criteria});

%%
% initialize figure count
numSubplots = 20;
figureBounds = [1 numSubplots];

numFigs = ceil(numTracks/numSubplots);
for sp = 1:numFigs-1
    figureBounds = [figureBounds; figureBounds(end,:) + numSubplots];
end


%%
for f = 1:numFigs
    
    cla
    subplot_counter = 0;
    %filename = strcat('dynamicOutlines-subplotGroup',num2str(f),'-n',num2str(n),'.png');
    
    for i = figureBounds(f,1):figureBounds(f,2)
        
        subplot_counter = subplot_counter + 1;
        % designate subplot position
        subplot(4, 5, subplot_counter)
        
        % plot
        figure(f)
        %plot(subset{criteria,1}(i).MajAx,'Linewidth',2)
        plot(subset(i).MajAx,'Linewidth',2)
        
        % label
        %title(subset{criteria,1}(i).TrackID(1));
        title(subset(i).TrackID(1));
        
    end
    
    %saveas(gcf,filename)
    
end 