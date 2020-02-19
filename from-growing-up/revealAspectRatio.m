%% revealAspectRatio

% Goal: Plot mean birth aspect ratios over time

%       Like revealBirthSizes.m, this script looks for an evolution of cell cycle behavior.
%       Unlike it, this calculates the average aspect ratio (width:length) at birth across time.
   


%  Last edit: Jen Nguyen, 2017 Oct 11



% Strategy:
%
%      0.  initialize data and binning parameters
%      1.  specify current condition of interest
%               2.  isolate length and width data from current condition
%               3.  accumulate aspect ratios at birth data by timebin
%               4.  convert bin # to absolute time
%               5.  calculate average and s.e.m. per timebin
%               6.  plot!
%      7.  repeat for all conditions


% OK! Lez go!

%%
%   Initialize.

% 0. initialze data
clc
clear

% trimmed dataset
load('lb-fluc-2017-10-10-window5-jiggle-0p5-bigger1p8.mat','D5','M','T');
dataMatrix = buildDM(D5,M,T);

% 0. initialize binning parameters
expHours = 10;          % duration of experiment in hours                      
binFactor = 6;         % bins per hour
hrPerBin = 1/binFactor; % hour fraction per bin

%%
% 1.  specify current condition of interest
totalCond = max(dataMatrix(:,35)); % col 35 = condition value

for condition = 1:totalCond
    
    % 2.  isolate data from current condition
    interestingData = dataMatrix(dataMatrix(:,35) == condition,:);
    
    % 3.  accumulate aspect ratio at birth data by timebin
    
    % i. isolate length, width and time data
    allLengths = interestingData(:,3); % col 3 = measured lengths (um)
    allWidths = interestingData(:,12);     % col 12 = measured widths (um)
    timestamps = interestingData(:,2)/3600; % time in seconds converted to hours
    
    % ii. select rows where isDrop = 1
    isDrop = interestingData(:,5);
    birthTimes = timestamps(isDrop == 1);
    birthLengths = allLengths(isDrop == 1);
    birthWidths = allWidths(isDrop == 1);
    
    % iii. calculate aspect ratios
    aspectRatios = birthWidths./birthLengths;
    
    % iv. convert birthTimes into timebins
    timeBins = ceil(birthTimes*binFactor);                
    binnedRatios = accumarray(timeBins,aspectRatios,[],@(x) {x});

    
    % 4.  convert bin # to absolute time
    timeVector = linspace(1, max(timeBins), max(timeBins));
    timeVector = hrPerBin*timeVector'; 
    
    
    % 5.  calculate average and s.e.m. per timebin
    meanRatio = cellfun(@mean,binnedRatios);
    countRatio = cellfun(@length,binnedRatios);
    stdRatio = cellfun(@std,binnedRatios);
    semRatio = stdRatio./sqrt(countRatio);
    
    
    % 6.  plot 
    figure(1)
    errorbar(timeVector,meanRatio,semRatio)
    axis([0,10.5,0,0.7])
    hold on
    xlabel('Time (hr)')
    ylabel('Aspect atio at birth + s.e.m.')
    legend('fluc','1/1000 LB','ave','1/50 LB');
    
    
    
    % 7. plot pdfs from steady-state
    
    % i. isolate data from stabilized timepoints
    stableAspectRatios = aspectRatios(birthTimes > 3);
    
    % ii. bin birth volumes
    bins = ceil(stableAspectRatios*100);
    binnedRatios = accumarray(bins,stableAspectRatios,[],@(x) {x});
    binCounts = cellfun(@length,binnedRatios);
    
    % iii. normalize bin quantities by total births 
    stable_counts = length(stableAspectRatios);
    normalizedRatios = binCounts/stable_counts;
   
    
    figure(2)
    subplot(totalCond,1,condition)
    bar(normalizedRatios,0.4)
    axis([0,100,0,0.12])
    hold on
    xlabel('aspect ratio at birth (um)')
    ylabel('population fraction')
    legend(num2str(condition));
    
   
    % 8. repeat for all conditions
end

               






