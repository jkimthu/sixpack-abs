%% revealBirthSize

% Goal: Plot mean birth size over time

%       Like revealDurations.m, this script looks for an evolution of cell cycle behavior.
%       Unlike it, this calculates an average of the size at birth across time.
   


%  Last edit: Jen Nguyen, 2017 Oct 13



% As of Sept 2017, this script uses function, buildDM, instead of specific
% matrixBuilder.m data outputs. The overall concept is still similar.



% Strategy:
%
%      0.  initialize data and binning parameters
%      1.  specify current condition of interest
%               2.  isolate data from current condition
%               3.  trim data to only account for stabilized growth
%               4.  accumulate size at birth data by timebin
%               5.  convert bin # to absolute time
%               6.  calculate average and s.e.m. per timebin
%               7.  plot!
%      8.  repeat for all conditions


% OK! Lez go!

%%
%   Initialize.

% 0. initialze data
clc
clear

% trimmed dataset
load('lb-fluc-2017-10-10-window5-width1p4v1p7-jiggle-0p5-bigger1p8.mat','D5','M','T');
dataMatrix = buildDM(D5,M,T);
load('meta.mat');
meta = meta_2017oct10;

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
    
    % i. isolate length and time data
    allLengths = interestingData(:,3); % col 3 = measured lengths (um)
    allWidths = interestingData(:,12); % col 12 = measured widths (um)
    timestamps = interestingData(:,2)/3600; % time in seconds converted to hours
    
    % ii. select rows where isDrop = 1
    isDrop = interestingData(:,5);
    birthTimes = timestamps(isDrop == 1);
    birthLengths = allLengths(isDrop == 1);
    birthWidths = allWidths(isDrop == 1);
    
    
    % 3.  trim data to only account for stabilized growth
    
    % i. remove data not in stabilized region
    minTime = meta(condition,3);  % hr
    maxTime = meta(condition,4);
    
    birthLengths_trim1 = birthLengths(birthTimes >= minTime);
    birthWidths_trim1 = birthWidths(birthTimes >= minTime);
    birthTimes_trim1 = birthTimes(birthTimes >= minTime);
    
    birthLengths_trim2 = birthLengths_trim1(birthTimes_trim1 <= maxTime);
    birthWidths_trim2 = birthWidths_trim1(birthTimes_trim1 <= maxTime);
    birthTimes_trim2 = birthTimes_trim1(birthTimes_trim1 <= maxTime);
    
    
    % 4.  accumulate size at birth data by timebin
    
    % iii. convert birthTimes into timebins
    timeBins = ceil(birthTimes_trim2*binFactor);                
    binnedSizes = accumarray(timeBins,birthLengths_trim2,[],@(x) {x});
    binnedWidths = accumarray(timeBins,birthWidths_trim2,[],@(x) {x});
    
    
    % 4.  convert bin # to absolute time
    timeVector = linspace(1, max(timeBins), max(timeBins));
    timeVector = hrPerBin*timeVector'; 
    
    
    % 5.  calculate average and s.e.m. per timebin
    meanVector = cellfun(@mean,binnedSizes);
    countVector = cellfun(@length,binnedSizes);
    stdVector = cellfun(@std,binnedSizes);
    semVector = stdVector./sqrt(countVector);
    
    meanWidth = cellfun(@mean,binnedWidths);
    countWidth = cellfun(@length,binnedWidths);
    stdWidth = cellfun(@std,binnedWidths);
    semWidth = stdWidth./sqrt(countWidth);
    
    
    % 6.  plot!
    figure(1)
    errorbar(timeVector,meanVector,semVector)
    axis([0,10.3,1,8])
    hold on
    xlabel('Time (hr)')
    ylabel('Length at birth + s.e.m. (min)')
    legend('fluc','1/1000 LB','ave','1/50 LB'); 
    
%     figure(2)
%     errorbar(timeVector,meanVector,stdVector)
%     axis([0,10.3,0,10])
%     hold on
%     xlabel('Time (hr)')
%     ylabel('Length at birth + standard dev (min)')
%     legend('fluc','1/1000 LB','ave','1/50 LB'); 
    
    % 7. plot pdfs from steady-state
    stableBirthLengths = birthLengths(birthTimes > 3);
    
    figure(3)
    histogram(stableBirthLengths,'Normalization','pdf')
    axis([0,12,0,1.5])
    hold on
    xlabel('Length at birth (um)')
    ylabel('pdf')
    legend('fluc','1/1000 LB','ave','1/50 LB'); 
    
    figure(4)
    subplot(totalCond,1,condition)
    histogram(stableBirthLengths,'Normalization','pdf')
    axis([0,20,0,1.5])
    hold on
    xlabel('Length at birth (um)')
    ylabel('pdf')
    legend(num2str(condition));
    
    figure(5)
    errorbar(timeVector,meanWidth,semWidth)
    axis([0,10.3,0.75,2.25])
    hold on
    xlabel('Time (hr)')
    ylabel('Width at birth + s.e.m. (min)')
    legend('fluc','1/1000 LB','ave','1/50 LB');
    
%     figure(6)
%     errorbar(timeVector,meanWidth,stdWidth)
%     axis([0,10.3,0.75,2.25])
%     hold on
%     xlabel('Time (hr)')
%     ylabel('Width at birth + standard dev (min)')
%     legend('fluc','1/1000 LB','ave','1/50 LB');
    
    % 7. plot pdfs from steady-state
    stableBirthWidths = birthWidths(birthTimes > 3);
    
    figure(7)
    histogram(stableBirthWidths,'Normalization','pdf')
    axis([0.8,1.8,0,8])
    hold on
    xlabel('Width at birth (um)')
    ylabel('pdf')
    legend('fluc','1/1000 LB','ave','1/50 LB');
    
    figure(8)
    subplot(totalCond,1,condition)
    histogram(stableBirthWidths,'Normalization','pdf')
    axis([0.8,1.8,0,9])
    hold on
    xlabel('Width at birth (um)')
    ylabel('pdf')
    legend(num2str(condition));
    
    % 8. repeat for all conditions
end

               






