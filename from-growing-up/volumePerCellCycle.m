%% volumePerCellCyle

% Goal: Plot the ratio of birth volume per cell cycle over time.

% Note: slower growing conditions are biased, as this script only measures
%       data at birth for complete cell cycles. Individual growth curves
%       that never end (drop) are excluded from this analysis.



%  Last edit: Jen Nguyen, 2017 Oct 22


% Strategy:
%
%      0.  initialize data and binning parameters
%      1.  specify current condition of interest
%               2.  isolate all data from current condition
%               3.  isolate volume, cell cycle duration, drop (birth event) and time data
%               4.  remove inappropriate duration data:
%                       i. from incomplete curves, where duration == 0
%                      ii. from "curves" that are too short to be physiological
%               5.  select rows where isDrop = 1, isolating data from birth events
%               6.  calculate birth volume per cell cycle
%               7.  bin birth volume/tau values by time
%               8.  calculate average and s.e.m. per timebin
%               9.  create a vector that converts timebin value to real time 
%              10.  plot birth volume/tau over time for all timepoints
%              11.  isolate data to stabilized regions of growth
%              12.  bin Vnot/tau production rates by timestamp
%              13.  normalize bin quantities by total birth events
%              14.  plot distribution of birth volumes/tau, post-stabilization
%              15.  plot birth volumes vs tau, using only post-stabilization values
%     16.  repeat for all conditions


% OK! Lez go!

%%
%   Initialize.

% 0. initialze data
clc
clear

% trimmed dataset
load('lb-fluc-2017-10-10-window5-width1p4v1p7-jiggle-0p5-bigger1p8.mat','D5','M','T');
load('lb-fluc-2017-10-10-window5va-width1p4v1p7-jiggle-0p5-bigger1p8.mat','M_va');
dataMatrix = buildDM(D5,M,M_va,T);
clear D5 M M_va T;

load('meta.mat');
meta = meta_2017oct10;

% 0. initialize binning parameters
expHours = 10;          % duration of experiment in hours                      
binFactor = 4;          % bins per hour
hrPerBin = 1/binFactor; % hour fraction per bin

% 0. define "too short to be physiological" in terms of curve duration
tooRapid = 11;          % in min

%%
% 1.  specify current condition of interest
totalCond = max(dataMatrix(:,35)); % col 35 = condition value

for condition = 1:totalCond
    
    % 2. isolate all data from current condition
    interestingData = dataMatrix(dataMatrix(:,35) == condition,:);
    
    % 3. isolate volume, cell cycle duration, drop (birth event) and time data
    va_vals = interestingData(:,15);        % col 15 = calculated va_vals (cubic um)
    durations = interestingData(:,8)/60;    % col 8 = curve (cell cycle) duration in sec converted to min
    timestamps = interestingData(:,2)/3600; % time in seconds converted to hours
    isDrop = interestingData(:,5);          % col 5 = isDrop
    
    % 4. remove inappropriate duration data:
    %    i. from incomplete curves, where duration == 0
    completeDurations = durations(durations > 0);
    completeVas = va_vals(durations > 0);
    completeTimes = timestamps(durations > 0);
    completeDrops = isDrop(durations > 0);  % now all drop==1 correspond to births for full curves
    
    %   ii. from "curves" that are too short to be physiological
    taus = completeDurations(completeDurations >= tooRapid);
    Va_nots = completeVas(completeDurations >= tooRapid);
    times = completeTimes(completeDurations >= tooRapid);
    drops = completeDrops(completeDurations >= tooRapid);
    
    % 5. select rows where isDrop = 1, isolating data from birth events
    birthTimes = times(drops == 1);
    birthVa = Va_nots(drops == 1);
    birthDurations = taus(drops ==1);
    
    % 6. calculate birth volume per cell cycle
    Vnot_over_tau = birthVa./birthDurations;
    
    % 7. bin birth volume/tau values by time
    timeBins = ceil(birthTimes*binFactor);
    binned = accumarray(timeBins,Vnot_over_tau,[],@(x) {x});
    
    % 8. calculate average and s.e.m. per timebin
    mean_VoPerCC = cellfun(@mean,binned);
    count_VoPerCC = cellfun(@length,binned);
    std_VoPerCC = cellfun(@std,binned);
    sem_VoPerCC = std_VoPerCC./sqrt(count_VoPerCC);
    
    % 9. create a vector that converts timebin value to real time
    rtVector = linspace(1, max(timeBins), max(timeBins));
    rtVector = hrPerBin*rtVector'; 
    
    % 10. plot birth volume/tau over time for all timepoints 
    figure(1)
    errorbar(rtVector,mean_VoPerCC,sem_VoPerCC)
    axis([0,10.5,0,0.3])
    hold on
    xlabel('Time (hr)')
    ylabel('Vo per tau (cubic um/hr)')
    legend('fluc','1/1000 LB','ave','1/50 LB');
    
    
    % 11.  isolate data to stabilized regions of growth
    minTime = meta(condition,3);  % hr
    maxTime = meta(condition,4);
    
    % for figure(2): histogram
    Vnot_over_tau_trim1 = Vnot_over_tau(birthTimes >= minTime);
    birthTimes_trim1 = birthTimes(birthTimes >= minTime);
    
    Vnot_over_tau_trim2 = Vnot_over_tau_trim1(birthTimes_trim1 <= maxTime);
    birthTimes_trim2 = birthTimes_trim1(birthTimes_trim1 <= maxTime);
    
    % for figure(3): Vo vs tau
    birthVa_trim1 = birthVa(birthTimes >= minTime);
    birthDurations_trim1 = birthDurations(birthTimes >= minTime);
    
    birthVa_trim2 = birthVa_trim1(birthTimes_trim1 <= maxTime);
    birthDurations_trim2 = birthDurations_trim1(birthTimes_trim1 <= maxTime);
    
    
    % 12.  bin Vnot/tau values by timestamp
    binStable_VoPerCC = ceil(Vnot_over_tau_trim2*100);
    binnedTrimmed = accumarray(binStable_VoPerCC,Vnot_over_tau_trim2,[],@(x) {x});
    binCounts_VoPerCC = cellfun(@length,binnedTrimmed);
    
    % 13. normalize bin quantities by total birth events
    countStable = length(birthTimes_trim2);
    pdf_VnotPerCC = binCounts_VoPerCC/countStable;
    
    
    % 14. plot distribution of birth volumes/tau, post-stabilization
    figure(2)
    subplot(totalCond,1,condition)
    bar(pdf_VnotPerCC,0.4)
    axis([0,100,0,0.4])
    hold on
    xlabel('Vo per tau (cubic um/hr)')
    ylabel('pdf')
    legend(num2str(condition));
    
    
    % 15. plot birth volumes vs tau, using only post-stabilization values
    figure(3)
    plot(birthDurations_trim2,birthVa_trim2,'o')
    hold on
    xlabel('Length of cell cycle (min)')
    ylabel('Volume at birth (cubic um)')
    legend('fluc','1/1000 LB','ave','1/50 LB');
   
    % 9. repeat for all conditions
end

               






