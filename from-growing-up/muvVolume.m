% muvVolume.m

% goal: plot mean volume at birth against mean growth rate (stabilized)

%  Last edit: Jen Nguyen, 2017 Oct 25



% Strategy:
%
%      0.  initialize data and binning parameters
%      1.  specify current condition of interest
%               2.  isolate all data from current condition
%               3.  isolate volume, mu_va, drop and timestamp data
%               4.  keep volume data from rows where isDrop = 1,
%                   constraining analysis to birth events
%               5.  keep only non-zero values of mu
%               6.  bin volumes at birth by time of birth
%               7.  bin mus by time (full)
%               8.  calculate average and s.e.m. of Vo and mu per timebin
%               9.  convert bin # to absolute time (for plotting)
%              10.  plot time evolution of mean Vo and mean mu
%              11.  trim data to limit analysis for time frames of stabilized growth
%              12.  find the mean and s.e.m. of birth volume, during stabilized periods
%              13.  find the mean and s.e.m. of mu_va, during stable periods
%              14.  plot mean birth volume vs mean mu
%     15.  repeat for all conditions

%% Initialize 2017-09-26 data

% 0. initialze data
clc
clear
experiment = '2017-09-26';

% 0. open folder for experiment of interest
newFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',experiment);%,'  (t300)');
cd(newFolder);

% trimmed dataset
load('lb-monod-2017-09-26-window5-jiggle-c12-0p1-c3456-0p5-bigger1p8.mat','D5','M','T');
load('lb-monod-2017-09-26-window5-va-jiggle-c12-0p1-c3456-0p5-bigger1p8.mat','M_va');
dataMatrix = buildDM(D5,M,M_va,T);
clear D5 M M_va T;

load('meta.mat');
meta = meta_2017sep26;

% 0. initialize binning parameters
expHours = 10;          % duration of experiment in hours                      
binFactor = 4;         % bins per hour
hrPerBin = 1/binFactor; % hour fraction per bin


%% Plot 2017-09-26 data

% 1.  specify current condition of interest
totalCond = max(dataMatrix(:,35)); % col 35 = condition value

for condition = 1:totalCond
    
    % 2. isolate all data from current condition
    interestingData = dataMatrix(dataMatrix(:,35) == condition,:);
    
    % 3. isolate volume, drop and timestamp data from current condition
    allVa = interestingData(:,15); % col 15 = calcalated va_vals (cubic um)
    mus_va = interestingData(:,18); % col 18 = calculated mu_vals
    timestamps = interestingData(:,2)/3600; % time in seconds converted to hours
    isDrop = interestingData(:,5);

    % 4. keep volume data from rows where isDrop = 1, constraining analysis to birth events
    birthTimes = timestamps(isDrop == 1);
    birthVa = allVa(isDrop == 1);
    
    % 5. keep only non-zero values of mu
    trueMus = mus_va(mus_va > 0);
    trueTimes = timestamps(mus_va > 0);
    
    % 6. bin volumes at birth by time of birth
    timeBins = ceil(birthTimes*binFactor);                
    binnedVa = accumarray(timeBins,birthVa,[],@(x) {x});
    
    % 7. bin mus by time (full)
    timeBins_mu = ceil(trueTimes*binFactor);
    binnedMus = accumarray(timeBins_mu,trueMus,[],@(x) {x});
    
    % 8.  calculate average and s.e.m. of Vo and mu per timebin
    meanVa = cellfun(@mean,binnedVa);
    countVa = cellfun(@length,binnedVa);
    stdVa = cellfun(@std,binnedVa);
    semVa = stdVa./sqrt(countVa);
    
    meanMu = cellfun(@mean,binnedMus);
    countMu = cellfun(@length,binnedMus);
    stdMu = cellfun(@std,binnedMus);
    semMu = stdMu./sqrt(countMu);
    
    % 9.  convert bin # to absolute time
    timeVector = linspace(1, max(timeBins), max(timeBins));
    timeVector = hrPerBin*timeVector';
    
    timeVector_mu = linspace(1, max(timeBins_mu), max(timeBins_mu));
    timeVector_mu = hrPerBin*timeVector_mu';
    
    
    % 10. plot time evolution of mean Vo and mean mu
    figure(1)
    errorbar(timeVector,meanVa,semVa,'--')
    axis([0,10.5,0,12])
    hold on
    xlabel('Time (hr)')
    ylabel('Volume at birth + s.e.m. (cubic um)')
    legend('full LB','1/8 LB','1/32 LB','1/100 LB','1/1000 LB','1/10000 LB');
    
    figure(2)
    errorbar(timeVector_mu,meanMu,semMu,'-')
    axis([0,10.5,0,12])
    hold on
    xlabel('Time (hr)')
    ylabel('Mu + s.e.m. (1/hr)')
    legend('full LB','1/8 LB','1/32 LB','1/100 LB','1/1000 LB','1/10000 LB');
    
    
    % 11. trim data to limit analysis for time frames of stabilized growth
    minTime = meta(condition,3);  % hr
    maxTime = meta(condition,4);
    
    %    i. Vo
    birthVa_trim1 = birthVa(birthTimes >= minTime);
    birthTimes_trim1 = birthTimes(birthTimes >= minTime);
    
    birthVa_trim2 = birthVa_trim1(birthTimes_trim1 <= maxTime);
    birthTimes_trim2 = birthTimes_trim1(birthTimes_trim1 <= maxTime);
    
    %   ii. mu
    trueMus_trim1 = trueMus(trueTimes >= minTime);
    trueTimes_trim1 = trueMus(trueTimes >= minTime);
    
    trueMus_trim2 = trueMus_trim1(trueTimes_trim1 <= maxTime);
    trueTimes_trim2 = trueTimes_trim1(trueTimes_trim1 <= maxTime);
    
    % 12. find the mean and s.e.m. of birth volume, during stabilized periods
    meanVa_stable = mean(birthVa_trim2);
    countVa_stable = length(birthVa_trim2);
    stdVa_stable = std(birthVa_trim2);
    semVa_stable = stdVa_stable./sqrt(countVa_stable);
    
    % 13. find the mean and s.e.m. of mu_va, during stable periods
    meanMu_stable = mean(trueMus_trim2);
    countMu_stable = length(trueMus_trim2);
    stdMu_stable = std(trueMus_trim2);
    semMu_stable = stdMu_stable./sqrt(countMu_stable);
    
    
    % 14. plot mean birth volume vs mean mu
    figure(3)
    % errorbar(x,y,yneg,ypos,xneg,xpos,'o')
    errorbar(meanMu_stable,meanVa_stable,semVa_stable,semVa_stable,'o','MarkerSize',10)
    axis([0,4,1,9])
    hold on
    xlabel('mu (1/hr)')
    ylabel('volume at birth (cubic um)')
    legend('full LB','1/8 LB','1/32 LB','1/100 LB','1/1000 LB','1/10000 LB');
    
    % 13. repeat for all conditions
end

%%  Initialize 2017-10-10 data

% 0. initialze data
clc
clear
experiment = '2017-10-10';

% 0. open folder for experiment of interest
newFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',experiment);%,'  (t300)');
cd(newFolder);

% trimmed dataset
load('lb-fluc-2017-10-10-window5-width1p4v1p7-jiggle-0p5-bigger1p8.mat','D5','M','T');
load('lb-fluc-2017-10-10-window5va-width1p4v1p7-jiggle-0p5-bigger1p8.mat','M_va');
dataMatrix = buildDM(D5,M,M_va,T);
clear D5 M M_va T;

load('meta.mat');
meta = meta_2017oct10;

% 0. initialize binning parameters
expHours = 10;          % duration of experiment in hours                      
binFactor = 4;         % bins per hour
hrPerBin = 1/binFactor; % hour fraction per bin

%% Plot 2017-10-10

% 1.  specify current condition of interest
totalCond = max(dataMatrix(:,35)); % col 35 = condition value

for condition = 1:totalCond
    
    % 2. isolate all data from current condition
    interestingData = dataMatrix(dataMatrix(:,35) == condition,:);
    
    % 3. isolate volume, drop and timestamp data from current condition
    allVa = interestingData(:,15); % col 15 = calcalated va_vals (cubic um)
    mus_va = interestingData(:,18); % col 18 = calculated mu_vals
    timestamps = interestingData(:,2)/3600; % time in seconds converted to hours
    isDrop = interestingData(:,5);

    % 4. keep volume data from rows where isDrop = 1, constraining analysis to birth events
    birthTimes = timestamps(isDrop == 1);
    birthVa = allVa(isDrop == 1);
    
    % 5. keep only non-zero values of mu
    trueMus = mus_va(mus_va > 0);
    trueTimes = timestamps(mus_va > 0);
    
    % 6. bin volumes at birth by time of birth
    timeBins = ceil(birthTimes*binFactor);                
    binnedVa = accumarray(timeBins,birthVa,[],@(x) {x});
    
    % 7. bin mus by time (full)
    timeBins_mu = ceil(trueTimes*binFactor);
    binnedMus = accumarray(timeBins_mu,trueMus,[],@(x) {x});
    
    % 8.  calculate average and s.e.m. of Vo and mu per timebin
    meanVa = cellfun(@mean,binnedVa);
    countVa = cellfun(@length,binnedVa);
    stdVa = cellfun(@std,binnedVa);
    semVa = stdVa./sqrt(countVa);
    
    meanMu = cellfun(@mean,binnedMus);
    countMu = cellfun(@length,binnedMus);
    stdMu = cellfun(@std,binnedMus);
    semMu = stdMu./sqrt(countMu);
    
    % 9.  convert bin # to absolute time
    timeVector = linspace(1, max(timeBins), max(timeBins));
    timeVector = hrPerBin*timeVector';
    
    timeVector_mu = linspace(1, max(timeBins_mu), max(timeBins_mu));
    timeVector_mu = hrPerBin*timeVector_mu';
    
    
    % 10. plot time evolution of mean Vo and mean mu
    figure(1)
    errorbar(timeVector,meanVa,semVa,'--')
    axis([0,10.5,0,12])
    hold on
    xlabel('Time (hr)')
    ylabel('Volume at birth + s.e.m. (cubic um)')
    legend('fluc','1/1000 LB','ave','1/50 LB');
    
    figure(2)
    errorbar(timeVector_mu,meanMu,semMu,'-')
    axis([0,10.5,0,12])
    hold on
    xlabel('Time (hr)')
    ylabel('Mu + s.e.m. (1/hr)')
    legend('fluc','1/1000 LB','ave','1/50 LB');
    
    
    % 11. trim data to limit analysis for time frames of stabilized growth
    minTime = meta(condition,3);  % hr
    maxTime = meta(condition,4);
    
    %    i. Vo
    birthVa_trim1 = birthVa(birthTimes >= minTime);
    birthTimes_trim1 = birthTimes(birthTimes >= minTime);
    
    birthVa_trim2 = birthVa_trim1(birthTimes_trim1 <= maxTime);
    birthTimes_trim2 = birthTimes_trim1(birthTimes_trim1 <= maxTime);
    
    %   ii. mu
    trueMus_trim1 = trueMus(trueTimes >= minTime);
    trueTimes_trim1 = trueMus(trueTimes >= minTime);
    
    trueMus_trim2 = trueMus_trim1(trueTimes_trim1 <= maxTime);
    trueTimes_trim2 = trueTimes_trim1(trueTimes_trim1 <= maxTime);
    
    % 12. find the mean and s.e.m. of birth volume, during stabilized periods
    meanVa_stable = mean(birthVa_trim2);
    countVa_stable = length(birthVa_trim2);
    stdVa_stable = std(birthVa_trim2);
    semVa_stable = stdVa_stable./sqrt(countVa_stable);
    
    % 13. find the mean and s.e.m. of mu_va, during stable periods
    meanMu_stable = mean(trueMus_trim2);
    countMu_stable = length(trueMus_trim2);
    stdMu_stable = std(trueMus_trim2);
    semMu_stable = stdMu_stable./sqrt(countMu_stable);
    
    
    % 14. plot mean birth volume vs mean mu
    figure(3)
    % errorbar(x,y,yneg,ypos,xneg,xpos,'o')
    errorbar(meanMu_stable,meanVa_stable,semVa_stable,semVa_stable,'o','MarkerSize',10)
    axis([0,4,1,9])
    hold on
    xlabel('mu (1/hr)')
    ylabel('volume at birth (cubic um)')
    legend('fluc','1/1000 LB','ave','1/50 LB');
    
    % 13. repeat for all conditions
end