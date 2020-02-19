%% revealBirthVolumes

% Goal: Plot mean birth volumes over time

%       Like revealBirthSizes.m, this script looks for an evolution of cell cycle behavior.
%       Unlike it, this calculates an average of the volume (vs length or width) at birth across time.
   


%  Last edit: Jen Nguyen, 2017 Oct 27


% Strategy:
%
%      0.  initialize data and binning parameters
%      1.  specify current condition of interest
%               2.  isolate all data from current condition
%               3.  isolate volume, drop and timestamp data from current condition
%               4.  keep data from rows where isDrop = 1,
%                   constraining analysis to birth events
%               5.  bin volumes at birth by time of birth
%               6.  calculate average and s.e.m. of Vo per timebin
%               7.  convert bin # to absolute time (for plotting)
%               8.  plot time evolution of mean Vo
%               9.  trim data to limit analysis for time frames of stabilized growth
%              10.  bin birth volumes by size
%              11.  normalize bin quantities by total births
%              12.  plot pdf of birth volume
%              13.  plot pdf of all volumes (a check)
%     14.  repeat for all conditions


% OK! Lez go!

%% Initialize.   

% 0. initialze data
clc
clear

% trimmed dataset
load('lb-monod-2017-09-26-window5-jiggle-c12-0p1-c3456-0p5-bigger1p8.mat','D5','M','T');
load('lb-monod-2017-09-26-window5-va-jiggle-c12-0p1-c3456-0p5-bigger1p8.mat','M_va');
%load('lb-fluc-2017-10-10-window5-width1p4v1p7-jiggle-0p5-bigger1p8.mat','D5','M','T');
%load('lb-fluc-2017-10-10-window5va-width1p4v1p7-jiggle-0p5-bigger1p8.mat','M_va');
dataMatrix = buildDM(D5,M,M_va,T);
clear D5 M M_va T;

load('meta.mat');
meta = meta_2017sep26;

% 0. initialize binning parameters
expHours = 10;          % duration of experiment in hours                      
binFactor = 4;         % bins per hour
hrPerBin = 1/binFactor; % hour fraction per bin

%%
% 1.  specify current condition of interest
totalCond = max(dataMatrix(:,35)); % col 35 = condition value

for condition = 1:totalCond
    
    % 2. isolate all data from current condition
    interestingData = dataMatrix(dataMatrix(:,35) == condition,:);
    
    % 3. isolate volume, drop and timestamp data from current condition
    allVa = interestingData(:,15); % col 15 = calcalated va_vals (cubic um)
    timestamps = interestingData(:,2)/3600; % time in seconds converted to hours
    isDrop = interestingData(:,5);

    % 4. keep data from rows where isDrop = 1, constraining analysis to birth events
    birthTimes = timestamps(isDrop == 1);
    birthVa = allVa(isDrop == 1);
    
    % 5. bin volumes at birth by time of birth
    timeBins = ceil(birthTimes*binFactor);                
    binnedVa = accumarray(timeBins,birthVa,[],@(x) {x});
    
    % 6.  calculate average and s.e.m. of Vo per timebin
    meanVa = cellfun(@mean,binnedVa);
    countVa = cellfun(@length,binnedVa);
    stdVa = cellfun(@std,binnedVa);
    semVa = stdVa./sqrt(countVa);
    
    % 7.  convert bin # to absolute time
    timeVector = linspace(1, max(timeBins), max(timeBins));
    timeVector = hrPerBin*timeVector';
    
    % 8. plot time evolution of mean Vo 
    figure(1)
    errorbar(timeVector,meanVa,semVa)
    axis([0,10.5,0,12])
    hold on
    xlabel('Time (hr)')
    ylabel('Volume at birth + s.e.m. (cubic um)')
    legend('full LB','1/8 LB','1/32 LB','1/100 LB','1/1000 LB','1/10000 LB');
    
    
    % 9. trim data to limit analysis for time frames of stabilized growth
    minTime = meta(condition,3);  % hr
    maxTime = meta(condition,4);
    
    birthVa_trim1 = birthVa(birthTimes >= minTime);
    birthTimes_trim1 = birthTimes(birthTimes >= minTime);
    
    birthVa_trim2 = birthVa_trim1(birthTimes_trim1 <= maxTime);
    birthTimes_trim2 = birthTimes_trim1(birthTimes_trim1 <= maxTime);
    
    % 10. bin birth volumes by size
    binStable_Va = ceil(birthVa_trim2*10);
    binnedVa = accumarray(binStable_Va,birthVa_trim2,[],@(x) {x});
    binCounts_Va = cellfun(@length,binnedVa);
    
    % 11. normalize bin quantities by total births 
    stableVa_counts = length(birthVa_trim2);
    normalizedVa = binCounts_Va/stableVa_counts;
    
    
    % 12. plot pdf of birth volume 
    figure(2)
    subplot(totalCond,1,condition)
    bar(normalizedVa,0.4)
    axis([0,200,0,0.15])
    hold on
    xlabel('volume at birth (cubic um)')
    ylabel('pdf')
    legend(num2str(condition));
    
    % 13. plot pdf of all volumes
    vaBins = ceil(allVa*10);
    binned_va_all = accumarray(vaBins,allVa,[],@(x) {x});
    binCounts_va_all = cellfun(@length,binned_va_all);
    totalCount = length(allVa);
    normalizedCounts_all = binCounts_va_all./totalCount;
    
    figure(3)
    subplot(totalCond,1,condition)
    bar(normalizedCounts_all,0.4)
    xlabel('all volumes (cubic um)')
    ylabel('pdf')
    axis([0,200,0,0.1])
    
    
    
    % 14. repeat for all conditions
end

               






