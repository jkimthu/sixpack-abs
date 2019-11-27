%% distribute Normally

%  Goal: plot distributions of cell cycle duration and size at birth,
%  normalized by the average value in corresponding stable environment



%  Last edit: Jen Nguyen, August 14th 2016


%  Sections:
%       one. cell cycle duration
%       two. size at birth


%  Strategy:
%       0. initialize data
%       1. trim data points to times of interest (stabilized)
%       2. find mean of constant
%       3. normalize fluc data by constant mean
%       4. plot distributions of 15 min and 60 min envir on same histogram


% OK! Lez go!

%% ONE. CELL CYCLE DURATION

% The intended input for these scripts is a data structure, spinOffs,
% created using dataMatrix.m

% dF_MMDD.mat
% dC_MMDD.mat




%       0. initialize all data: durations and timestamps

% from 60 min period experiment, 2015-08-10
load('dC_0810.mat')
dC = spinOffs;
load('dF_0810.mat')
dF = spinOffs;
clear spinOffs;

p60_cDurations = dC.allDurations;
p60_fDurations = dF.allDurations;

p60_cTime = dC.allTimestamps;
p60_fTime = dF.allTimestamps;
clear spinOffs dC dF;

%%
% from 15 min period experiment, 2015-08-18
load('dC_0818.mat')
dC = spinOffs;
load('dF_0818.mat')
dF = spinOffs;

p15_cDurations = dC.allDurations;
p15_fDurations = dF.allDurations;

p15_cTime = dC.allTimestamps;
p15_fTime = dF.allTimestamps;
clear spinOffs dC dF;


%%
%       1. trim data points to times of interest (stabilized)

% designate time window, in hours
firstTimepoint = 5; 
lastTimepoint = 10;

% 1. trim timepoints earlier than first
p60_cDurations = p60_cDurations(p60_cTime >= firstTimepoint);
lowTrimmed_p60_cTime = p60_cTime(p60_cTime >= firstTimepoint);

p60_fDurations = p60_fDurations(p60_fTime >= firstTimepoint);
lowTrimmed_p60_fTime = p60_fTime(p60_fTime >= firstTimepoint);

p15_cDurations = p15_cDurations(p15_cTime >= firstTimepoint);
lowTrimmed_p15_cTime = p15_cTime(p15_cTime >= firstTimepoint);

p15_fDurations = p15_fDurations(p15_fTime >= firstTimepoint);
lowTrimmed_p15_fTime = p15_fTime(p15_fTime >= firstTimepoint);

%clear firstTimepoint p60_cTime p60_fTime p15_cTime p15_fTime;



% 1. trim off timepoints later than last
p60_cDurations = p60_cDurations(lowTrimmed_p60_cTime <= lastTimepoint);
finalTrimmed_p60_cTime = lowTrimmed_p60_cTime(lowTrimmed_p60_cTime <= lastTimepoint);

p60_fDurations = p60_fDurations(lowTrimmed_p60_fTime <= lastTimepoint);
finalTrimmed_p60_fTime = lowTrimmed_p60_fTime(lowTrimmed_p60_fTime <= lastTimepoint);

p15_cDurations = p15_cDurations(lowTrimmed_p15_cTime <= lastTimepoint);
finalTrimmed_p15_cTime = lowTrimmed_p15_cTime(lowTrimmed_p15_cTime <= lastTimepoint);

p15_fDurations = p15_fDurations(lowTrimmed_p15_fTime <= lastTimepoint);
finalTrimmed_p15_fTime = lowTrimmed_p15_fTime(lowTrimmed_p15_fTime <= lastTimepoint);

%clear lastTimepoint lowTrimmed_p60_cTime lowTrimmed_p60_fTime lowTrimmed_p15_cTime lowTrimmed_p15_fTime;

    

% 2. calculate mean of constant environment
mean_p60_cDur = mean(p60_cDurations);
mean_p15_cDur = mean(p15_cDurations);

median_p60_cDur = median(p60_cDurations);
median_p15_cDur = median(p15_cDurations);
    
%clear p60_cDurations p15_cDurations;
    


% 3. normalize fluc data by constant mean
norm_p60_fDur = p60_fDurations./mean_p60_cDur;
norm_p15_fDur = p15_fDurations./mean_p15_cDur;



% 4. normalize fluc data by constant median
median_norm_p60_fDur = p60_fDurations./median_p60_cDur;
median_norm_p15_fDur = p15_fDurations./median_p15_cDur;

%clear p60_fDurations p15_fDurations;



% 5. plot distributions of 15 min and 60 min envir on same histogram

figure(1)
histogram(norm_p60_fDur,'Normalization', 'probability', 'BinWidth',0.1)
hold on
histogram(norm_p15_fDur,'Normalization', 'probability', 'BinWidth',0.1)
xlabel('cell cycle duration, normalized by constant mean')
legend('60 min period','15 min period')

%figure(2)
%histogram(median_norm_p60_fDur,'Normalization', 'probability', 'BinWidth',0.1)
%hold on
%histogram(median_norm_p15_fDur,'Normalization', 'probability', 'BinWidth',0.1)
%xlabel('cell cycle duration, normalized by constant median')
%legend('60 min period','15 min period')

figure(3)
histogram(p60_fDurations,'Normalization', 'probability', 'BinWidth', 0.1)
hold on
histogram(p15_fDurations,'Normalization', 'probability', 'BinWidth', 0.1)
hold on
histogram(p60_cDurations,'Normalization', 'probability', 'BinWidth', 0.1)
hold on
histogram(p15_cDurations,'Normalization', 'probability', 'BinWidth', 0.1)
xlabel('cell cycle duration, without normalization')
legend('60 min, fluc', '15 min, fluc', '60 min, const', '15 min, fluc')


% yay!



%% TWO. SIZE AT BIRTH



%       0. initialize all data: durations and timestamps

% from 60 min period experiment, 2015-08-10
load('dC_0810.mat')
dC = spinOffs;
load('dF_0810.mat')
dF = spinOffs;
clear spinOffs;

p60_cSizes = dC.birthSizes;
p60_fSizes = dF.birthSizes;

p60_cTime = dC.birthTimes;
p60_fTime = dF.birthTimes;
clear spinOffs dC dF;

%%
% from 15 min period experiment, 2015-08-18
load('dC_0818.mat')
dC = spinOffs;
load('dF_0818.mat')
dF = spinOffs;

p15_cSizes = dC.birthSizes;
p15_fSizes = dF.birthSizes;

p15_cTime = dC.birthTimes;
p15_fTime = dF.birthTimes;
clear spinOffs dC dF;


%%
%       1. trim data points to times of interest (stabilized)

% designate time window, in hours
firstTimepoint = 5; 
lastTimepoint = 10;

% 1. trim timepoints earlier than first
p60_cSizes = p60_cSizes(p60_cTime >= firstTimepoint);
lowTrimmed_p60_cTime = p60_cTime(p60_cTime >= firstTimepoint);

p60_fSizes = p60_fSizes(p60_fTime >= firstTimepoint);
lowTrimmed_p60_fTime = p60_fTime(p60_fTime >= firstTimepoint);

p15_cSizes = p15_cSizes(p15_cTime >= firstTimepoint);
lowTrimmed_p15_cTime = p15_cTime(p15_cTime >= firstTimepoint);

p15_fSizes = p15_fSizes(p15_fTime >= firstTimepoint);
lowTrimmed_p15_fTime = p15_fTime(p15_fTime >= firstTimepoint);

%clear firstTimepoint p60_cTime p60_fTime p15_cTime p15_fTime;



% 1. trim off timepoints later than last
p60_cSizes = p60_cSizes(lowTrimmed_p60_cTime <= lastTimepoint);
finalTrimmed_p60_cTime = lowTrimmed_p60_cTime(lowTrimmed_p60_cTime <= lastTimepoint);

p60_fSizes = p60_fSizes(lowTrimmed_p60_fTime <= lastTimepoint);
finalTrimmed_p60_fTime = lowTrimmed_p60_fTime(lowTrimmed_p60_fTime <= lastTimepoint);

p15_cSizes = p15_cSizes(lowTrimmed_p15_cTime <= lastTimepoint);
finalTrimmed_p15_cTime = lowTrimmed_p15_cTime(lowTrimmed_p15_cTime <= lastTimepoint);

p15_fSizes = p15_fSizes(lowTrimmed_p15_fTime <= lastTimepoint);
finalTrimmed_p15_fTime = lowTrimmed_p15_fTime(lowTrimmed_p15_fTime <= lastTimepoint);

%clear lastTimepoint lowTrimmed_p60_cTime lowTrimmed_p60_fTime lowTrimmed_p15_cTime lowTrimmed_p15_fTime;

    

% 2. calculate mean of constant environment
mean_p60_cSize = mean(p60_cSizes);
mean_p15_cSize = mean(p15_cSizes);
    
%clear p60_cSizes p15_cSizes;
    


% 3. normalize fluc data by constant mean
norm_p60_fSize = p60_fSizes./mean_p60_cSize;
norm_p15_fSize = p15_fSizes./mean_p15_cSize;

%clear p60_fDurations p15_fDurations;



% 4. plot distributions of 15 min and 60 min envir on same histogram

figure(4)
histogram(norm_p60_fSize(norm_p60_fSize<1.8),'Normalization', 'probability', 'BinWidth',0.02)
hold on
histogram(norm_p15_fSize(norm_p15_fSize<1.8),'Normalization', 'probability', 'BinWidth',0.02)
xlabel('size at birth, normalized by constant mean')
legend('60 min period','15 min period')


figure(5)
histogram(p60_fSizes(p60_fSizes<5),'Normalization', 'probability', 'BinWidth', 0.1)
hold on
histogram(p15_fSizes(p15_fSizes<5),'Normalization', 'probability', 'BinWidth', 0.1)
hold on
histogram(p60_cSizes(p60_cSizes<5),'Normalization', 'probability', 'BinWidth', 0.1)
hold on
histogram(p15_cSizes(p15_cSizes<5),'Normalization', 'probability', 'BinWidth', 0.1)
xlabel('absolute size at birth')
legend('60 min, fluc', '15 min, fluc', '60 min, const', '15 min, fluc')


