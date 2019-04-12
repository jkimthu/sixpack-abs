%% figure 56: quantification of growth rate fluctuations
%             bins growth rates into 1hr bins, regardless of timescale

%  Goals: quantify growth rate fluctuations
%           1) time to stabilization (autocorrelation of successive periods)
%           2) mean growth rate at stabilization
%           3) mean high and low growth rates (during high and low phases)
%           4) amplitude (difference between high/low and mean)


%  Strategy:
%



%  last updated: jen, 2019 April 12

%  commit: plot all conditions as mean of replicates with shaded error (standard dev)


% OK let's go!

%% initialize

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));


% 0. define growth rate of interest
specificGrowthRate = 'log2';
specificColumn = 3;


% 0. initialize new storate structure for quantified data 
meanData = cell(size(storedMetaData));


%% Part 1. loop through experiments and record data vectors of binned mu and nutrient signal

% 1. create array of experiments of interest, then loop through each:
exptArray = [2:4,5:7,9:12,13:15]; % use corresponding dataIndex values

for e = 4:6
    
    
    % 2. initialize experiment meta data
    index = exptArray(e);
    date = storedMetaData{index}.date;
    expType = storedMetaData{index}.experimentType;
    bubbletime = storedMetaData{index}.bubbletime;
    timescale = storedMetaData{index}.timescale;
    disp(strcat(date, ': analyze!'))
    
    
    
    % 3. load measured experiment data    
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    if strcmp(date,'2017-11-09') == 1
        filename = strcat('lb-control-',date,'-width1p4-jiggle-0p5.mat');
    elseif strcmp(date,'2017-09-26') == 1
        filename = 'lb-monod-2017-09-26-jiggle-c12-0p1-c3456-0p5-bigger1p8.mat';
    elseif strcmp(date,'2018-12-04') == 1
        filename = 'lb-monod-2018-12-04-c12-width1p7-c34-width1p4-jiggle-0p5.mat';
    elseif strcmp(expType,'origFluc') == 1
        filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');

        % note: this code is not written for single shift or fluc-to-stable data
        
    end
    load(filename,'D5','T');
    
    

    % 4. build data matrix from specified condition
    for condition = 1:length(bubbletime)
        
        
        xy_start = storedMetaData{index}.xys(condition,1);
        xy_end = storedMetaData{index}.xys(condition,end);
        conditionData = buildDM(D5, T, xy_start, xy_end,index,expType);
        
        
        
        % 5. isolate volume (Va), timestamp, drop, curve, and trackNum data
        volumes = getGrowthParameter(conditionData,'volume');             % volume = calculated va_vals (cubic um)
        timestamps_sec = getGrowthParameter(conditionData,'timestamp');   % ND2 file timestamp in seconds
        isDrop = getGrowthParameter(conditionData,'isDrop');              % isDrop == 1 marks a birth event
        curveFinder = getGrowthParameter(conditionData,'curveFinder');    % col 5  = curve finder (ID of curve in condition)
        trackNum = getGrowthParameter(conditionData,'trackNum');          % track number, not ID from particle tracking
        
        
        
        % 6. calculate growth rate
        growthRates = calculateGrowthRate(volumes,timestamps_sec,isDrop,curveFinder,trackNum);
        clear volumes isDrop curveFinder trackNum
        
        

        % 7. truncate data to non-erroneous (e.g. bubbles) timestamps
        maxTime = bubbletime(condition);
        timestamps_hr = conditionData(:,2)/3600; % time in seconds converted to hours
        
        if maxTime > 0
            conditionData_bubbleTrimmed = conditionData(timestamps_hr <= maxTime,:);
            growthRates_bubbleTrimmed = growthRates(timestamps_hr <= maxTime,:);
        else
            conditionData_bubbleTrimmed = conditionData;
            growthRates_bubbleTrimmed = growthRates;
        end
        clear timestamps_hr timestamps_sec
        
        
        
        
        % 8. isolate selected specific growth rate and timestamp (corrected if condition 1)
        growthRt = growthRates_bubbleTrimmed(:,specificColumn);
        if strcmp(date,'2017-10-10')
            time_sec = getGrowthParameter(conditionData_bubbleTrimmed,'timestamp');
        elseif condition == 1
            time_sec = getGrowthParameter(conditionData_bubbleTrimmed,'correctedTime'); % ND2 file timestamp in seconds
        else
            time_sec = getGrowthParameter(conditionData_bubbleTrimmed,'timestamp');
        end
        timeInHours = time_sec./3600;
        clear conditionData time_sec
        
        
        

        % 9. remove NaNs from growth rate data
        growthRt_noNaNs = growthRt(~isnan(growthRt),:);
        timeInHours_noNans = timeInHours(~isnan(growthRt),:);
        clear growthRt timeInHours growthRates
        
        
        
        % 10. calculate binary nutrient signals
        if condition == 1
            binsPerHour = 1; % one bin per hour
            periodsPerHour = 3600/timescale;
            quartersPerHour = periodsPerHour*4;
            timeFractions = timeInHours_noNans - floor(timeInHours_noNans);
            timeQuarters = ceil(timeFractions*quartersPerHour);
            
            bns = nan(length(timeQuarters),1);
            if timescale == 300
                
                highs = [1,4,5,8,9,12,13,16,17,20,21,24,25,28,29,32,33,36,37,40,41,44,45,48];
                bns = ismember(timeQuarters,highs);
                
            elseif timescale == 900
                
                highs = [1,4,5,8,9,12,13,16];
                bns = ismember(timeQuarters,highs);
                
            elseif timescale == 3600
                
                bns(timeQuarters == 1) = 1;
                bns(timeQuarters == 2) = 0;
                bns(timeQuarters == 3) = 0;
                bns(timeQuarters == 4) = 1;
            end
            
        end
        
        
        
        % 11. bin growth rate and binary nutrient signal into period-sized bins
        bins = ceil(timeInHours_noNans*binsPerHour);
        binned_growthRt = accumarray(bins,growthRt_noNaNs,[],@(x) {x});
        if condition == 1
            binned_signal = accumarray(bins,bns,[],@(x) {x});
        end
        
        
        
        % 12. calculate mean across bins
        bin_means = cellfun(@mean,binned_growthRt);
        if condition == 1
            meanData{exptArray(e)}.fluc = bin_means;
            meanData{exptArray(e)}.fluc_growthRates = binned_growthRt;
            meanData{exptArray(e)}.fluc_signal = binned_signal;
        elseif condition == 2
            meanData{exptArray(e)}.low = bin_means;
        elseif condition == 3
            meanData{exptArray(e)}.ave = bin_means;
        else
            meanData{exptArray(e)}.high = bin_means;
        end
        

    end
end

save('meanData_5MIN.mat','meanData')


%% Part 2. plot mean mu over time of individual replicates

clear
clc

% 0. initialize mean data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('meanData_B.mat')
exptArray = [5:7,9:12,13:15];
colorArray = [1,1,1,2,2,2,2,3,3,3];

palette = {'GoldenRod','ForestGreen','MediumPurple'};
exp_data = cell(size(meanData));


for e = 1:length(exptArray)
    
    index = exptArray(e);
    
    fluc = meanData{index}.fluc;
    growthRates = meanData{index}.fluc_growthRates;
    signal = meanData{index}.fluc_signal;
    
    if index == 9
        fluc = meanData{index}.fluc(1:end-1);
        growthRates = meanData{index}.fluc_growthRates(1:end-1);
        signal = meanData{index}.fluc_signal(1:end-1);
    end
    
    % calculate mean high and mean low per bin as a percent difference from total mean
    pDif = [];
    for b = 1:length(fluc)
        
        currentMean = fluc(b);
        currentMus = growthRates{b};
        currentSignals = signal{b};
        
        bins = double(currentSignals);
        bins(bins == 0) = 2;
        
        splitMus = accumarray(bins,currentMus,[],@mean);
        percentDif = (splitMus-currentMean)/currentMean;
        pDif = [pDif; percentDif(1), percentDif(2)];
        
    end
    exp_data{index} = pDif;
    
    
    color = rgb(palette(colorArray(e)));
    
    figure(1)
    subplot(1,2,1)
    plot(pDif(:,1),'Color',color)
    hold on
    subplot(1,2,2)
    plot(pDif(:,2),'Color',color)
    hold on
    
    figure(2)
    plot(fluc,'Color',color)
    hold on
    
end
figure(1)
subplot(1,2,1)
ylabel('% difference from mean growth rate')
title('high nutrient phase')
xlabel('hours from start')
subplot(1,2,2)
title('low nutrient phase')

figure(2)
ylabel('mean growth rate')
xlabel('time (h)')
title('mean growth rate binned every hour')
axis(1,8,0.5,1.8)


%% Part 3. plot mean and error of replicates, including steady low, ave and high

clear
clc

% 0. initialize mean data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('meanData_B.mat')
exptArray = [2:4,5:7,9:12,13:15];
colorArray = [1,1,1,2,2,2,3,3,3,3,4,4,4,4];
numReps = [3,3,4,3];


% 0. initialize color designations
palette_fluc = {'LightSkyBlue','SteelBlue','DeepSkyBlue','Navy'};
palette_stable = {'Indigo','DarkGoldenRod','DarkRed'};


counter = 0;
for e = 1:length(exptArray)
    
    
    counter = counter + 1;
    timescale = colorArray(counter);
    index = exptArray(e);
    
    if timescale == 3
        maxTime = 4;
    else
        maxTime = 5;
    end
    
    fluc = meanData{index}.fluc;
    low = meanData{index}.low;
    ave = meanData{index}.ave;
    high = meanData{index}.high;
    
    if length(fluc) > maxTime
        fluc = meanData{index}.fluc(1:maxTime);
    elseif length(fluc) < maxTime
        f2 = nan(maxTime,1);
        f2(1:length(fluc),1) = fluc;
        fluc = f2;
    end
    
    if length(low) > 5
        low = meanData{index}.low(1:5);
    elseif length(low) < 5
        l2 = nan(5,1);
        l2(1:length(low),1) = low;
        low = l2;
    end
    
    if length(ave) > 5
        ave = meanData{index}.ave(1:5);
    elseif length(ave) < 5
        a2 = nan(5,1);
        a2(1:length(ave),1) = ave;
        ave = a2;
    end
    
    if length(high) > 5
        high = meanData{index}.high(1:5);
    elseif length(high) < 5
        h2 = nan(5,1);
        h2(1:length(high),1) = high;
        high = h2;
    end
    clear l2 a2 h2
    
    % compile 2d matrix of mean replicate data, with rows as replicate and columns being time (period bin)
    if timescale == 1
        rep = e;
        replicate_30_fluc(rep,:) = fluc;
%         replicate_30_low(rep,:) = low;
%         replicate_30_ave(rep,:) = ave;
%         replicate_30_high(rep,:) = high;
%         
    elseif timescale == 2
        rep = e-3;
        rep_5_fluc(rep,:) = fluc;
%         rep_5_low(rep,:) = low;
%         rep_5_ave(rep,:) = ave;
%         rep_5_high(rep,:) = high;
%         
    elseif timescale == 3
        maxTime = 4; % maxTime for 15 min condition is 4 hrs due to data quality
        rep = e-6;
        rep_15_fluc(rep,:) = fluc;
%         rep_15_low(rep,:) = low;
%         rep_15_ave(rep,:) = ave;
%         rep_15_high(rep,:) = high;
        
    elseif timescale == 4
        rep = e-10;
        rep_60_fluc(rep,:) = fluc;
%         rep_60_low(rep,:) = low;
%         rep_60_ave(rep,:) = ave;
%         rep_60_high(rep,:) = high;
        
    end
    
    rep_low(e,:) = low;
    rep_ave(e,:) = ave;
    rep_high(e,:) = high;
    
    clear rep fluc high ave low
    
end


% calculate mean and stdev across replicates (rows)
r30_means = nanmean(replicate_30_fluc);
r30_std = nanstd(replicate_30_fluc);

r5_means = nanmean(rep_5_fluc);
r5_std = nanstd(rep_5_fluc);

r15_means = nanmean(rep_15_fluc);
r15_std = nanstd(rep_15_fluc);

r60_means = nanmean(rep_60_fluc);
r60_std = nanstd(rep_60_fluc);


% plot

figure(1) % fluc
ss = shadedErrorBar([1, 2, 3, 4, 5],replicate_30_fluc,{@nanmean,@nanstd},'lineprops',{'Color',rgb(palette_fluc(1))},'patchSaturation',0.3);
ss.mainLine.LineWidth = 3;
hold on
ss = shadedErrorBar([1, 2, 3, 4, 5],rep_5_fluc,{@nanmean,@nanstd},'lineprops',{'Color',rgb(palette_fluc(2))},'patchSaturation',0.3);
ss.mainLine.LineWidth = 3;
hold on
ss = shadedErrorBar([1, 2, 3, 4],rep_15_fluc,{@nanmean,@nanstd},'lineprops',{'Color',rgb(palette_fluc(3))},'patchSaturation',0.3);
ss.mainLine.LineWidth = 3;
hold on
ss = shadedErrorBar([1, 2, 3, 4, 5],rep_60_fluc,{@nanmean,@nanstd},'lineprops',{'Color',rgb(palette_fluc(4))},'patchSaturation',0.3);
ss.mainLine.LineWidth = 3;

title('fluctuating: time-averaged growth rate binned every hour')
ylabel('Growth rate (1/h)')
xlabel('Time (h)')



figure(2) % stable
ss = shadedErrorBar([1, 2, 3, 4, 5],rep_low,{@nanmean,@nanstd},'lineprops',{'Color',rgb(palette_stable(1))},'patchSaturation',0.3);
ss.mainLine.LineWidth = 3;
hold on
ss = shadedErrorBar([1, 2, 3, 4, 5],rep_ave,{@nanmean,@nanstd},'lineprops',{'Color',rgb(palette_stable(2))},'patchSaturation',0.3);
ss.mainLine.LineWidth = 3;
hold on
ss = shadedErrorBar([1, 2, 3, 4, 5],rep_high,{@nanmean,@nanstd},'lineprops',{'Color',rgb(palette_stable(3))},'patchSaturation',0.3);
ss.mainLine.LineWidth = 3;
title('steady: time-averaged growth rate binned every hour')
ylabel('Growth rate (1/h)')
xlabel('Time (h)')



figure(3) % all
ss = shadedErrorBar([1, 2, 3, 4, 5],replicate_30_fluc,{@nanmean,@nanstd},'lineprops',{'Color',rgb(palette_fluc(1))},'patchSaturation',0.3);
ss.mainLine.LineWidth = 3;
hold on
ss = shadedErrorBar([1, 2, 3, 4, 5],rep_5_fluc,{@nanmean,@nanstd},'lineprops',{'Color',rgb(palette_fluc(2))},'patchSaturation',0.3);
ss.mainLine.LineWidth = 3;
hold on
ss = shadedErrorBar([1, 2, 3, 4],rep_15_fluc,{@nanmean,@nanstd},'lineprops',{'Color',rgb(palette_fluc(3))},'patchSaturation',0.3);
ss.mainLine.LineWidth = 3;
hold on
ss = shadedErrorBar([1, 2, 3, 4, 5],rep_60_fluc,{@nanmean,@nanstd},'lineprops',{'Color',rgb(palette_fluc(4))},'patchSaturation',0.3);
ss.mainLine.LineWidth = 3;
hold on
ss = shadedErrorBar([1, 2, 3, 4, 5],rep_low,{@nanmean,@nanstd},'lineprops',{'Color',rgb(palette_stable(1))},'patchSaturation',0.3);
ss.mainLine.LineWidth = 3;
hold on
ss = shadedErrorBar([1, 2, 3, 4, 5],rep_ave,{@nanmean,@nanstd},'lineprops',{'Color',rgb(palette_stable(2))},'patchSaturation',0.3);
ss.mainLine.LineWidth = 3;
hold on
ss = shadedErrorBar([1, 2, 3, 4, 5],rep_high,{@nanmean,@nanstd},'lineprops',{'Color',rgb(palette_stable(3))},'patchSaturation',0.3);
ss.mainLine.LineWidth = 3;

title('time-averaged growth rate binned every hour')
ylabel('Growth rate (1/h)')
xlabel('Time (h)')
    
axis([1,5,0.5,3])
%%
  








% save only single shift data
figure(3)
plotName = strcat('figure51-upshift-',specificGrowthRate,'-mean&std-singleShiftOnly');
%saveas(gcf,plotName,'epsc')


% save mean signal (across replicates)
save('response_singleUpshift.mat','upshift_means','upshift_times')
