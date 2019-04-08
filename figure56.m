%% figure 56: quantification of growth rate fluctuations
%             bins growth rates into period-sized bins
%             (56_B uses hr bins across the board)

%  Goals: quantify growth rate fluctuations
%           1) time to stabilization (autocorrelation of successive periods)
%           2) mean growth rate at stabilization
%           3) mean high and low growth rates (during high and low phases)
%           4) amplitude (difference between high/low and mean)


%  Strategy:
%



%  last updated: jen, 2019 April 8

%  commit: first commit, first quantifications


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

for e = 1:length(exptArray)
    
    
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
            binsPerHour = 3600/timescale; % one bin = one period
            quartersPerHour = binsPerHour*4;
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

save('meanData.mat','meanData')


%% Part 2. plot mean mu over time

clear
clc

% 0. initialize mean data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('meanData.mat')
exptArray = [2:4,5:7,9:12,13:15];

palette = {'DodgerBlue','DeepSkyBlue','Navy','DarkCyan'};
exp_data = cell(size(meanData));

counter = 0;
for index = 9:12
    
    counter = counter + 1;
    
    fluc = meanData{index}.fluc;
    growthRates = meanData{index}.fluc_growthRates;
    signal = meanData{index}.fluc_signal;
    
    if index == 9
        fluc = meanData{index}.fluc(1:end-1);
        growthRates = meanData{index}.fluc_growthRates(1:end-1);
        signal = meanData{index}.fluc_signal(1:end-1);
    end
    
    % calculate mean high and mean low per bin as a percent difference from total mean
    pdifs = [];
    for b = 1:length(fluc)
        
        currentMean = fluc(b);
        currentMus = growthRates{b};
        currentSignals = signal{b};
        
        bins = double(currentSignals);
        bins(bins == 0) = 2;
        
        splitMus = accumarray(bins,currentMus,[],@mean);
        percentDif = (splitMus-currentMean)/currentMean;
        pdifs = [pdifs; percentDif(1), percentDif(2)];
        
    end
    exp_data{index} = pdifs;
    
    color = rgb(palette(counter));
    
    figure(1)
    plot(pdifs(:,1),'Color',color)
    hold on
    plot(pdifs(:,2),'Color',color)
    
end



