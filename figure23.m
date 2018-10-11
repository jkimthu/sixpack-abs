% figure 23

%  Goal: compute the T-statistic to determine similarity of replicates.
%        Does statistic support our belief that fluctuating populations are 
%  	     the same between stable replicates, but different across timescales?

%        H0: sub-populations are similar
%        H1: T-stat falls within "critical region" set by degress of freedom. 
%            if so, reject that the contrasts of the means is zero.
        

%  Conceptual strategy:
%       
%       1. Combine replicates by corresponding timescale dataset
%               e.g. all stable high mus belonging to 30 sec replicates will belong
%                    to the same aggregate population.
%               i.e. there will be four aggregate populations per condition, 
%                    one for each of the four timescales (k_timescale)
%
%       2. Fix "contrast coefficients", c_k , for each timescale dataset, such that
%          their sum equals zero.
%               i.e. These are the weights for each dataset, based on
%                    cleanliness / reliability.
%
%       3. Calculate T-statistic:
%               i.e. the sum of weighted aggregate population mean, over
%                    a term that accounts for the different standard deviations
%                    across replicates. see below for more details.
%
%       4. Calculate the number of degrees of freedom:
%               See below for details.
%
%
%          BIG thank yous to Natasha for helping formulate this test!
%




%  Last edit: jen, 2018 Oct 11
%  commit: T-statistic to compare replicates within environments



% OK let's go!

%% 0. initialize data and growth rate parameters

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));


% 0. define indeces of experiment to loop through?
index_i = 2; % do not count 2017-10-31, an outlier
prompt = 'Enter final experiment index as double (meta data index of final 60 min rep): ';
index_f = input(prompt);


% 0. define growth rates of interest, see comments below for details
prompt = 'Enter specific growth rate definition as string (raw / norm / log2 / lognorm): ';
specificGrowthRate = input(prompt);
clear prompt

% 0. initialize data collection
growthRT_aggregates{4,4} = [];

%% 1. aggregate data for calculations

% create array of experiments of interest
% i.e. index numbers between initial and final indeces

intermed = dataIndex(dataIndex >= index_i); % use corresponding dataIndex values
exptArray = intermed(intermed<= index_f); 
clear intermed


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
    filename = strcat('lb-fluc-',date,'-width1p7-jiggle-0p5.mat');
    load(filename,'D5','T');
    
    
    
    % 4. build data matrix from specified condition
    for condition = 1:length(bubbletime)
        
        xy_start = storedMetaData{index}.xys(condition,1);
        xy_end = storedMetaData{index}.xys(condition,end);
        conditionData = buildDM(D5, T, xy_start, xy_end,index,expType);
        clear xy_start xy_end
        
        
        % 5. isolate condition data to those with full cell cycles
        curveIDs = conditionData(:,5);           % col 5 = curve ID
        conditionData_fullOnly = conditionData(curveIDs > 0,:);
        clear curveFinder curveIDs
        
        
        
        % 6. isolate volume (Va), timestamp, mu, drop and curveID data
        volumes = conditionData_fullOnly(:,11);        % col 11 = calculated va_vals (cubic um)
        timestamps_sec = conditionData_fullOnly(:,2);  % col 2  = timestamp in seconds
        isDrop = conditionData_fullOnly(:,4);          % col 4  = isDrop, 1 marks a birth event
        curveFinder = conditionData_fullOnly(:,5);     % col 5  = curve finder (ID of curve in condition)
        trackNum = conditionData_fullOnly(:,20);       % col 20 = track number (not ID from particle tracking)
        
        
        
        % 7. calculate growth rate
        growthRates = calculateGrowthRate(volumes,timestamps_sec,isDrop,curveFinder,trackNum);
        clear volumes isDrop curveFinder trackNum timestamps_sec
        
        
        % 8. truncate data to non-erroneous (e.g. bubbles) timestamps
        maxTime = bubbletime(condition);
        timestamps_hr = conditionData_fullOnly(:,2)/3600; % time in seconds converted to hours
        
        if maxTime > 0
            conditionData_bubbleTrimmed = conditionData_fullOnly(timestamps_hr <= maxTime,:);
            growthRates_bubbleTrimmed = growthRates(timestamps_hr <= maxTime,:);
        else
            conditionData_bubbleTrimmed = conditionData_fullOnly;
            growthRates_bubbleTrimmed = growthRates;
        end
        clear timestamps_hr maxTime
        
        
        
        % 9. truncate data to stabilized regions
        minTime = 3;
        timestamps_hr = conditionData_bubbleTrimmed(:,2)/3600; % time in seconds converted to hours
        
        conditionData_lagTrimmed = conditionData_bubbleTrimmed(timestamps_hr >= minTime,:);
        growthRates_lagTrimmed = growthRates_bubbleTrimmed(timestamps_hr >= minTime,:);
        clear timestamps_hr
        
        
        
        % 10. isolate selected specific growth rate and remove nans from data analysis
        if strcmp(specificGrowthRate,'raw') == 1
            specificColumn = 1;         % for selecting appropriate column in growthRates
        elseif strcmp(specificGrowthRate,'norm') == 1
            specificColumn = 2;
        elseif strcmp(specificGrowthRate,'log2') == 1
            specificColumn = 3;
        elseif strcmp(specificGrowthRate,'lognorm') == 1
            specificColumn = 4;
        end
        
        growthRt = growthRates_lagTrimmed(:,specificColumn);
        growthRt_noNaNs = growthRt(~isnan(growthRt),:); % nans occur at drops and track-to-track transitions
        %length(growthRt_noNaNs)
        clear gr_median gr_std gr_temp

        
        
        % 11. compile growth rates by timescale and condition into a 4x4 cell structure
        
        %                fluc    low     ave    high
        %       k_30s    {u1}    {u1}    {u1}   {u1}
        %       k_5m     {u2}    {u2}    {u2}   {u2}
        %       k_15m    {u3}    {u3}    {u3}   {u3}
        %       k_60m    {u4}    {u4}    {u4}   {u4}
        
        if timescale == 30
            row = 1;
        elseif timescale == 300
            row = 2;
        elseif timescale == 900
            row = 3;
        elseif timescale == 3600
            row = 4;
        end
        
        % concatenate current condition growth rates to existing array of rates (if any)
        interm = growthRT_aggregates{row,condition}; % condition = column
        growthRT_aggregates{row,condition} = [interm; growthRt_noNaNs]; 
        clear interm row
        
    end
    
end

clear growthRt_noNaNs growthRt growthRates_lagTrimmed growthRates_bubbleTrimmed growthRates
clear conditionData conditionData_bubbleTrimmed conditionData_fullOnly conditionData_lagTrimmed
clear D5 T timescale filename experimentFolder date e condition expType

save('growthRate_aggregates.mat', 'growthRT_aggregates')

%% 2. Fix "contrast coefficients", c_k , for each timescale dataset, such that their sum equals zero.

clear
clc

cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
load('growthRate_aggregates.mat')

ratio = 2;
signs = [-1;1;-1;1];
contrasts = [ratio; 1; 1; ratio].*signs;



%% 3. Calculate T-statistic, per condition (fluc, low, ave, high)

% T = top / bottom

% top = sum(contrast coeff * pop ave mu) ... of each aggregate timescale data set
% bottom = accounts for different standard deviations between k

% for each condition (column in growthRT_aggregates)
%for aggregate = 1:4
    
    % (top)
    % compute mean of each aggregate
    aggregate_means = cellfun(@mean,growthRT_aggregates);
    aggregate_counts = cellfun(@length,growthRT_aggregates);
    
    % weight aggregate means by contrast coefficient
    for k = 1:4
        cu(:,k) = aggregate_means(:,k).*contrasts;
    end
    clear k
    
    % sum weighted means
    k_sums = sum(cu);
    
    
    % (bottom)
    % compute sample variance of each aggregate
    aggregate_variance = nan(4,4);
    for index = 1:16
        differences = growthRT_aggregates{index} - aggregate_means(index);
        differences_sq = differences.^2;
        div_by_sampleSize = differences_sq./(aggregate_counts(index)-1);
        aggregate_variance(index) = sum(div_by_sampleSize);
    end
    
    
    % weight aggregate sample variances by contrast coefficients and sample size
    weights = [contrasts, contrasts, contrasts, contrasts];
    weights_sq = weights.^2;
    div_by_n_i = weights_sq./aggregate_counts;
    weighted_variance = div_by_n_i.*aggregate_variance;
    
    % sum weighted sample variances
    weighted_variance_summed = sum(weighted_variance);
    
    % take square root of sum
    squareroot_wv_sum = sqrt(weighted_variance_summed);
    
    % divide top by bottom
    T = k_sums./squareroot_wv_sum;

    
    
   
    
    % CALCULATE DEGREES OF FREEDOM
    % (top)
    % square sum(Ci squared * Si squared / ni)
    square_wv_sum = weighted_variance_summed.^2;

    
    % (bottom)
    % weights^4 * variance^2 / counts^2 * (counts -1)
    weights_quad = weights.^4;
    variance_sq = aggregate_variance.^2;
    counts_sq = aggregate_counts.^2;
    counts_minus1 = aggregate_counts-1;
    bottoms_top = weights_quad.*variance_sq;
    bottoms_bottom = counts_sq.*counts_minus1;
    bottoms = bottoms_top./bottoms_bottom;
    
    % sum bottoms
    bottoms_summed = sum(bottoms);
   
    % divide top and bottom
    degrees_free = square_wv_sum./bottoms_summed;
%

%% 4. Calculate T and degrees of freedom







