% figure 43

%  Goal: compute the T-statistic to determine similarity between stable conditions.
%        Does statistic support that mean of stable average is different from
%  	     mean of stable high, for instance?

%        H0: sub-populations are similar
%        H1: T-stat falls within "critical region" set by degress of freedom. 
%            if so, reject that the contrasts of the means is zero.
        

%  Conceptual strategy:
%       
%       1. Gather replicates by corresponding condition, WITHOUT AGGREGATION
%               e.g. all stable high mus belonging to each replicate are a cell.
%               i.e. there will be twelve aggregate populations per condition, 
%                    three from each of the four timescales
%
%       2. Fix "contrast coefficients", c_k , for each replicate, such that
%          their sum equals zero. Try a few different weights to make sure
%          an arbitrary one changes the results.
%               i.e. [1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1], alt each point
%                    [1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1], alt each pair
%                    [1, 1, 1, -1, -1, -1, 1, 1, 1, -1, -1, -1], alt each triplet
%                    [1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1], alt every six
%
%       3. Calculate T-statistic:
%               i.e. the sum of weighted means, over a term that accounts
%                    for the different standard deviations across replicates.
%                    see below for more details.
%
%       4. Calculate the number of degrees of freedom:
%               See below for details.
%
%
%          BIG thank yous to Natasha and Vicente for discussions and helpful to formulate this test!
%




%  Last edit: jen, 2018 Oct 11
%  commit: T-statistic to compare replicates between environments



% OK let's go!

%% 0. initialize data and growth rate parameters

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));


% 0. define growth rates of interest, see comments below for details
prompt = 'Enter specific growth rate definition as string (raw / norm / log2 / lognorm): ';
specificGrowthRate = input(prompt);
clear prompt


% 0. initialize data collection
growthRT_replicates{12,4} = [];

%% 1. aggregate data for calculations

% create array of experiments of interest
% i.e. index numbers between initial and final indeces

replicate = 0;
exptArray = [2:4,5:7,9,11,12,13:15];  % use corresponding dataIndex values

for e = 1:length(exptArray)
    
    replicate = replicate + 1;
    
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

        
        
        % 11. compile growth rates by condition into a 12 x 4 cell structure
        
        %                fluc    low     ave    high
        %       rep1     {u1}    {u1}    {u1}   {u1}
        %       rep2     {u2}    {u2}    {u2}   {u2}
        %       rep3     {u3}    {u3}    {u3}   {u3}
        %       rep4     {u4}    {u4}    {u4}   {u4}
        %       ...      ...      ...     ...    ...
        %       rep12    {u12}   {u12}   {u12}  {u12} 

        growthRT_replicates{replicate,condition} = growthRt_noNaNs; % condition = column
        
    end
    
end

clear growthRt_noNaNs growthRt growthRates_lagTrimmed growthRates_bubbleTrimmed growthRates
clear conditionData conditionData_bubbleTrimmed conditionData_fullOnly conditionData_lagTrimmed
clear D5 T timescale filename experimentFolder date e condition expType

save('growthRate_replicates.mat', 'growthRT_replicates')

%% 2. Fix "contrast coefficients", c_k , for each timescale dataset, such that their sum equals zero.

clear
clc

cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
load('growthRate_replicates.mat')


alternative_contrasts(:,1) = [1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1];
alternative_contrasts(:,2) = [1; 1; -1; -1; 1; 1; -1; -1; 1; 1; -1; -1]; 
alternative_contrasts(:,3) = [1; 1; 1; -1; -1; -1; 1; 1; 1; -1; -1; -1];
alternative_contrasts(:,4) = [1; 1; 1; 1; 1; 1; -1; -1; -1; -1; -1; -1];



%% 3. Calculate T-statistic, per pair of conditions (fluc, low, ave, high)

% T = top / bottom

% top = sum(contrast coeff * pop ave mu) ... of each aggregate timescale data set
% bottom = accounts for different standard deviations between k

% for each condition (column in growthRT_aggregates)


% (top)
% compute mean of each aggregate
replicate_means = cellfun(@mean,growthRT_replicates); % 12 x 4, replicates x conditions
replicate_counts = cellfun(@length,growthRT_replicates);

% weight aggregate means by contrast coefficient
%contrast = nan(12,4);
%for a = 1:4
%    contrast(:,a) = [alternative_contrasts(:,a); alternative_contrasts(:,a)];
%end

for conditions = 1:4
    % 12 x 4, replicates x conditions
    cu_alt1(:,conditions) = replicate_means(:,conditions).*alternative_contrasts(:,1);
    cu_alt2(:,conditions) = replicate_means(:,conditions).*alternative_contrasts(:,2);
    cu_alt3(:,conditions) = replicate_means(:,conditions).*alternative_contrasts(:,3);
    cu_alt4(:,conditions) = replicate_means(:,conditions).*alternative_contrasts(:,4);
end
clear conditions a

% sum weighted means
k_sums(1,:) = sum(cu_alt1);
k_sums(2,:) = sum(cu_alt2);
k_sums(3,:) = sum(cu_alt3);
k_sums(4,:) = sum(cu_alt4);
% k_sums is a 4x4 matrix with columns as a condition
% rows are alternative contrast coefficient schemes


% (bottom)
% compute sample variance of each aggregate
replicate_variance = nan(12,4); % 12 x 4, replicates x conditions
for index = 1:(12*4)
    differences = growthRT_replicates{index} - replicate_means(index);
    differences_sq = differences.^2;
    div_by_sampleSize = differences_sq./(replicate_counts(index)-1);
    replicate_variance(index) = sum(div_by_sampleSize);
end
% cellfun(@var,growthRT_replicates) is the same :)

% weight aggregate sample variances by contrast coefficients and sample size

weights_alt1 = [alternative_contrasts(:,1), alternative_contrasts(:,1), alternative_contrasts(:,1), alternative_contrasts(:,1)];
weights_alt2 = [alternative_contrasts(:,2), alternative_contrasts(:,2), alternative_contrasts(:,2), alternative_contrasts(:,2)];
weights_alt3 = [alternative_contrasts(:,3), alternative_contrasts(:,3), alternative_contrasts(:,3), alternative_contrasts(:,3)];
weights_alt4 = [alternative_contrasts(:,4), alternative_contrasts(:,4), alternative_contrasts(:,4), alternative_contrasts(:,4)];

weights_alt1_sq = weights_alt1.^2;
weights_alt2_sq = weights_alt2.^2;
weights_alt3_sq = weights_alt3.^2;
weights_alt4_sq = weights_alt4.^2;

div_by_n_i1 = weights_alt1_sq./12;
div_by_n_i2 = weights_alt2_sq./12;
div_by_n_i3 = weights_alt3_sq./12;
div_by_n_i4 = weights_alt4_sq./12;

weighted_variance_alt1 = div_by_n_i1.*replicate_variance;
weighted_variance_alt2 = div_by_n_i2.*replicate_variance;
weighted_variance_alt3 = div_by_n_i3.*replicate_variance;
weighted_variance_alt4 = div_by_n_i4.*replicate_variance;

% sum weighted sample variances
weighted_variance_summed(1,:) = sum(weighted_variance_alt1);
weighted_variance_summed(2,:) = sum(weighted_variance_alt2);
weighted_variance_summed(3,:) = sum(weighted_variance_alt3);
weighted_variance_summed(4,:) = sum(weighted_variance_alt4);
% weighted_variance_summed is a 4x4 matrix with columns as a condition
% rows are alternative contrast coefficient schemes

% take square root of sum
squareroot_wv_sum = sqrt(weighted_variance_summed);

% divide top by bottom
T = k_sums./squareroot_wv_sum;





% CALCULATE DEGREES OF FREEDOM
% (top)
% square sum(Ci squared * Si squared / ni)
square_wv_sum = weighted_variance_summed.^2; % 4x4 (alt scheme x condition)


% (bottom)
% weights^4 * variance^2 / counts^2 * (counts -1)
weights_quad1 = weights_alt1.^4;
weights_quad2 = weights_alt2.^4;
weights_quad3 = weights_alt3.^4;
weights_quad4 = weights_alt4.^4;

variance_sq = replicate_variance.^2;
counts_sq = 12.^2;
counts_minus1 = 12-1;

bottoms_top_alt1 = weights_quad1.*variance_sq;
bottoms_top_alt2 = weights_quad2.*variance_sq;
bottoms_top_alt3 = weights_quad3.*variance_sq;
bottoms_top_alt4 = weights_quad4.*variance_sq;

bottoms_bottom = counts_sq.*counts_minus1;

bottoms_alt1 = bottoms_top_alt1./bottoms_bottom;
bottoms_alt2 = bottoms_top_alt2./bottoms_bottom;
bottoms_alt3 = bottoms_top_alt3./bottoms_bottom;
bottoms_alt4 = bottoms_top_alt4./bottoms_bottom;

% sum bottoms
bottoms_summed(1,:) = sum(bottoms_alt1);
bottoms_summed(2,:) = sum(bottoms_alt2);
bottoms_summed(3,:) = sum(bottoms_alt3);
bottoms_summed(4,:) = sum(bottoms_alt4);
% like weighted_variance_summed & therefore square_wv_sum
% bottoms_summed is a 4x4 matrix with columns as a condition
% rows are alternative contrast coefficient schemes

% divide top and bottom
degrees_free = square_wv_sum./bottoms_summed;
%


clear bottoms_alt1 bottoms_alt2 bottoms_alt3 bottoms_alt4
clear bottoms_top_alt1 bottoms_top_alt2 bottoms_top_alt3 bottoms_top_alt4
clear cu_alt1 cu_alt2 cu_alt3 cu_alt4
clear weights_alt1 weights_alt2 weights_alt3 weights_alt4
clear weights_alt1_sq weights_alt2_sq weights_alt3_sq weights_alt4_sq
clear weights_quad1 weights_quad2 weights_quad3 weights_quad4
clear div_by_n_i1 div_by_n_i2 div_by_n_i3 div_by_n_i4
clear weighted_variance_alt1 weighted_variance_alt2 weighted_variance_alt3 weighted_variance_alt4


%%