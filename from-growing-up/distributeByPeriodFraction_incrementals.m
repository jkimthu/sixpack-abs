% timestamp-dependent distributions of dataMatrix parameters



%  Goal: plot distributions of instantaneous added mass,
%        based on period fraction

%  Last edit: Jen Nguyen, August 13th 2016



%  Sample code pulled from distributions.m and nsyncInFlux.m
%  Strategy:
%
%     0. initialize experiment and analysis parameters
%     1. isolate data of interest
%     2. accumulate data points by time bin (period fraction)
%     3. isolate periods of interest
%               a. nutrient high (fluc/const)
%               b. nutrient low   (fluc/const)
%               c. all constant
%     4. plot distributions for all isolated groups


%%
% The intended input for this script is the following data matrix,
% saved with the naming convention of:

% dmMMDD-const.mat
% dmMMDD-fluc.mat

% Initialize data.
load('dm0818-const.mat');
dC = dataMatrix;

load('dm0818-fluc.mat');
dF = dataMatrix;
clear dataMatrix;

% OK! Lez go!
%%


% 0. initiaize mass data
cumulativeMass_c = dC(:,3);
cumulativeMass_f = dF(:,3);

% 0. initiaize timestamp data
time_c = dC(:,2);
time_f = dF(:,2);

dataz{1} = cumulativeMass_c;
dataz{2} = cumulativeMass_f;
dataz{3} = time_c;
dataz{4} = time_f;

for i = 1:2 % Remove zeros
        currentVar = dataz{i};
        currentVar(currentVar <= 0) = NaN;
        nanFilter = find(~isnan(currentVar));
        currentVar = currentVar(nanFilter);
        dataz_trimmed{i} = currentVar;
        
        if i == 1
            dataz_trimmed{3} = time_c(nanFilter);
        else
            dataz_trimmed{4} = time_f(nanFilter);    
        end
end
clear duration_c duration_f addedMass_c addedMass_f timeStamp_c timeStamp_f currentVar i nanFilter;


% 0. initialize time binning parameters
periodDuration = .25;                             % duration of nutrient period in hours                 
binsPerHour = 1/periodDuration;                                % self-explanatory                       

% 0. initialize time vector for plotting
binsPerPeriod = 4;
timePerBin = 1/binsPerPeriod; 
periodTime = linspace(1, binsPerPeriod, binsPerPeriod);
periodTime = timePerBin*periodTime';                                       

% 0. initialize looping parameters for analysis
firstHour = 5;                                  % time at which to initate analysis
finalHour = 10;                                 % time at which to terminate analysis

% 1. isolate data of interest based on condition
timeStamps = dataz_trimmed{4}; %fluc
instantaneousMass_f = diff(dataz_trimmed{2});
instantaneousMass_c = diff(dataz_trimmed{1});


% trim off lower timepoints
interestingPoints = find(timeStamps >= firstHour);
interestingTimes = timeStamps(interestingPoints);
interestingMass = instantaneousMass_f(interestingPoints);

% trim off higher timepoints
interestingPoints_2 = find(interestingTimes <= finalHour);
interestingTimes_2 = interestingTimes(interestingPoints_2);
interestingMass_2 = interestingMass(interestingPoints_2);

% trim off negatives
negFilter_f = find(interestingMass_2 >= 0);
interestingTimes_3 = interestingTimes_2(negFilter_f);
interestingMass_3 = interestingMass_2(negFilter_f);

negFilter_c = find(instantaneousMass_c >= 0);
instantaneousMass_c_3 = instantaneousMass_c(negFilter_c);

% trim off crazy add-ons
posFilter_f = find(interestingMass_3 <= .2);
interestingTimes_trimmed = interestingTimes_3(posFilter_f);
interestingMass_trimmed = interestingMass_3(posFilter_f);

posFilter_c = find(instantaneousMass_c_3 <= .2);
instantaneousMass_c_trimmed = instantaneousMass_c_3(posFilter_c);


% 2. accumulate data by associated time & concentration bin
hourFractions = interestingTimes_trimmed - floor(interestingTimes_trimmed);
seqPeriods = hourFractions / timePerBin;
periodFractions = seqPeriods - floor(seqPeriods);
fractionBins = ceil(periodFractions*binsPerPeriod);


concentrationBins = zeros(length(fractionBins),1);

% for i = 1:length(fractionBins)
%     if fractionBins(i) == 1
%         concentrationBins(i) = 1;
%     elseif fractionBins(i) == 4
%         concentrationBins(i) = 1;
%     else
%         concentrationBins(i) = 2; 
%     end
% end
% clear i;

for i = 1:length(fractionBins)
    if fractionBins(i) == 1
        concentrationBins(i) = 1;
    elseif fractionBins(i) == 2
        concentrationBins(i) = 1;
    else
        concentrationBins(i) = 2; 
    end
end
clear i;

binnedByConcentration_mass = accumarray(concentrationBins,interestingMass_trimmed,[],@(x) {x});
%binnedByConcentration_mass = accumarray(fractionBins,interestingMass_trimmed,[],@(x) {x});




% 4. Plot distribution of instantaneous added sizes, low vs high
figure(2)
histogram(binnedByConcentration_mass{1},'Normalization', 'probability', 'BinWidth',0.005) %low
hold on
histogram(binnedByConcentration_mass{2},'Normalization', 'probability', 'BinWidth',0.005) %high
xlabel('instantaneous added size')
legend('low','high')

% 5. Plot distribution of instantaneous added sizes, fluc vs stable
figure(1)
histogram(interestingMass_trimmed,'Normalization', 'probability', 'BinWidth',0.005) %low
hold on
histogram(instantaneousMass_c_trimmed,'Normalization', 'probability', 'BinWidth',0.005) %high
xlabel('instantaneous added size')
legend('fluc','const')

% Yeah!