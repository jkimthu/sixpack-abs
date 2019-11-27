% timestamp-dependent distributions



%  Goal: plot distributions of cell cycle duration and added mass,
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

% dFMMDD.mat
% dCMMDD.mat

% Initialize data.
load('dC0818_mit.mat');
load('dF0818_mit.mat');
dC = dC0818_mit;
dF = dF0818_mit;
clear dC0818_mit dF0818_mit;

% OK! Lez go!

%%



% 0. initiaize duration data
duration_c = dC(:,1);
duration_f = dF(:,1);

% 0. initiaize mass data
addedMass_c = dC(:,2);
addedMass_f = dF(:,2);

% 0. initiaize timestamp data
timeStamp_c = dC(:,3);
timeStamp_f = dF(:,3);

dataz{1} = duration_c;
dataz{2} = duration_f;
dataz{3} = addedMass_c;
dataz{4} = addedMass_f;
dataz{5} = timeStamp_c;
dataz{6} = timeStamp_f;

for i = 1:2 % Remove zeros
        currentVar = dataz{i};
        currentVar(currentVar <= 0) = NaN;
        nanFilter = find(~isnan(currentVar));
        currentVar = currentVar(nanFilter);
        dataz_trimmed{i} = currentVar;
        
        if i == 1
            dataz_trimmed{3} = addedMass_c(nanFilter);
            dataz_trimmed{5} = timeStamp_c(nanFilter);
        else
            dataz_trimmed{4} = addedMass_f(nanFilter);
            dataz_trimmed{6} = timeStamp_f(nanFilter);    
        end
end
clear duration_c duration_f addedMass_c addedMass_f timeStamp_c timeStamp_f currentVar i nanFilter;


% 0. initialize time binning parameters
periodDuration = 1;                             % duration of nutrient period in hours                 
binsPerHour = 4;                                % self-explanatory 
hrPerBin = 1/binsPerHour;                       % time bins of 0.005 hr

% 0. initialize time vector for plotting
binsPerPeriod = periodDuration/hrPerBin;
periodTime = linspace(1, binsPerPeriod, binsPerPeriod);
periodTime = hrPerBin*periodTime';                                       

% 0. initialize looping parameters for analysis
firstHour = 5;                                  % time at which to initate analysis
finalHour = 10;                                 % time at which to terminate analysis

% 1. isolate data of interest based on condition
timeStamps = dataz_trimmed{6}; %fluc
addedMass = dataz_trimmed{4};
cycleDuration = dataz_trimmed{2};

% trim off lower timepoints
interestingPoints = find(timeStamps >= firstHour);
interestingTimes = timeStamps(interestingPoints);
interestingMass = addedMass(interestingPoints);
interestingDurations = cycleDuration(interestingPoints);

% trim off higher timepoints
interestingPoints_2 = find(interestingTimes <= finalHour);
interestingTimes_trimmed = interestingTimes(interestingPoints_2);
interestingMass_trimmed = interestingMass(interestingPoints_2);
interestingDurations_trimmed = interestingDurations(interestingPoints_2);


% 2. accumulate data by associated time & concentration bin
periodFractions = interestingTimes_trimmed - floor(interestingTimes_trimmed);
fractionBins = ceil(periodFractions*binsPerHour);


concentrationBins = zeros(length(fractionBins),1);

for i = 1:length(fractionBins)
    if fractionBins(i) == 1
        concentrationBins(i) = 1;
    elseif fractionBins(i) == 4
        concentrationBins(i) = 1;
    else
        concentrationBins(i) = 2; 
    end
end
clear i;

binnedByConcentration_mass = accumarray(concentrationBins,interestingMass_trimmed,[],@(x) {x});
binnedByConcentration_duration = accumarray(concentrationBins,interestingDurations_trimmed,[],@(x) {x});

% Yeah!

% 4. Plot distribution of added sizes per cell cycle
figure(1)
histogram(binnedByConcentration_mass{1},'Normalization', 'probability', 'BinWidth',0.2) %low
hold on
histogram(binnedByConcentration_mass{2},'Normalization', 'probability', 'BinWidth',0.2) %high
%hold on
%histogram(dataz_trimmed{3},'Normalization', 'probability', 'BinWidth',0.1) %high
xlabel('added size per cell cycle')
legend('low','high','stable ave')



% 6. Plot distribution of cell cycle durations
figure(2)
histogram(binnedByConcentration_duration{1},'Normalization', 'probability','BinWidth',0.2) %low
hold on
histogram(binnedByConcentration_duration{2},'Normalization', 'probability', 'BinWidth',0.2) %high
%hold on
%histogram(dataz_trimmed{1},'Normalization', 'probability', 'BinWidth',0.1) %high
xlabel('cell cycle durations')
legend('low','high','stable ave')

