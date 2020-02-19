% defineStable.m

% goal: function that defines the start and end of a stabilized growth
%       region, for applications that require selection of steady-state
%       growth physiology, etc.
%
%       definition of stable is when the derivative of mean population
%       curve is within the largest s.e.m. value of that curve


% strategy:

%   0. three input vectors: mean of population metric, s.e.m., and time
%   1. determine median s.e.m. value
%   2. find the derivative of mean curve
%   3. threshold derivative, such that values under threshold == 0
%   4. find longest stretch of consecutive zeros
%   5. identify bounding timepoints of this stretch
%   6. return timepoints as variables, startStable and endStable
%   7. visualize bounded stable region, to check


% last update: jen, 2017 oct 23

%%

%  0. three input vectors: mean of population metric, s.e.m., and time
%function [startStable, endStable] = defineStable(meanVector,semVector,timeVector)

meanVector = meanVa;
semVector = semVa;

% 0. smooth meanVector, with the mean of sliding window of 5 points
smoothedMean = slidefun(@mean,3,meanVector);

figure()
plot(timeVector,meanVector)
hold on
plot(timeVector,smoothedMean)
hold on
plot(timeVector,smoothedMean)

% 1. determine largest s.e.m. value
range = median(semVector);

% 2. find the derivative of mean curve
derivative = abs(diff(smoothedMean));

% 3. threshold derivative, such that values under threshold == 0
threshold = derivative <= range;

% 4. find longest stretch of consecutive ones
%       i. stretch identity (0 or 1)
ii = abs(diff(threshold));
id = threshold(ii == 1);

%      ii. determine stretch length
switches = find(diff(threshold));
ic = [0; switches];
stretchLengths = diff(ic);

%     iii. identify longest stretch of ones
theOne = max(stretchLengths .* id);
location = find(stretchLengths == theOne);

% 5. identify bounding timepoints of this stretch
endFrame = switches(location);
startFrame = switches(location-1) + 1;

% 6. return timepoints as variables, startStable and endStable
startStable = timeVector(startFrame);
endStable = timeVector(endFrame);

% 7. plot visual check
%   i. trim timepoints earlier than start
mean_test = smoothedMean(timeVector >=startStable);
time_test = timeVector(timeVector>=startStable);

%  ii. trim timepoints later than end
mean_test2 = mean_test(time_test<=endStable);
time_test2 = time_test(time_test<=endStable);

plot(time_test2,mean_test2)

%end
