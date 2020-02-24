% experimentStats

% goal: plot a bar graph visualizing ratio of kept and discarded tracks

% strategy:
%
%       0. initialize experiment data
%       1. for each experiment...load reject data
%               2. convert cells into numbers
%               3. for each condition... 
%                       4. combine appropriate xys for full condition stats
%                       5. plot bar graph of condition breakdown
%               6. repeat for all conditions
%       7.  repeat for all experiments



% last updated: 2017 Dec 4

%% 0. initialize experiment data
clear
clc

cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')

dataIndex = find(~cellfun(@isempty,storedMetaData));
bioProdRateData = cell(size(storedMetaData));

% initialize summary vectors for calculated data
experimentCount = length(dataIndex);

%% 1. for each experiment, move to folder and load data
counter = 0;

for e = 1:experimentCount
    
    % identify experiment by date
    index = dataIndex(e);
    date = storedMetaData{index}.date;
    
    % move directory to experiment data
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    
    % load data
    timescale = storedMetaData{index}.timescale;
    if ischar(timescale) == 0
        filename = strcat('lb-fluc-',date,'-window5-width1p4-1p7-jiggle-0p5.mat');
    elseif strcmp(date,'2017-09-26') == 1
        filename = 'lb-monod-2017-09-26-window5-va-jiggle-c12-0p1-c3456-0p5-bigger1p8.mat';
    elseif strcmp(date, '2017-11-09') == 1
        filename = 'lb-control-2017-11-09-window5-width1p4-jiggle-0p5.mat';
    end
    load(filename,'D','D5','rejectD')
    
    % print current experiment
    display(strcat('Experiment (', num2str(e),') of (', num2str(length(dataIndex)),')'))
    
    % 2. convert cells into numbers
    totalTracks = cellfun(@length,D);
    finalTracks = cellfun(@length,D5);
    rejectTracks = cellfun(@length,rejectD);
    clear filename experimentFolder D D5 rejectD
   
    % 3. for each condition...
    xys = storedMetaData{index}.xys;
    xy_dimensions = size(xys);
    totalConditions = xy_dimensions(1);
    
    for c = 1:totalConditions
        
        counter = counter + 1;
        
        % 4. combine appropriate xys for full condition stats
        conditionXYs = xys(c,:);
        conditionTotal = sum(totalTracks(conditionXYs));
        conditionFinal = sum(finalTracks(conditionXYs));
        conditionRejects = sum(rejectTracks(:,conditionXYs),2);
        conditionStats = [conditionFinal, conditionRejects'];
        
        stats(counter,:) = conditionStats;
    end
    
    % 5. plot bar graph of condition breakdown
    % Display one bar for each row of the matrix.
    % The height of each bar is the sum of the elements in the row.
    figure(e)
    bar(stats,'stacked')
    hold on
    legend('kept','jump','short','noisy','small')
    
end




%%
