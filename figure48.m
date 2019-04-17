%% figure 48 - size as a function of the cell cycle


% Goal: Is width, length and volume changing over the course of a cell cycle?

%       This script accumulates instantaneous width by cell cycle fraction,
%       and plots a mean line with error (s.e.m.)
%
%       x axis: fraction of cell cycle
%       y axis: mean width from condition


% last edit: Jen, 2019 April 10
% commit: first commit, size over the cell cycle


% OK! Lez go!

%% Part 0. initialize data

clc
clear

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')

dataIndex = find(~cellfun(@isempty,storedMetaData));
experimentCount = length(dataIndex);
palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};


% 0. define experiments of interest
exptArray = 2:4;


% 0. define binning parameters
binsPerCC = 5;


%% Part 1. 

% 1. loop through designated experiments
for e = 1:length(exptArray)
    
    % 2. collect experiment date
    index = exptArray(e);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    expType = storedMetaData{index}.experimentType;
    disp(strcat(date, ': analyze!'))
    
    % 3. initialize experiment meta data
    xys = storedMetaData{index}.xys;
    bubbletime = storedMetaData{index}.bubbletime;
    
     % 4. load measured data
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-width1p7-jiggle-0p5.mat');
    load(filename,'D5','T');
    xy_start = min(min(xys));
    xy_end = max(max(xys));
    exptData = buildDM(D5, T, xy_start, xy_end,index,expType);

    
    % 5. for each condition,
    for condition = 1:4
        
        
        % 6. isolate condition data
        conditionData = exptData(exptData(:,21) == condition,:);
        
        
        % 7. trim data to after 3h and before bubbles
        timestamps_sec = getGrowthParameter(conditionData,'timestamp');   % ND2 file timestamp in seconds
        timestamps_hr = timestamps_sec/3600;
        
        minTime = 3;
        conditionData_post3 = conditionData(timestamps_hr >= minTime,:);
        timestamps_post3 = timestamps_hr(timestamps_hr >= minTime,:);
        
        maxTime = bubbletime(condition);
        if maxTime > 0
            conditionData_final = conditionData_post3(timestamps_post3 <= maxTime,:);
        else
            conditionData_final = conditionData_post3;
        end
        clear timestamps_hr timestamps_sec timestamps_post3
        
        
            
        % 8. remove nans from analysis (tracks without full cycles)
        ccFrac = getGrowthParameter(conditionData_final,'ccFraction');
        fulls = ~isnan(ccFrac);
        theData = conditionData_final(fulls == 1,:);
        clear ccFrac fulls 
        
        
        
        % 9. isolate width, length, volume, and cell cycle fraction data
        ww = getGrowthParameter(theData,'width');
        ll = getGrowthParameter(theData,'length');
        vv = getGrowthParameter(theData,'volume');
        ccFraction = getGrowthParameter(theData,'ccFraction');
        clear conditionData conditionData_post3 conditionData_final
        
        
        
        % 10. convert cell cycle fraction into bins
        binnedCC = ceil(ccFraction*binsPerCC);
        binnedCC(ccFraction == 0) = 1; % assign widths at birth to first bin
        clear ccFraction        
        
    
        
        % 11. bin size parameters by cell cycle fraction
        binnedW = accumarray(binnedCC,ww,[],@mean);
        binnedL = accumarray(binnedCC,ll,[],@mean);
        binnedV = accumarray(binnedCC,vv,[],@mean);
        clear ll ww vv 

        
        
        % 12. plot width, length and volume vs cell cycle fraction
        color = rgb(palette(condition));
        
        figure(1) % width
        plot(binnedW,'Color',color,'LineWidth',2)
        hold on
        
        figure(2) % length
        plot(binnedL,'Color',color,'LineWidth',2)
        hold on
        
        figure(3) % volume
        plot(binnedV,'Color',color,'LineWidth',2)
        hold on
        
        
    end
    
    
end

%% Part 2. save plots

figure(1)
xlabel('fraction of cell cycle')
ylabel('mean width')
title('width as a function of the cell cycle')

figure(2)
xlabel('fraction of cell cycle')
ylabel('mean length')
title('length as a function of the cell cycle')

figure(3)
xlabel('fraction of cell cycle')
ylabel('mean volume')
title('volume as a function of the cell cycle')

