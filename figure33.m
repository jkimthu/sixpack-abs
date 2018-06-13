%% figures 1 and 2


%  Goal: determine whether growth in low, ave, high and fluc environments
%        are characterizable as adder, sizer, or timer

%        plot (1) population-averaged division width vs birth width
%        plot (2) population-averaged division length vs birth length



%  Strategy: 
%
%       a) initialize experimental data
%       b) identify complete cell cycles > 10 min long within each condition 
%       c) compile birth lengths & widths, division lengths & widths
%       d) calculate mean and plot condition data
%



%  Last edit: jen, 2018 June 13

%  Commit: plot population-averaged division W or L vs birth W or L


%  OK let's go!

%% initialize

clc
clear

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')

dataIndex = find(~cellfun(@isempty,storedMetaData));
experimentCount = length(dataIndex);


% 1. for all experiments in dataset
exptCounter = 0;
for e = 1:experimentCount
    
    
    % 1. collect experiment meta data
    index = dataIndex(e);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    bubbletime = storedMetaData{index}.bubbletime;
    
    % exclude outliers from analysis (2017-10-31 and monod experiments)
    if strcmp(date, '2017-10-31') == 1 || strcmp (timescale, 'monod') == 1
        disp(strcat(date,': excluded from analysis'))
        continue
    end
    disp(strcat(date, ': analyze!'))
    exptCounter = exptCounter + 1;
    
    
    % 2. load measured data
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-window5-width1p4-1p7-jiggle-0p5.mat');
    load(filename,'D5','M','M_va','T');
    
    
    % 3. compile experiment data matrix
    xy_start = min(min(storedMetaData{index}.xys));
    xy_end = max(max(storedMetaData{index}.xys));
    exptData = buildDM(D5, M, M_va, T, xy_start, xy_end,e);
    clear D5 M M_va T xy_start xy_end e
    
    
    % 4. initialize colors for plotting
    palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
    shapes = {'o','x','square','*'};
    
    for condition = 1:length(bubbletime)
        
        
        % 5. isolate condition specific data
        conditionData = exptData(exptData(:,23) == condition,:);  % col 23 = cond vals
        
        
        % 6. trim data to full cell cycles ONLY
        curveFinder = conditionData(:,6);        % col 6 = curveFinder, ID of full cell cycles
        conditionData_fullOnly = conditionData(curveFinder > 0,:);
        clear curveFinder
        
        
        % 7. isolate interdivision time
        interdivisionTime = conditionData_fullOnly(:,8)/60;     % col 8  = curve duration (sec converted to min)
        
        
        % 8. trim data to include only full cell cycles longer than 10 min
        conditionData_trim = conditionData_fullOnly(interdivisionTime > 10,:);
        clear interdivisionTime
        
        
        % 9. isolate volume and interdiv time
        lengthVals = conditionData_trim(:,3);     % col 3  = measured length (um)
        widths = conditionData_trim(:,11);        % col 11 = measured width (um)
        curveFinder = conditionData_trim(:,6);    % col 6 = curveFinder, ID of full cell cycles
        
        
        % 10. identify unique unique cell cycles
        curveIDs = unique(curveFinder);
        
        
        % 11. for each unique, store width and length at birth and divsion
        W_division = nan(length(curveIDs),1);
        W_birth = nan(length(curveIDs),1);
        
        L_division = nan(length(curveIDs),1);
        L_birth = nan(length(curveIDs),1);
        
        for cc = 1:length(curveIDs)
            
            % widths at birth and division
            currentWidths = widths(curveFinder == curveIDs(cc));
            W_division(cc,1) = currentWidths(end);
            W_birth(cc,1) = currentWidths(1);
            
            % lengths at birth and division
            currentLengths = lengthVals(curveFinder == curveIDs(cc));
            L_division(cc,1) = currentLengths(end);
            L_birth(cc,1) = currentLengths(1);
            
        end
        
        
        % 11. initialize plotting parameters
        color = rgb(palette(condition));
        
        if condition == 1 && timescale == 300
            xmark = shapes{2};
        elseif condition == 1 && timescale == 900
            xmark = shapes{3};
        elseif condition == 1 && timescale == 3600
            xmark = shapes{4};
        else
            xmark = shapes{1};
        end
        
        
        % 12. calculate mean and standard deviation of:
        
        % (i) width at birth, width at division
        W_birth_means = mean(W_birth);
        W_div_means = mean(W_division);
        W_birth_stds = std(W_birth);
        W_div_stds = std(W_division);

        % (ii) length at birth, length at division
        L_birth_means = mean(L_birth);
        L_div_means = mean(L_division);
        L_birth_stds = std(L_birth);
        L_div_stds = std(L_division);
        

        % 13. plot!
        
        % (i) division width vs. birth width
        figure(1)
        errorbar(W_birth_means,W_div_means,W_div_stds,'Color',color)
        hold on
        plot(W_birth_means,W_div_means,'Marker',xmark,'Color',color)
        hold on
        xlabel('birth width (um)')
        ylabel('division width (um)')
        title('population averages from all experiments')
        axis([1.1 1.5 .8 1.8])
        
        % (ii) division length vs. birth length
        figure(2)
        errorbar(L_birth_means,L_div_means,L_div_stds,'Color',color)
        hold on
        plot(L_birth_means,L_div_means,'Marker',xmark,'Color',color)
        hold on
        xlabel('birth length (um)')
        ylabel('division length (um)')
        title('population averages from all experiments')
        axis([1.5 4.5 2 9])

        
    end
    
    
end



