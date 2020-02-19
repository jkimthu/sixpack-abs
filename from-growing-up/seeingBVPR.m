%%  VISUALIZING BVPR
%
%   For working with instantaneous growth rates after SlidingFits.m
%   
%   Goal: plot mean biovolume production rate (with s.e.m.) in each
%         condition over time
            

%   Strategy:
%
%       0. initialize data
%       0. initialize conditions
%       1. build experiment data matrix
%       2. for each condition, bin and plot data
%               3. isolate data from current condition


% last updated: jen, 2018 Mar 8
% commit: plot bvpr vs time for experiment 2018-01-29, 60 min timescale
 

%% plot average growth rate (bvpr) over time
  

% 0. initialize data
clear;
clc;
date = '2018-01-29';
experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
cd(experimentFolder)

%load('lb-monod-2017-09-26-window5-jiggle-varied.mat','D5','M','T');
load(strcat('lb-fluc-',date,'-window5-width1p4-1p7-jiggle-0p5.mat'),'D5','M','M_va','T');

% 0. initialize conditions: col1 = first xy; col2 = final xy; col3 = time (hr) cutoff
conditions = [1 10; 11 20; 21 30; 31 40];% 41 50; 51 60];
binsPerHour = 30;


% 1. build experiment data matrix
exptData = buildDM(D5,M,M_va,T,1,40);

%%
% 2. for each condition, bin and plot data
for c = 1:length(conditions) %number of conditions
    
    % 3. isolate data from current condition
    conditionData = exptData(exptData(:,28) == c,:);    % col 28 = condVals 
    
    % 4. isolate bvpr and time data
    bvpr = conditionData(:,29);                  % col 29 = biovol production rate (cubic um/hr)
    timestamps = conditionData(:,2)/3600;       % time in seconds converted to hours
    clear conditionData
    
    
    % 5. eliminate zero and negative growth rates
    true_bvpr = bvpr(bvpr >= 0);
    true_time = timestamps(bvpr >= 0);
    
    
    % 6. assign timestamps to time bins
    Bins = ceil(true_time*binsPerHour);          
    
    
    % 7. accumulate growth rates by bin, and calculate mean and std dev
    bvpr_Means = accumarray(Bins,true_bvpr,[],@mean);
    bvpr_STDs = accumarray(Bins,true_bvpr,[],@std);
    bvpr_Counts = accumarray(Bins,true_bvpr,[],@length);
    

    % 8. divide standard dev by square root of tracks per bin
    bvpr_sems = bvpr_STDs./sqrt(bvpr_Counts);

    
    % 9. plot mean and s.e.m. over time
    figure(1)
    errorbar(bvpr_Means(1:end-2),bvpr_sems(1:end-2))
    hold on
    axis([0,binsPerHour*10+1,0,16])
    xlabel('Time (hr)')
    ylabel('biomass prod rate (cubic um/hr)')
    legend('fluc','1/1000 LB','ave', '1/50 LB')
    title(date)

    
end
