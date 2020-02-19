%%  VISUALIZING MU
%
%   For working with instantaneous growth rates after SlidingFits.m
%   
%   Goals:
%
%       1. Plot raw data and mu for selected trajectories in each condition
%               - Any qualitative differences?
%
%       2. Plot average mu (per cell cycle) vs. cell cycle #
%               - Does a steady-state emerge?
%
%       3. Average mu (with standard deviation) per timepoint in each condition 
%               - Another way of looking at steady state
%
% 

% last updated: jen, 2018 Mar 8
% commit: edit to automatically change directory into applicable experiment folder
 

%% CHECK FOUR: plot average growth rate over time
%
%       - generates a single plot with all conditions 
%       - options to plot standard deviation or standard error
%       - saves average mu, standard dev, s.e.m., and number of tracks per bin per condition 
%       

% Initialize
clear;
clc;
date = '2018-02-01';
experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
cd(experimentFolder)

%load('lb-monod-2017-09-26-window5-jiggle-varied.mat','D5','M','T');
load(strcat('lb-fluc-',date,'-window5-width1p4-1p7-jiggle-0p5.mat'),'D5','M','M_va','T');

% defining conditions: col1 = first xy; col2 = final xy; col3 = time (hr) cutoff
conditions = [1 10; 11 20; 21 30; 31 40];% 41 50; 51 60];
binsPerHour = 30;

%%
for i = 1:length(conditions) %number of conditions
    
    %  initialize data concatenation
    mu_l = [];
    mu_va = [];
    Time_cond = [];
    
    
    for n = conditions(i,1):conditions(i,2)
        for m = 1:length(M_va{n})
            
            %  assemble all instantaneous growth rates into a single vector
            mu_l = [mu_l; M{n}(m).mu];
            mu_va = [mu_va; M_va{n}(m).mu_va];

            
            %  assemble a corresponding timestamp vector
            vectorLength = length(M_va{n}(m).mu_va);
            trackFrames = D5{n}(m).Frame(3:vectorLength+2);
            Time_cond = [Time_cond; T{n}(trackFrames)];
            
        end
    end
    
    %  convert all timestamps from seconds to hours
    Time_cond = Time_cond/3600;
    
    %  eliminate zero and negative growth rates
    true_mu_l = mu_l(mu_va>=0);
    true_mu_va = mu_va(mu_va>=0);
    true_time = Time_cond(mu_va>=0);
    
    %  assign timestamps to time bins
    Bins = ceil(true_time*binsPerHour);            % multiplying by 200 gives time bins of 0.005 hr
    
    
    %  accumulate growth rates by bin, and calculate mean and std dev
    mu_l_Means = accumarray(Bins,true_mu_l,[],@mean);
    mu_l_STDs = accumarray(Bins,true_mu_l,[],@std);
    mu_l_Counts = accumarray(Bins,true_mu_l,[],@length);
    
    mu_va_Means = accumarray(Bins,true_mu_va,[],@mean);
    mu_va_STDs = accumarray(Bins,true_mu_va,[],@std);
    mu_va_Counts = accumarray(Bins,true_mu_va,[],@length);
    

    %   2. divide standard dev by square root of tracks per bin
    mu_l_sems = mu_l_STDs./sqrt(mu_l_Counts);
    mu_va_sems = mu_va_STDs./sqrt(mu_va_Counts);

    figure(1)
    errorbar(mu_l_Means,mu_l_sems)
    hold on
    axis([0,binsPerHour*10+1,0,4])
    xlabel('Time (hr)')
    ylabel('doubling rate of length (1/hr)')
    legend('fluc','1/1000 LB','ave', '1/50 LB')
    title(date)
    
    figure(2)
    errorbar(mu_va_Means,mu_va_sems)
    hold on
    axis([0,binsPerHour*10+1,0,4])
    xlabel('Time (hr)')
    ylabel('doubling rate of volume (1/hr)')
    legend('fluc','1/1000 LB','ave', '1/50 LB')
    title(date)


    
end
