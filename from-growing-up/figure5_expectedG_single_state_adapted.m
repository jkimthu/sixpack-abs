%% figure 6 - measured G_fluc is not predicted by steady-state null models


% Output: plot of time-averaged growth rate (G) vs nutrient timescale

% Input: (1) data structure from figure 3A, growthRates_monod_curve,
%        (2) meta data structure of all experiments
%        (3) growth rate vs time data from single up- and downshifts


%        (1) fluctuating experiment data
%        contains summary stats from each condition of each experiment
%
%        (2) meta data
%        psuedo-manually compiled structure, storedMetaData.mat
%        structure helps ensure corrent sorting of growthRates_monod_curve data
%
%        (3) single shift data
%        upshift data: average signal between two replicates
%        downshift data: only one replicate reaches G_low, uses this signal
%        uses signal until calculated stabililzation time from figure 5,
%        after this stabilization, use G_low or G_high to simulate growth
%        rate for the rest of the low or high nutrient phase


%        this code is written as if all experiment files (.mat, containing D5
%        and T) and meta data file were in the same folder



% Strategy: 

%  A. Initialize 
%
%       0. input complete meta data
%       0. input fluctuating experiment data
%       0. input single shift experiment data
%       0. define period lengths for slow calculations of predicted G (12,24,48,96 hours)
%       0. define period lengths for fast calculations of predicted G (30 sec,5,15,60 min)



%  B. Assemble fluctuating data according to timescale or boundary conditions
%
%       1. for each experiment, collect time-averaged growth rates for each condition
%       2. calculate mean and standard dev G_low, G_ave, G_high and G_jensens across all data
%       3. calculate mean and standard dev G_fluc for each fluctuating timescale
%       4. plot G_data by timescale, with the following normalizations:
%               i. none
%              ii. normalized by G_monod
%             iii. normalized by G_jensens
%       5. display raw values of G_ave, G_jensens, G_fluc, G_high and G_low



%  C. Slow calculations of predicted G from single shift data
%
%       0. initialize time to stabilize for both shifts
%          note: these manually set values are calculated means from Figure 5
%                (response quantifications)
%       1. trim growth rate signals to times until stabilized
%       2. trim growth signals to times postshift
%       3. time-average transitions
%       4. calculate time-averaged growth rates during high nutrient phase
%       5. calculate time-averaged growth rates during low nutrient phase
%       6. calculate time-averaged growth rates from entire period



%  D. Add slow predictions to timescale plot


%  E. Fast calculations of predicted G from single shift data
%
%     Method 1: low-adapted (using single upshift data)
%     Method 2: high-adapted


%  F. Add fast predictions to timescale plot


% Last edit: jen, 2019 November 27
% Commit: compares measured G_fluc to null models that predict G_fluc from steady states


% OK let's go!

%% A. Initialize


clear
clc


% 0. input complete meta data
%cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
exptArray = [2:4,5:7,9,10,11,12,13,14,15]; % list experiments by index


% 0. input fluctuating experiment data
%cd('/Users/jen/Documents/StockerLab/Writing/manuscript 1/figures v15')
load('growthRates_monod_curve.mat')
dataIndex = find(~cellfun(@isempty,growthRates_monod_curve));


% 0. input single shift experiment data
load('response_singleUpshift.mat')
load('response_singleDownshift.mat')



% 0. define period lengths for hypothetical calculations (period in hours)
periods_slow = [12,24,48,96]; % hours
periods_fast = [30,300,900,3600]; % sec


%% B. Assemble fluctuating data according to timescale or boundary conditions

% 1. for each experiment, collect time-averaged growth rates for each condition
counter = 0;
for e = 1:length(exptArray)
    
    % identify experiment by date
    index = exptArray(e);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    
    disp(strcat(date, ': collecting data!'))
    counter = counter + 1;
    
    % collect growth rate data
    fluc = 1;
    low = 2;
    ave = 3;
    high = 4;
    
    flucRates(counter,fluc) = growthRates_monod_curve{index}{1,fluc}.mean;
    stableRates(counter,low) = growthRates_monod_curve{index}{1,low}.mean;
    stableRates(counter,ave) = growthRates_monod_curve{index}{1,ave}.mean;
    stableRates(counter,high) = growthRates_monod_curve{index}{1,high}.mean;
    
    timescales_perG(counter) = timescale;
    dates_perG{counter} = date;
    
end
clear timescale


% 2. calculate mean and standard dev G_low, G_ave, G_high and G_jensens across all data
stableRates_mean = nanmean(stableRates);
stableRates_std = nanstd(stableRates);
G_jensens = (stableRates_mean(high) + stableRates_mean(low))/2;
G_monod = stableRates_mean(ave);



% 3. normalize each G_fluc by it's corresponding G_ave
G_fluc_norm = flucRates./stableRates(:,ave);



% 4. calculate mean and standard dev G_fluc for each fluctuating timescale
t30 = 1;
t300 = 2;
t900 = 3;
t3600 = 4;

Gfluc_means(t30) = mean(G_fluc_norm(timescales_perG==30));
Gfluc_means(t300) = mean(G_fluc_norm(timescales_perG==300));
Gfluc_means(t900) = mean(G_fluc_norm(timescales_perG==900));
Gfluc_means(t3600) = mean(G_fluc_norm(timescales_perG==3600));

Gfluc_std(t30) = std(G_fluc_norm(timescales_perG==30));
Gfluc_std(t300) = std(G_fluc_norm(timescales_perG==300));
Gfluc_std(t900) = std(G_fluc_norm(timescales_perG==900));
Gfluc_std(t3600) = std(G_fluc_norm(timescales_perG==3600));



% 4. plot G_data by timescale

% growth rate normalized
figure(6)
plot([2 3 4 5],Gfluc_means,'o','Color',rgb('DarkTurquoise'),'MarkerSize',10,'LineWidth',2);
hold on
errorbar([2 3 4 5],Gfluc_means, Gfluc_std,'o','Color',rgb('DarkTurquoise'));
hold on
plot(1,G_monod./G_monod,'o','Color',rgb('DarkCyan'),'MarkerSize',10,'LineWidth',2);
hold on
errorbar(1,G_monod./G_monod,stableRates_std(ave)./G_monod,'o','Color',rgb('DarkCyan'),'LineWidth',2);
hold on
plot(10, G_jensens./G_monod,'o','Color',rgb('SlateGray'),'MarkerSize',10,'LineWidth',2)


axis([1 10 0.1 1.4])
title('growth, relative to average nutrient concentration')
xlabel('fluctuating timescale')
ylabel('mean growth rate, normalized by G_ave)') % growth rate = d(logV)/(dt*ln(2))


% 5. display raw values of G_ave, G_jensens, G_fluc, G_high and G_low
%G_monod
%G_jensens
%Gfluc_means
%Gfluc_std
G_high = stableRates_mean(high);
G_low = stableRates_mean(low);


clear counter index
clear e experimentCount t30 t300 t900 t3600
clear date 


%% C. Slow calculations of predicted G from single shift data


% 0. initialize time to stabilize for both shifts
ts_up = 116.3; % min, mean value from Figure 4 (response quantifications)
ts_down = 297.5;


% 0. initialize shorter names for data from both shifts
upshift_means = upshift_means_single;
upshift_times = upshift_times_single;
downshift_means = downshift_means_single;
downshift_times = downshift_times_single;
clear upshift_means_single upshift_times_single downshift_means_single downshift_times_single


% 1. trim growth rate signals to times until stabilized
trimmed_upshift_gr = upshift_means(upshift_times < ts_up);
trimmed_upshift_times = upshift_times(upshift_times < ts_up);

trimmed_downshift_gr = downshift_means(downshift_times < ts_down);
trimmed_downshift_times = downshift_times(downshift_times < ts_down);
clear upshift_means upshift_times downshift_means downshift_times


% 2. trim growth signals to times postshift
transition_upshift = trimmed_upshift_gr(trimmed_upshift_times >= 0);
transition_upshift_times = trimmed_upshift_times(trimmed_upshift_times >= 0);

transition_downshift = trimmed_downshift_gr(trimmed_downshift_times >= 0);
transition_downshift_times = trimmed_downshift_times(trimmed_downshift_times >= 0);
clear trimmed_upshift_gr trimmed_upshift_times trimmed_downshift_gr trimmed_downshift_times


% 3. time-average transitions
G_transit_up = mean(transition_upshift);
G_transit_down = mean(transition_downshift);


% 4. calculate time-averaged growth rates during high nutrient phase
for p = 1:length(periods_slow)
    
    tau = periods_slow(p) * 60; % convert hr to min
    phase = tau/2;
    
    fractionTransition = ts_up/phase;
    weightedTransition = G_transit_up * fractionTransition;
    weightedStable = G_high * (1-fractionTransition);
    
    G_phase_high(p) = mean(weightedTransition+weightedStable);
    
end
clear fractionTransition weightedTransition weightedStable tau phase p 



% 5. calculate time-averaged growth rates during low nutrient phase
for p = 1:length(periods_slow)
    
    tau = periods_slow(p) * 60; % convert hr to min
    phase = tau/2;
    
    fractionTransition = ts_down/phase;
    weightedTransition = G_transit_down * fractionTransition;
    weightedStable = G_low * (1-fractionTransition);
    
    G_phase_low(p) = mean(weightedTransition+weightedStable);
    
end
clear fractionTransition weightedTransition weightedStable tau phase p
clear high ave low fluc


% 6. calculate time-averaged growth rates from entire period
G_periods = (G_phase_high + G_phase_low)./2;
G_periods_normalized = G_periods./G_monod;
clear G_phase* periods_slow


% 7. plot slow predictions to timescale plot
figure(6)
hold on
plot([6 7 8 9],G_periods_normalized,'o','Color',rgb('SlateGray'),'MarkerSize',10,'LineWidth',2);
hold on
axis([1,10,0.3,1.2])


%% D. Fast calculations of predicted G from single shift data
%
%  1) high-adapted, from single-shifts
%  2) low-adapted, from single-shifts


%  Single-state adapated cells, adapted to either high or low 

%  1. for each timescale, initialize the value of 1/2 period
G_low_adapted_in_high = zeros(length(periods_fast),1);
G_high_adapted_in_low = zeros(length(periods_fast),1);


%  2. because first single-downshift and upshift measurement is 7.5 min after the shift,
%     approximate a linear slope between tpt #2 and t=0 value in case of upshifts and
%     growth rate = 0 in case of downshifts
first_tpt_up = transition_upshift_times(2);
first_tpt_down = transition_downshift_times(2);
first_gr_up = transition_upshift(2);
first_gr_down = transition_downshift(2);
t_zero_gr_up = transition_upshift(1);
t_zero_gr_down = 0;

slope_up = (first_gr_up - t_zero_gr_up)/first_tpt_up;
slope_down = (first_gr_down - t_zero_gr_down)/first_tpt_down;
resolved_time_up = linspace(0,first_tpt_up,(first_tpt_up/0.25)+1);
resolved_time_down = linspace(0,first_tpt_down,(first_tpt_down/0.25)+1);

resolved_upshift = slope_up*resolved_time_up + t_zero_gr_up;
resolved_downshift = slope_down*resolved_time_down + t_zero_gr_down;


%  3. for each timescale, 
for tscale = 1:length(periods_fast)
    
    timescale = periods_fast(tscale);
    half_T = timescale/2;
    half_T_minutes = half_T/60; % convert sec to minutes
    
    %  i) calculate mean growth rate across half period from single upshift data
    if half_T_minutes <= resolved_time_up(end) % first_tpt_up
        
        tpt = find(resolved_time_up == half_T_minutes);
        G_low_adapted_in_high(tscale) = mean(resolved_upshift(1:tpt));
        
    else
        % when using later data, the time between timepoints are greater
        % and should thus have greater weight in the overall mean
        resolved_fraction = resolved_time_up(end)/half_T_minutes;
        
        % calculate mean growth rate after resolved time
        extended_upshift_data = transition_upshift(3:end);
        extended_upshift_time = transition_upshift_times(3:end);
        
        tpt = find(extended_upshift_time == half_T_minutes);
        extended_mean_gr = mean(extended_upshift_data(1:tpt));
        
        G_low_adapted_in_high(tscale) = (mean(resolved_upshift)*resolved_fraction) + (extended_mean_gr*(1-resolved_fraction));
    end
    clear tpt
    
    
    % ii) calculate mean growth rate across half period from single downshift data
    if half_T_minutes <= resolved_time_down(end)
        
        tpt = find(resolved_time_down == half_T_minutes);
        G_high_adapted_in_low(tscale) = mean(resolved_downshift(1:tpt));
        
    else
        % when using later data, the time between timepoints are greater
        % and should thus have greater weight in the overall mean
        resolved_fraction = resolved_time_down(end)/half_T_minutes;
        
        % calculate mean growth rate after resolved time
        extended_downshift_data = transition_downshift(3:end);
        extended_downshift_time = transition_downshift_times(3:end);
        
        tpt = find(extended_downshift_time == half_T_minutes);
        extended_mean_gr = mean(extended_downshift_data(1:tpt));
        
        G_high_adapted_in_low(tscale) = (mean(resolved_downshift)*resolved_fraction) + (extended_mean_gr*(1-resolved_fraction));
        
    end
    clear tpt
end
clear tscale


    
    
% 4. combine high and low phases to calculate G_restart
%    G_restart = ( G_in_high(half period) + G_in_low(half period) )/ 2
G_adapted_low = (G_low_adapted_in_high + G_low)./2;
G_adapted_low_normalized = G_adapted_low./G_monod;
    
G_adapted_high = (G_high_adapted_in_low + G_high)./2;
G_adapted_high_normalized = G_adapted_high./G_monod;



% 5. plot both fast prediction onto timescale plot
figure(6)
hold on

% Method 1. low adapted (plot like slow predictions)
plot([2 3 4 5],G_adapted_low_normalized,'o','Color',rgb('SlateGray'),'MarkerSize',10,'LineWidth',2);
hold on

% Method 2. high adapted 
plot([2 3 4 5],G_adapted_high_normalized,'o','Color',rgb('DarkSlateGray'),'MarkerSize',10,'LineWidth',2);
axis([1,10,0.3,1.2])


%% E. Caluclate difference from measured fluc and predicted fluc

% 1. predicted by low-adapted simulations

% i) create vector of predicted values that matches replicate measured vals
G_fluc_predicted_low(1:3,1) = G_adapted_low_normalized(1);
G_fluc_predicted_low(4:6,1) = G_adapted_low_normalized(2);
G_fluc_predicted_low(7:10,1) = G_adapted_low_normalized(3);
G_fluc_predicted_low(11:13,1) = G_adapted_low_normalized(4);

% ii) calculate relative change from measured to predicted
relative_G_low = (G_fluc_norm-G_fluc_predicted_low)./G_fluc_norm;

% iii) calculate mean and std for each timescale
relative_G_low_mean(1) = mean(relative_G_low(1:3));
relative_G_low_mean(2) = mean(relative_G_low(4:6));
relative_G_low_mean(3) = mean(relative_G_low(7:10));
relative_G_low_mean(4) = mean(relative_G_low(11:13))

relative_G_low_std(1) = std(relative_G_low(1:3));
relative_G_low_std(2) = std(relative_G_low(4:6));
relative_G_low_std(3) = std(relative_G_low(7:10));
relative_G_low_std(4) = std(relative_G_low(11:13));



% 2. predicted by high-adapted simulations

% i) create vector of predicted values that matches replicate measured vals
G_fluc_predicted_high(1:3,1) = G_adapted_high_normalized(1);
G_fluc_predicted_high(4:6,1) = G_adapted_high_normalized(2);
G_fluc_predicted_high(7:10,1) = G_adapted_high_normalized(3);
G_fluc_predicted_high(11:13,1) = G_adapted_high_normalized(4);

% ii) calculate relative change from measured to predicted
relative_G_high = (G_fluc_norm-G_fluc_predicted_high)./G_fluc_norm;

% iii) calculate mean and std for each timescale
relative_G_high_mean(1) = mean(relative_G_high(1:3));
relative_G_high_mean(2) = mean(relative_G_high(4:6));
relative_G_high_mean(3) = mean(relative_G_high(7:10));
relative_G_high_mean(4) = mean(relative_G_high(11:13));

relative_G_high_std(1) = std(relative_G_high(1:3));
relative_G_high_std(2) = std(relative_G_high(4:6));
relative_G_high_std(3) = std(relative_G_high(7:10));
relative_G_high_std(4) = std(relative_G_high(11:13));

