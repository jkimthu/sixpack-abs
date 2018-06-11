% figure 22

% Goal: plot time-averaged growth from all experiments as a function of
%       fluctuating timescale. plot values for Monod and Jensen's
%       expectations, as well.
%
%       five output plots:
%
%           1. G_high vs G_low across experiments (part one)

%           2. expectations, plotting raw dV/dt values
%           3. expectations, plotting raw dV/dt values normalized to G_jensens
%           4. expectations, plotting dV/dt values normalized by initial volume
%           5. expectations, plotting dV/dt values normalized by initial
%              volume and G_jensens



% Strategy:
%
%       - do we normalize average each fluctuating data point to average value from stable experiments? 
%       - or do we normalize each experiment with its own stable values?
%
%       1. plot G_high vs G_low: if correlated, then suggestive of important day-to-day variability
%                                if not correlated, average of all experiments is fine
%
%       2. plot G_fluc (normalized as decided from step 1) vs timescale,
%          including G_monod and G_jensen's for comparison to expectable
%          values



%  Last edit: Jen Nguyen, 2018 Jun 5
%  Commit: compare G_fluc to G_monod and G_jensens with raw and initial vol
%  normalized plots, using all four replicates in 15 min timescale



% OK let's go!


%% ONE. correlation between inter-experiment G_lows and G_highs?

% Strategy:

%       0. initialize meta data and dV/dt stats
%       0. initialize summary vectors for data
%       1. for each experiment, collect values of G_low and G_high
%       2. plot G_data 


% 0. initialize meta data and dV/dt stats
clear
clc

cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
load('dVdtData_fullOnly_newdVdt.mat') % this data uses instantaneous points present after 3hrs, from full curves only


% 0. initialize summary vectors for data
dataIndex = find(~cellfun(@isempty,storedMetaData));
experimentCount = length(dataIndex);


% 1. for each experiment, collect values of G_low and G_high
counter = 0;
for e = 1:experimentCount
    
    % identify experiment by date
    index = dataIndex(e);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    
    % exclude outlier from analysis
    if strcmp(date, '2017-10-31') == 1 || strcmp (timescale, 'monod') == 1
        disp(strcat(date,': excluded from analysis'))
        continue
    end
    disp(strcat(date, ': collecting data!'))
    counter = counter + 1;
    
    % collect data
    low = 2;
    high = 4;
    
    G_low(counter) = dVdtData_fullOnly_newdVdt{index}{1,low}.mean;
    G_high(counter) = dVdtData_fullOnly_newdVdt{index}{1,high}.mean;
    
    timescales_perG(counter) = timescale;
    dates_perG{counter} = date;
    
end


% 2. plot G_data
figure(1)
plot(G_low*3600,G_high*3600,'o','MarkerSize',10,'LineWidth',2);
%axis([0 5 0 12])
title('test for strong day-to-day correlations in dV/dt')
xlabel('mean dV/dt in low')
ylabel('mean dV/dt in high')


%% TWO. growth expectations: G_fluc vs timescale

% Goal: plot G_fluc (normalized as decided from step 1) vs timescale,
%       including lines for G_monod, G_jensen's, G_high and G_low
%
%       no correlation seen in part one, thus each fluctuating data point
%       will be normalized to average dV/dt from entire stable data set


% Strategy: 
%
%       0. initialize meta data and dV/dt stats
%       0. initialize summary vectors for data
%       1. for each experiment, collect values of G_fluc, G_low, G_ave, G_high
%       2. calculate mean G_low, mean G_ave, mean G_high and G_jensens
%       3. generate hypothetical points, where on the timescale vector:
%             -  0 represents infinitely fast
%             -  5 represents infinitely slow
%       4. plot G_data by timescale


% 0. initialize meta data and dV/dt stats
clear
clc

cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
load('dVdtData_fullOnly_newdVdt.mat') % this data uses instantaneous points present after 3hrs, from full curves only
                                      % units are in cubic um per sec
load('dVdtData_fullOnly_normalized_newdVdt.mat') % both datasets are generated and saved in figure21.m script


% 0. initialize summary vectors for data
dataIndex = find(~cellfun(@isempty,storedMetaData));
experimentCount = length(dataIndex);


% 1. for each experiment, collect values of G_low and G_high
counter = 0;
for e = 1:experimentCount
    
    % identify experiment by date
    index = dataIndex(e);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    
    % exclude outliers from analysis
    if strcmp(date, '2017-10-31') == 1 || strcmp (timescale, 'monod') == 1
        disp(strcat(date,': excluded from analysis'))
        continue
    end
    
%     if strcmp(date, '2018-01-12') == 1
%         disp(strcat(date,': excluded from analysis'))
%         continue
%     end
%     
%     if strcmp(date, '2017-11-13') == 1
%         disp(strcat(date,': excluded from analysis'))
%         continue
%     end


    disp(strcat(date, ': collecting data!'))
    counter = counter + 1;
    
    % collect data
    fluc = 1;
    low = 2;
    ave = 3;
    high = 4;
    
    % raw data
    flucRates_raw(counter,fluc) = dVdtData_fullOnly_newdVdt{index}{1,fluc}.mean*3600;
    stableRates_raw(counter,low) = dVdtData_fullOnly_newdVdt{index}{1,low}.mean*3600;
    stableRates_raw(counter,ave) = dVdtData_fullOnly_newdVdt{index}{1,ave}.mean*3600;
    stableRates_raw(counter,high) = dVdtData_fullOnly_newdVdt{index}{1,high}.mean*3600;
    
    % normalized data
    flucRates_norm(counter,fluc) = dVdtData_fullOnly_normalized_newdVdt{index}{1,fluc}.mean*3600;
    stableRates_norm(counter,low) = dVdtData_fullOnly_normalized_newdVdt{index}{1,low}.mean*3600;
    stableRates_norm(counter,ave) = dVdtData_fullOnly_normalized_newdVdt{index}{1,ave}.mean*3600;
    stableRates_norm(counter,high) = dVdtData_fullOnly_normalized_newdVdt{index}{1,high}.mean*3600;
    
    timescales_perG(counter) = timescale;
    dates_perG{counter} = date;
    
end


% 2. calculate mean G_low, mean G_ave, mean G_high and G_jensens
stableRates_raw_mean = nanmean(stableRates_raw);
stableRates_raw_std = nanstd(stableRates_raw);
G_jensens_raw = (stableRates_raw_mean(high) + stableRates_raw_mean(low))/2;
G_monod_raw = stableRates_raw_mean(ave);

stableRates_norm_mean = nanmean(stableRates_norm);
stableRates_norm_std = nanstd(stableRates_norm);
G_jensens_norm = (stableRates_norm_mean(high) + stableRates_norm_mean(low))/2;
G_monod_norm = stableRates_norm_mean(ave);


% 3. calculate mean dV/dt for each fluctuating timescale
t30 = 1;
t300 = 2;
t900 = 3;
t3600 = 4;

Gfluc_raw_means(t30) = mean(flucRates_raw(timescales_perG==30));
Gfluc_raw_means(t300) = mean(flucRates_raw(timescales_perG==300));
Gfluc_raw_means(t900) = mean(flucRates_raw(timescales_perG==900));
Gfluc_raw_means(t3600) = mean(flucRates_raw(timescales_perG==3600));

Gfluc_raw_std(t30) = std(flucRates_raw(timescales_perG==30));
Gfluc_raw_std(t300) = std(flucRates_raw(timescales_perG==300));
Gfluc_raw_std(t900) = std(flucRates_raw(timescales_perG==900));
Gfluc_raw_std(t3600) = std(flucRates_raw(timescales_perG==3600));

Gfluc_norm_means(t30) = mean(flucRates_norm(timescales_perG==30));
Gfluc_norm_means(t300) = mean(flucRates_norm(timescales_perG==300));
Gfluc_norm_means(t900) = mean(flucRates_norm(timescales_perG==900));
Gfluc_norm_means(t3600) = mean(flucRates_norm(timescales_perG==3600));

Gfluc_norm_std(t30) = std(flucRates_norm(timescales_perG==30));
Gfluc_norm_std(t300) = std(flucRates_norm(timescales_perG==300));
Gfluc_norm_std(t900) = std(flucRates_norm(timescales_perG==900));
Gfluc_norm_std(t3600) = std(flucRates_norm(timescales_perG==3600));


% 4. plot G_data by timescale
% raw values
figure(2)
plot([1 2 3 4],Gfluc_raw_means,'o','Color',rgb('DarkTurquoise'),'MarkerSize',10,'LineWidth',2);
hold on
errorbar([1 2 3 4],Gfluc_raw_means,Gfluc_raw_std,'Color',rgb('DarkTurquoise'));
hold on
plot(-1, G_monod_raw,'o','Color',rgb('DarkCyan'),'MarkerSize',10,'LineWidth',2);
hold on
errorbar(-1,G_monod_raw,stableRates_raw_std(ave),'Color',rgb('DarkCyan'),'LineWidth',2);
hold on
plot(6, G_jensens_raw,'o','Color',rgb('SlateGray'),'MarkerSize',10,'LineWidth',2)
hold on
axis([-1 6 0 12])
title('growth expectations')
xlabel('fluctuating timescale')
ylabel('mean dV/dt (cubic um/hr)')


% raw, normalized by G_jensens
figure(3)
plot([1 2 3 4],Gfluc_raw_means./G_jensens_raw,'o','Color',rgb('DarkTurquoise'),'MarkerSize',10,'LineWidth',2);
hold on
errorbar([1 2 3 4],Gfluc_raw_means./G_jensens_raw,Gfluc_raw_std./G_jensens_raw,'Color',rgb('DarkTurquoise'));
hold on
plot(-1, G_monod_raw/G_jensens_raw,'o','Color',rgb('DarkCyan'),'MarkerSize',10,'LineWidth',2);
hold on
errorbar(-1,G_monod_raw/G_jensens_raw,stableRates_raw_std(ave)./G_jensens_raw,'Color',rgb('DarkCyan'),'LineWidth',2);
hold on
plot(6, G_jensens_raw/G_jensens_raw,'o','Color',rgb('SlateGray'),'MarkerSize',10,'LineWidth',2)

axis([-1 6 0 1.2])
title('growth, relative to Jensens expectations')
xlabel('fluctuating timescale')
ylabel('mean dV/dt, normalized to G_jensens')


% dVdt normalized by initial vol
figure(4)
plot([1 2 3 4],Gfluc_norm_means,'o','Color',rgb('DarkTurquoise'),'MarkerSize',10,'LineWidth',2);
hold on
errorbar([1 2 3 4],Gfluc_norm_means,Gfluc_norm_std,'Color',rgb('DarkTurquoise'));
hold on
plot(-1, G_monod_norm,'o','Color',rgb('DarkCyan'),'MarkerSize',10,'LineWidth',2);
hold on
errorbar(-1,G_monod_norm,stableRates_norm_std(ave),'Color',rgb('DarkCyan'),'LineWidth',2);
hold on
plot(6, G_jensens_norm,'o','Color',rgb('SlateGray'),'MarkerSize',10,'LineWidth',2)
hold on
axis([-1 6 0 2])
title('growth expectations')
xlabel('fluctuating timescale')
ylabel('mean dV/dt / V (1/hr)')


% dVdt normalized by initial vol, normalized by G_jensens
figure(5)
plot([1 2 3 4],Gfluc_norm_means./G_jensens_norm,'o','Color',rgb('DarkTurquoise'),'MarkerSize',10,'LineWidth',2);
hold on
errorbar([1 2 3 4],Gfluc_norm_means./G_jensens_norm,Gfluc_norm_std./G_jensens_norm,'Color',rgb('DarkTurquoise'));
hold on
plot(-1,G_monod_norm/G_jensens_norm,'o','Color',rgb('DarkCyan'),'MarkerSize',10,'LineWidth',2);
hold on
errorbar(-1,G_monod_norm/G_jensens_norm,stableRates_norm_std(ave)./G_jensens_norm,'Color',rgb('DarkCyan'),'LineWidth',2);
hold on
plot(6, G_jensens_norm/G_jensens_norm,'o','Color',rgb('SlateGray'),'MarkerSize',10,'LineWidth',2)

axis([-1 6 0 1.6])
title('growth, relative to Jensens expectations')
xlabel('fluctuating timescale')
ylabel('mean dV/dt / V, normalized to G_jensens')

