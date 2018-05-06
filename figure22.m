% figure 22

% Goal: plot time-averaged growth from all experiments as a function of
%       fluctuating timescale. draw lines for Monod and Jensen's
%       expectations, as well as measured high and low.


% Strategy:
%
%       - do we normalize average each fluctuating data point to average value from stable experiments? 
%       - or do we normalize each experiment with its own stable values?
%
%       1. plot G_high vs G_low: if correlated, then suggestive of important day-to-day variability
%                                if not correlated, average of all experiments is fine
%
%       2. plot G_fluc (normalized as decided from step 1) vs timescale,
%          including lines for G_monod, G_jensen's, G_high and G_low




%  Last edit: Jen Nguyen, 2018 May 6
%  Commit: plot growth under fluctuations, raw and normalized to Jensens
%          expectations



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
    fluc = 1;
    low = 2;
    ave = 3;
    high = 4;
    
    flucRates(counter,fluc) = dVdtData_fullOnly_newdVdt{index}{1,fluc}.mean*3600;
    stableRates(counter,low) = dVdtData_fullOnly_newdVdt{index}{1,low}.mean*3600;
    stableRates(counter,ave) = dVdtData_fullOnly_newdVdt{index}{1,ave}.mean*3600;
    stableRates(counter,high) = dVdtData_fullOnly_newdVdt{index}{1,high}.mean*3600;
    
    timescales_perG(counter) = timescale;
    dates_perG{counter} = date;
    
end


% 2. calculate mean G_low, mean G_ave, mean G_high and G_jensens
stableRates_mean = nanmean(stableRates);
stableRates_std = nanstd(stableRates);
G_jensens = (stableRates_mean(high) + stableRates_mean(low))/2;
G_monod = stableRates_mean(ave);


% 3. calculate mean dV/dt for each fluctuating timescale
t30 = 1;
t300 = 2;
t900 = 3;
t3600 = 4;

Gfluc_means(t30) = mean(flucRates(timescales_perG==30));
Gfluc_means(t300) = mean(flucRates(timescales_perG==300));
Gfluc_means(t900) = mean(flucRates(timescales_perG==900));
Gfluc_means(t3600) = mean(flucRates(timescales_perG==3600));

Gfluc_std(t30) = std(flucRates(timescales_perG==30));
Gfluc_std(t300) = std(flucRates(timescales_perG==300));
Gfluc_std(t900) = std(flucRates(timescales_perG==900));
Gfluc_std(t3600) = std(flucRates(timescales_perG==3600));

% 4. plot G_data by timescale
% raw values
figure(2)
plot([1 2 3 4],Gfluc_means,'o','Color',rgb('DarkTurquoise'),'MarkerSize',10,'LineWidth',2);
hold on
errorbar([1 2 3 4],Gfluc_means,Gfluc_std,'Color',rgb('DarkTurquoise'));
hold on
plot(-1, G_monod,'o','Color',rgb('SlateGray'),'MarkerSize',10,'LineWidth',2);
hold on
errorbar(-1,G_monod,stableRates_std(ave),'Color',rgb('SlateGray'),'LineWidth',2);
hold on
plot(6, G_jensens,'o','Color',rgb('SlateGray'),'MarkerSize',10,'LineWidth',2)
axis([-1 6 0 12])
title('growth expectations')
xlabel('fluctuating timescale')
ylabel('mean dV/dt (cubic um/hr)')

% normalized by G_jensens
figure(3)
plot([1 2 3 4],Gfluc_means./G_jensens,'o','Color',rgb('DarkTurquoise'),'MarkerSize',10,'LineWidth',2);
hold on
errorbar([1 2 3 4],Gfluc_means./G_jensens,Gfluc_std./G_jensens,'Color',rgb('DarkTurquoise'));
hold on
plot(-1, G_monod/G_jensens,'o','Color',rgb('SlateGray'),'MarkerSize',10,'LineWidth',2);
hold on
errorbar(-1,G_monod/G_jensens,stableRates_std(ave)./G_jensens,'Color',rgb('SlateGray'),'LineWidth',2);
hold on
plot(6, G_jensens/G_jensens,'o','Color',rgb('SlateGray'),'MarkerSize',10,'LineWidth',2)
axis([-1 6 0 1.2])
title('growth, relative to Jensens expectations')
xlabel('fluctuating timescale')
ylabel('mean dV/dt, normalized to G_jensens')

