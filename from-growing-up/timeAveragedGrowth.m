% timeAveragedGrowth

% goal: plot time-averaged growth from all experiments as a function of
%       fluctuating timescale. draw lines for Monod and Jensen's
%       expectations, as well as measured high and low.


%  Last edit: Jen Nguyen, 2018 Mar 6

%  Commit:  plot growth rate in fluc conditions alongside Monod and Jensen's based expectations,
%           (1) raw values, (2) raw G_high vs G_low
%           (3) values normalized by measured G_monod (ave stable)


% for strategies, see each section


%% ONE (A). G_jensen and G_fluc NOT normalized by G_monod (BVPR)


% Strategy:

%       0. initialize meta data and bvpr stats
%       0. initialize summary vectors for data
%       1. for each experiment, collect values of G_fluc, G_low, G_ave, G_high
%       2. calculate G_jensens
%       3. generate hypothetical points, where on the timescale vector:
%             -  0 represents infinitely fast
%             -  5 represents infinitely slow
%       4. plot G_data by timescale, varying brightness bewteen replicates



% 0. initialize meta data and bvpr stats

clear
clc

cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
load('bioProdRateData.mat')

% 0. initialize summary vectors for data
dataIndex = find(~cellfun(@isempty,storedMetaData));
experimentCount = length(dataIndex);


% 1. for each experiment, collect values of G_fluc, G_low, G_ave, G_high
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
    
    G_fluc(counter) = bioProdRateData{index}{1,fluc}.mean;
    G_low(counter) = bioProdRateData{index}{1,low}.mean;
    G_ave(counter) = bioProdRateData{index}{1,ave}.mean;
    G_high(counter) = bioProdRateData{index}{1,high}.mean;
    
    timescales_perG(counter) = timescale;
    dates_perG{counter} = date;
    
end
clear fluc low ave high counter date timescale index e


% 2. calculate G_jensens
G_jensens = (G_high+G_low)/2;


% 3. generate hypothetical points, where on the timescale vector:
%       0 represents infinitely fast
%       5 represents infinitely slow

% timescale
timescale_vector(timescales_perG==30) = 1;
timescale_vector(timescales_perG==300) = 2;
timescale_vector(timescales_perG==900) = 3;
timescale_vector(timescales_perG==3600) = 4;

% G_fluctuating
gf_30 = G_fluc(timescales_perG==30);
gf_300 = G_fluc(timescales_perG==300);
gf_900 = G_fluc(timescales_perG==900);
gf_3600 = G_fluc(timescales_perG==3600);

% hypothetical infinutes
hypothetical_0 = mean(G_ave);
hypothetical_30 = mean(gf_30);
hypothetical_300 = mean(gf_300);
hypothetical_900 = mean(gf_900);
hypothetical_3600 = mean(gf_3600);
hypothetical_5 = nanmean(G_jensens);

hypothetical_vector = [hypothetical_0 hypothetical_30 hypothetical_300 hypothetical_900 hypothetical_3600 hypothetical_5];
hypothetical_time = [0 1 2 3 4 5];


% 4. plot G_data by timescale, varying brightness bewteen replicates

figure(1)
color_Monod = rgb('DarkBlue');
color_Jensen = rgb('Teal');
color_fluc = rgb('Chocolate');
color_hyp = rgb('DimGray');

for i = 1:length(G_jensens)
    
    % plot G_monod
    plot(timescale_vector,G_ave,'o','Color',color_Monod,'MarkerSize',8,'LineWidth',2); % dark blue
    hold on
    
    % plot G_jensen
    plot(timescale_vector,G_jensens,'o','Color',color_Jensen,'MarkerSize',8,'LineWidth',2); % teal
    hold on
    
    % plot G_fluc
    plot(timescale_vector,G_fluc,'o','Color',color_fluc,'MarkerSize',8,'LineWidth',2); % chocolate
    hold on
    
    % plot hypothetical line
    plot(hypothetical_time, hypothetical_vector, '-.', 'Color',color_fluc,'LineWidth',2); % chocolate
    plot(hypothetical_time, [hypothetical_0 NaN NaN NaN NaN hypothetical_5], 'o', 'Color',color_hyp,'LineWidth',2); % dim gray
    
    % plot mean monod and jensens
    plot(hypothetical_time, ones(1,length(hypothetical_time))*hypothetical_0, ':', 'Color',color_Monod,'LineWidth',2); % dark blue
    plot(hypothetical_time, ones(1,length(hypothetical_time))*hypothetical_5, ':', 'Color',color_Jensen,'LineWidth',2); % teal
    
    
end
axis([0 5 0 12])
title('expected and measured mean growth rates')
xlabel('timescale of nutrient fluctuation')
ylabel('mean bvpr (cubic um/hr)')


%% ONE (B). G_jensen and G_fluc NOT normalized by G_monod (MU)


% Strategy:

%       0. initialize meta data and mu stats
%       0. initialize summary vectors for data
%       1. for each experiment, collect values of G_fluc, G_low, G_ave, G_high
%       2. calculate G_jensens
%       3. generate hypothetical points, where on the timescale vector:
%             -  0 represents infinitely fast
%             -  5 represents infinitely slow
%       4. plot G_data by timescale, varying brightness bewteen replicates



% 0. initialize meta data and mu stats

clear
clc

cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
load('muData.mat')

% 0. initialize summary vectors for data
dataIndex = find(~cellfun(@isempty,storedMetaData));
experimentCount = length(dataIndex);


% 1. for each experiment, collect values of G_fluc, G_low, G_ave, G_high
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
    
    G_fluc(counter) = muData{index}{1,fluc}.mean;
    G_low(counter) = muData{index}{1,low}.mean;
    G_ave(counter) = muData{index}{1,ave}.mean;
    G_high(counter) = muData{index}{1,high}.mean;
    
    timescales_perG(counter) = timescale;
    dates_perG{counter} = date;
    
end
clear fluc low ave high counter date timescale index e


% 2. calculate G_jensens
G_jensens = (G_high+G_low)/2;


% 3. generate hypothetical points, where on the timescale vector:
%       0 represents infinitely fast
%       5 represents infinitely slow

% timescale
timescale_vector(timescales_perG==30) = 1;
timescale_vector(timescales_perG==300) = 2;
timescale_vector(timescales_perG==900) = 3;
timescale_vector(timescales_perG==3600) = 4;

% G_fluctuating
gf_30 = G_fluc(timescales_perG==30);
gf_300 = G_fluc(timescales_perG==300);
gf_900 = G_fluc(timescales_perG==900);
gf_3600 = G_fluc(timescales_perG==3600);

% hypothetical infinutes
hypothetical_0 = mean(G_ave);
hypothetical_30 = mean(gf_30);
hypothetical_300 = mean(gf_300);
hypothetical_900 = mean(gf_900);
hypothetical_3600 = mean(gf_3600);
hypothetical_5 = nanmean(G_jensens);

hypothetical_vector = [hypothetical_0 hypothetical_30 hypothetical_300 hypothetical_900 hypothetical_3600 hypothetical_5];
hypothetical_time = [0 1 2 3 4 5];


% 4. plot G_data by timescale, varying brightness bewteen replicates

figure(2)
color_Monod = rgb('DarkBlue');
color_Jensen = rgb('Teal');
color_fluc = rgb('Chocolate');
color_hyp = rgb('DimGray');

for i = 1:length(G_jensens)
    
    % plot G_monod
    plot(timescale_vector,G_ave,'o','Color',color_Monod,'MarkerSize',8,'LineWidth',2); % dark blue
    hold on
    
    % plot G_jensen
    plot(timescale_vector,G_jensens,'o','Color',color_Jensen,'MarkerSize',8,'LineWidth',2); % teal
    hold on
    
    % plot G_fluc
    plot(timescale_vector,G_fluc,'o','Color',color_fluc,'MarkerSize',8,'LineWidth',2); % chocolate
    hold on
    
    % plot hypothetical line
    plot(hypothetical_time, hypothetical_vector, '-.', 'Color',color_fluc,'LineWidth',2); % chocolate
    plot(hypothetical_time, [hypothetical_0 NaN NaN NaN NaN hypothetical_5], 'o', 'Color',color_hyp,'LineWidth',2); % dim gray
    
    % plot mean monod and jensens
    plot(hypothetical_time, ones(1,length(hypothetical_time))*hypothetical_0, ':', 'Color',color_Monod,'LineWidth',2); % dark blue
    plot(hypothetical_time, ones(1,length(hypothetical_time))*hypothetical_5, ':', 'Color',color_Jensen,'LineWidth',2); % teal
    
    
end
axis([0 5 0 3.5])
title('expected and measured mean growth rates')
xlabel('timescale of nutrient fluctuation')
ylabel('mean mu (cubic um/hr)')


%% TWO (A). Is normalization even necessary? Verdict: NO, not with BVPR.


% assessed with the logic that scaling via normalization makes sense when
% there is a correlation between experiments. i.e. that experiments with
% higher growth rates will have higher values for high, ave, and low. 

% if the converse is true, that no correlation exists, then the variation
% in growth rates observed from the same condition between different
% experiments is simply noise.

% this section plots G_high as a function of G_low, such that each point is
% the mean growth rate in high LB over low LB. in seeing no correlation, our
% variation seems to be a result of noise. (to determine whether it is
% measurement noise or day-to-day variability, variation should then be
% compared for the same condition between days and the same condition
% within the same day.)


% Strategy:

%       0. initialize meta data and bvpr stats
%       0. initialize summary vectors for data
%       1. collect all values of G_high, G_low and s.e.m.s
%       2. plot G_high vs G_low


% 0. initialize meta data and bvpr stats

clear
clc

cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
load('bioProdRateData.mat')

% 0. initialize summary vectors for data
dataIndex = find(~cellfun(@isempty,storedMetaData));
experimentCount = length(dataIndex);

% 1. collect all values of G_high, G_low and s.e.m.s
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
    high = 4;
    low = 2;
    
    G_high_means(counter) = bioProdRateData{index}{1,high}.mean;
    G_high_sems(counter) = bioProdRateData{index}{1,high}.sem;
    G_low(counter) = bioProdRateData{index}{1,low}.mean;
    G_low_sems(counter) = bioProdRateData{index}{1,low}.sem;
    
    timescales_perG(counter) = timescale;
    dates_perG{counter} = date;
    
end

% 2. plot G_high vs G_low
shift = 0.05;

figure(3)
for i = 1:length(timescales_perG)
    if timescales_perG(i) == 30
        plot(G_low(i),G_high_means(i),'o','Color',[1 0 0],'MarkerSize',10); % red
        text(G_low(i)+shift,G_high_means(i)+shift, dates_perG(i),'FontSize',10);
        hold on
    elseif timescales_perG(i) == 300
        plot(G_low(i),G_high_means(i),'o','Color',[1 0.85 0.01],'MarkerSize',10); % sunflower yellow
        text(G_low(i)+shift,G_high_means(i)+shift, dates_perG(i),'FontSize',10);
        hold on
    elseif timescales_perG(i) == 900
        plot(G_low(i),G_high_means(i),'o','Color',[0 0.7 0.7],'MarkerSize',10); % green
        text(G_low(i)+shift,G_high_means(i)+shift, dates_perG(i),'FontSize',10);
        hold on
    elseif timescales_perG(i) == 3600
        plot(G_low(i),G_high_means(i),'o','Color',[0.5 0 0.9],'MarkerSize',10); % purple
        text(G_low(i)+shift,G_high_means(i)+shift, dates_perG(i),'FontSize',10);
        hold on
    end
end
title('no correlation between G_H and G_L')
xlabel('mean bvpr in low (cubic um/hr)')
ylabel('mean bvpr in high (cubic um/hr)')
axis([1.2 3.7 10 19])


%% TWO (B). Is normalization even necessary? Verdict: YES, positive correlation with MU


% assessed with the logic that scaling via normalization makes sense when
% there is a correlation between experiments. i.e. that experiments with
% higher growth rates will have higher values for high, ave, and low. 

% if the converse is true, that no correlation exists, then the variation
% in growth rates observed from the same condition between different
% experiments is simply noise.

% this section plots G_high as a function of G_low, such that each point is
% the mean growth rate in high LB over low LB. in seeing a correlation, our
% variation seems to be a result of day-to-day variability in growth.

% (to determine whether it is measurement noise or day-to-day variability,
% variation should then be compared for the same condition between days
% and the same condition within the same day.)


% Strategy:

%       0. initialize meta data and mu stats
%       0. initialize summary vectors for data
%       1. collect all values of G_high, G_low and s.e.m.s
%       2. plot G_high vs G_low


% 0. initialize meta data and mu stats

clear
clc

cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
load('muData.mat')

% 0. initialize summary vectors for data
dataIndex = find(~cellfun(@isempty,storedMetaData));
experimentCount = length(dataIndex);

% 1. collect all values of G_high, G_low and s.e.m.s
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
    high = 4;
    low = 2;
    
    G_high_means(counter) = muData{index}{1,high}.mean;
    G_high_sems(counter) = muData{index}{1,high}.sem;
    G_low(counter) = muData{index}{1,low}.mean;
    G_low_sems(counter) = muData{index}{1,low}.sem;
    
    timescales_perG(counter) = timescale;
    dates_perG{counter} = date;
    
end

% 2. plot G_high vs G_low
shift = 0.05;

figure(4)
for i = 1:length(timescales_perG)
    if timescales_perG(i) == 30
        plot(G_low(i),G_high_means(i),'o','Color',[1 0 0],'MarkerSize',10); % red
        text(G_low(i)+shift,G_high_means(i)+shift, dates_perG(i),'FontSize',10);
        hold on
    elseif timescales_perG(i) == 300
        plot(G_low(i),G_high_means(i),'o','Color',[1 0.85 0.01],'MarkerSize',10); % sunflower yellow
        text(G_low(i)+shift,G_high_means(i)+shift, dates_perG(i),'FontSize',10);
        hold on
    elseif timescales_perG(i) == 900
        plot(G_low(i),G_high_means(i),'o','Color',[0 0.7 0.7],'MarkerSize',10); % green
        text(G_low(i)+shift,G_high_means(i)+shift, dates_perG(i),'FontSize',10);
        hold on
    elseif timescales_perG(i) == 3600
        plot(G_low(i),G_high_means(i),'o','Color',[0.5 0 0.9],'MarkerSize',10); % purple
        text(G_low(i)+shift,G_high_means(i)+shift, dates_perG(i),'FontSize',10);
        hold on
    end
end
title('no correlation between G_H and G_L')
xlabel('mean mu in low (1/hr)')
ylabel('mean mu in high (1/hr)')
axis([0.7 1.7 2.6 3.4])


%% THREE (A). G_jensen and G_fluc NORMALIZED by G_monod (BVPR)

% Strategy:

%       0.  initialize meta data and bvpr stats
%       0.  initialize summary vectors for data
%       1.  for each experiment, collect values of G_fluc, G_low, G_ave, G_high
%       2.  calculate G_jensens = average growth between high and low
%       3.  normalize G_jensens and G_fluc by G_monod (G_ave)
%       4.  generate hypothetical points, where on the timescale vector:
%                -  0 represents infinitely fast
%                -  5 represents infinitely slow
%       5. plot G_data by timescale, varying brightness bewteen replicates 



% 0. initialize meta data and bvpr stats

clear
clc

cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
load('bioProdRateData.mat')

% 0. initialize summary vectors for data
dataIndex = find(~cellfun(@isempty,storedMetaData));
experimentCount = length(dataIndex);


% 1. for each experiment, collect values of G_fluc, G_low, G_ave, G_high
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
    
    G_fluc(counter) = bioProdRateData{index}{1,fluc}.mean;
    G_low(counter) = bioProdRateData{index}{1,low}.mean;
    G_ave(counter) = bioProdRateData{index}{1,ave}.mean;
    G_high(counter) = bioProdRateData{index}{1,high}.mean;
    
    timescales_perG(counter) = timescale;
    dates_perG{counter} = date;
    
end
clear fluc low ave high counter date timescale index e


% 2. calculate G_jensens = average growth between high and low
G_jensens = (G_high+G_low)/2;


% 3. normalize G_jensens and G_fluc by G_monod (G_ave)
G_jensens_normalized = G_jensens./G_ave;
G_fluc_normalized = G_fluc./G_ave;


% 4. generate hypothetical points, where on the timescale vector:
%       0 represents infinitely fast
%       5 represents infinitely slow

% timescale
timescale_vector(timescales_perG==30) = 1;
timescale_vector(timescales_perG==300) = 2;
timescale_vector(timescales_perG==900) = 3;
timescale_vector(timescales_perG==3600) = 4;

% G_fluctuating
gf_30 = G_fluc_normalized(timescales_perG==30);
gf_300 = G_fluc_normalized(timescales_perG==300);
gf_900 = G_fluc_normalized(timescales_perG==900);
gf_3600 = G_fluc_normalized(timescales_perG==3600);

% hypothetical infinutes
hypothetical_0 = 1;
hypothetical_30 = mean(gf_30);
hypothetical_300 = mean(gf_300);
hypothetical_900 = mean(gf_900);
hypothetical_3600 = mean(gf_3600);
hypothetical_5 = nanmean(G_jensens_normalized);

hypothetical_vector = [hypothetical_0 hypothetical_30 hypothetical_300 hypothetical_900 hypothetical_3600 hypothetical_5];
hypothetical_time = [0 1 2 3 4 5];


% 5. plot G_data by timescale, varying brightness bewteen replicates

figure(5)
color_Monod = rgb('DarkBlue');
color_Jensen = rgb('Teal');
color_fluc = rgb('Chocolate');
color_hyp = rgb('DimGray');

for i = 1:length(G_jensens)
    
    % plot G_monod
    plot(timescale_vector,G_ave./G_ave,'o','Color',color_Monod,'MarkerSize',8,'LineWidth',2); % dark blue
    hold on
    
    % plot G_jensen
    plot(timescale_vector,G_jensens_normalized,'o','Color',color_Jensen,'MarkerSize',8,'LineWidth',2); % teal
    hold on
    
    % plot G_fluc
    plot(timescale_vector,G_fluc_normalized,'o','Color',color_fluc,'MarkerSize',8,'LineWidth',2); % chocolate
    hold on
    
    % plot hypothetical line
    plot(hypothetical_time, hypothetical_vector, '-.', 'Color',color_fluc,'LineWidth',2); % chocolate
    plot(hypothetical_time, [hypothetical_0 NaN NaN NaN NaN hypothetical_5], 'o', 'Color',color_hyp,'LineWidth',2); % dim gray
    
    % plot mean monod and jensens
    plot(hypothetical_time, ones(1,length(hypothetical_time)), ':', 'Color',color_Monod,'LineWidth',2); % dark blue
    plot(hypothetical_time, ones(1,length(hypothetical_time))*hypothetical_5, ':', 'Color',color_Jensen,'LineWidth',2); % teal
    
    
end
axis([0 5 0 1.2])
title('expected and measured mean growth rates')
xlabel('timescale of nutrient fluctuation')
ylabel('mean bvpr as fraction of stable average')


%% THREE (B). G_jensen and G_fluc NORMALIZED by G_monod (MU)

% Strategy:

%       0.  initialize meta data and mu stats
%       0.  initialize summary vectors for data
%       1.  for each experiment, collect values of G_fluc, G_low, G_ave, G_high
%       2.  calculate G_jensens = average growth between high and low
%       3.  normalize G_jensens and G_fluc by G_monod (G_ave)
%       4.  generate hypothetical points, where on the timescale vector:
%                -  0 represents infinitely fast
%                -  5 represents infinitely slow
%       5. plot G_data by timescale, varying brightness bewteen replicates 



% 0. initialize meta data and mu stats

clear
clc

cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
load('muData.mat')

% 0. initialize summary vectors for data
dataIndex = find(~cellfun(@isempty,storedMetaData));
experimentCount = length(dataIndex);


% 1. for each experiment, collect values of G_fluc, G_low, G_ave, G_high
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
    
    G_fluc(counter) = muData{index}{1,fluc}.mean;
    G_low(counter) = muData{index}{1,low}.mean;
    G_ave(counter) = muData{index}{1,ave}.mean;
    G_high(counter) = muData{index}{1,high}.mean;
    
    timescales_perG(counter) = timescale;
    dates_perG{counter} = date;
    
end
clear fluc low ave high counter date timescale index e


% 2. calculate G_jensens = average growth between high and low
G_jensens = (G_high+G_low)/2;


% 3. normalize G_jensens and G_fluc by G_monod (G_ave)
G_jensens_normalized = G_jensens./G_ave;
G_fluc_normalized = G_fluc./G_ave;

% 4. generate hypothetical points, where on the timescale vector:
%       0 represents infinitely fast
%       5 represents infinitely slow

% timescale
timescale_vector(timescales_perG==30) = 1;
timescale_vector(timescales_perG==300) = 2;
timescale_vector(timescales_perG==900) = 3;
timescale_vector(timescales_perG==3600) = 4;

% G_fluctuating
gf_30 = G_fluc_normalized(timescales_perG==30);
gf_300 = G_fluc_normalized(timescales_perG==300);
gf_900 = G_fluc_normalized(timescales_perG==900);
gf_3600 = G_fluc_normalized(timescales_perG==3600);

% hypothetical infinutes
hypothetical_0 = 1;
hypothetical_30 = mean(gf_30);
hypothetical_300 = mean(gf_300);
hypothetical_900 = mean(gf_900);
hypothetical_3600 = mean(gf_3600);
hypothetical_5 = nanmean(G_jensens_normalized);

hypothetical_vector = [hypothetical_0 hypothetical_30 hypothetical_300 hypothetical_900 hypothetical_3600 hypothetical_5];
hypothetical_time = [0 1 2 3 4 5];


% 5. plot G_data by timescale, varying brightness bewteen replicates

figure(6)
color_Monod = rgb('DarkBlue');
color_Jensen = rgb('Teal');
color_fluc = rgb('Chocolate');
color_hyp = rgb('DimGray');

for i = 1:length(G_jensens)
    
    % plot G_monod
    plot(timescale_vector,G_ave./G_ave,'o','Color',color_Monod,'MarkerSize',8,'LineWidth',2); % dark blue
    hold on
    
    % plot G_jensen
    plot(timescale_vector,G_jensens_normalized,'o','Color',color_Jensen,'MarkerSize',8,'LineWidth',2); % teal
    hold on
    
    % plot G_fluc
    plot(timescale_vector,G_fluc_normalized,'o','Color',color_fluc,'MarkerSize',8,'LineWidth',2); % chocolate
    hold on
    
    % plot hypothetical line
    plot(hypothetical_time, hypothetical_vector, '-.', 'Color',color_fluc,'LineWidth',2); % chocolate
    plot(hypothetical_time, [hypothetical_0 NaN NaN NaN NaN hypothetical_5], 'o', 'Color',color_hyp,'LineWidth',2); % dim gray
    
    % plot mean monod and jensens
    plot(hypothetical_time, ones(1,length(hypothetical_time)), ':', 'Color',color_Monod,'LineWidth',2); % dark blue
    plot(hypothetical_time, ones(1,length(hypothetical_time))*hypothetical_5, ':', 'Color',color_Jensen,'LineWidth',2); % teal
    
    
end
axis([0 5 0 1.2])
title('expected and measured mean growth rates')
xlabel('timescale of nutrient fluctuation')
ylabel('mean mu as fraction of stable average')
