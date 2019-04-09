%% figure 58


%  Goals: test whether conditions run on the same day correlate
%         plots: (1) steady-state Low vs steady-state High
%                (2) steady-state Low vs steady-state Ave
%                (3) steady-state Ave vs steady-state High
%                (4) steady-state Ave vs fluc

%  Strategy:
%


%  last updated: jen, 2019 April 8

%  commit: first commit, correlation between conditions run on same day


% OK let's go!

%% initialize

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Writing/manuscript 1/figure3/')
load('growthRates_monod_curve.mat')


% 0. initialize experiments to plot
exptArray = [2:7,9:15];


%% compile data

% 1. initialize vectors for concatenating same condition data
fluc = [];
low = [];
ave = [];
high = [];

for e = 1:length(exptArray)
    
    % 2. access growth rate stats for each experiment
    index = exptArray(e);
    exptData = growthRates_monod_curve{1,index};
    
    % 3. loop through conditions and collect mean
    for c = 1:4
        
        condData = exptData{1,c};
        if isempty(condData) == 1
            condMean = NaN;
        else
            condMean = condData.mean;
        end
        
        
        if c == 1
            fluc = [fluc; condMean];
        elseif c == 2
            low = [low; condMean];
        elseif c == 3
            ave = [ave; condMean];
        elseif c == 4
            high = [high; condMean];
        end
        
    end


end

%% plot data - steady-state Low vs steady-state High
    

% 1. prep data
rr = find(~isnan(high)); % rows to plot, avoiding NaN
h = high(rr);
l = low(rr);


% 2. linear regression and correlation coefficient

% fit linear regression
p = polyfit(l,h,1); 
x = l;
y = p(1)*x + p(2);

% calculate correlation coefficient
r = corrcoef(l,h); 
R = r(1,2);


% 3. plot
figure(1)
plot(l,h,'o','MarkerSize',6)
hold on
plot(x,y,'Color',rgb('SlateGray'),'LineWidth',2)
hold on
txt = strcat('r=',num2str(R));
text(x(end),y(end),txt,'FontSize',14)
xlabel('growth rate in low')
ylabel('growth rate in high')
title('correlation of Glow and Ghigh')

%% plot data - steady-state Low vs steady-state Ave

% 1. prep data
rr = find(~isnan(low)); % rows to plot, avoiding NaN
a = ave(rr);
l = low(rr);


% 2. linear regression and correlation coefficient

% fit linear regression
p = polyfit(l,a,1); 
x = l;
y = p(1)*x + p(2);

% calculate correlation coefficient
r = corrcoef(l,a); 
R = r(1,2);


% 3. plot
figure(2)
plot(l,a,'o','MarkerSize',6)
hold on
plot(x,y,'Color',rgb('SlateGray'),'LineWidth',2)
hold on
txt = strcat('r=',num2str(R));
text(x(end),y(end),txt,'FontSize',14)
xlabel('growth rate in low')
ylabel('growth rate in ave')
title('correlation of Glow and Gave')

%% plot data - steady-state Ave vs steady-state High

% 1. prep data
rr = find(~isnan(high)); % rows to plot, avoiding NaN
h = high(rr);
a = ave(rr);


% 2. linear regression and correlation coefficient

% fit linear regression
p = polyfit(a,h,1); 
x = a;
y = p(1)*x + p(2);

% calculate correlation coefficient
r = corrcoef(a,h); 
R = r(1,2);


% 3. plot
figure(3)
plot(a,h,'o','MarkerSize',6)
hold on
plot(x,y,'Color',rgb('SlateGray'),'LineWidth',2)
hold on
txt = strcat('r=',num2str(R));
text(x(end),y(end),txt,'FontSize',14)
xlabel('growth rate in ave')
ylabel('growth rate in high')
title('correlation of Gave and Ghigh')

%% plot data - steady-state Ave vs G_fluc

% 1. prep data
rr = find(~isnan(ave)); % rows to plot, avoiding NaN
a = ave(rr);
f = fluc(rr);


% 2. linear regression and correlation coefficient

% fit linear regression
p = polyfit(a,f,1); 
x = a;
y = p(1)*x + p(2);

% calculate correlation coefficient
r = corrcoef(a,f); 
R = r(1,2);


% 3. plot
figure(4)
plot(a,f,'o','MarkerSize',6)
hold on
plot(x,y,'Color',rgb('SlateGray'),'LineWidth',2)
hold on
txt = strcat('r=',num2str(R));
text(x(end),y(end),txt,'FontSize',14)
xlabel('growth rate in ave')
ylabel('growth rate in fluc')
title('correlation of Gave and Gfluc')

%% save plots

% 1. save plots
cd('/Users/jen/Documents/StockerLab/Data_analysis/currentPlots/')

figure(1)
plotName = strcat('figure58-fig1-Glow_v_Ghigh');
saveas(gcf,plotName,'epsc')

figure(2)
plotName = strcat('figure58-fig2-Glow_v_Gave');
saveas(gcf,plotName,'epsc')

figure(3)
plotName = strcat('figure58-fig3-Gave_v_Ghigh');
saveas(gcf,plotName,'epsc')

figure(4)
plotName = strcat('figure58-fig4-Gave_v_Gfluc');
saveas(gcf,plotName,'epsc')
    



