%% figure 3C - expected and measured growth rates


% Output: Calculated percent change from reference growth rate, and
%         calculated errors (accounting for error combination)
%           

% Input:  Mean growth rate (i.e. G_fluc, G_low, G_ave, G_high) of each replicate


% Percent change calculations:
%           1. from G_ave
%           2. from G_Jensens
%
%    example: relative change = (G_fluc - G_ave)/G_ave * 100


% Error calculations:
%
%         Two parts: "top" accounts for subtraction step of % change calculation
%                    "final" additionally accounts for division step
%
%         error_top = sqrt( err_fluc^2 + err_ave^2 )
%
%         error_final = ( (G_fluc - G_ave)/G_ave ) * sqrt( (error_top/(G_fluc - G_ave))^2 + (err_ave/G_ave)^2 )
%
%         Both error_top and error_final are functions performing these
%         calculations. For details, see "How to combine errors" by Robin Hogan, 2006




% Last edit: jen, 2019 October 19
% Commit: G_fluc relative to G_ave and G_jensens with standard error


% OK let's go!

%% Part One. Initialize data

clear
clc


% 0. initialize column IDs of each nutrient condition
fluc = 1; low = 2; ave = 3; high = 4;


% 0. initialize row IDs of each timescale in mean growth rate data
t30s = 1:3; t5 = 4:6; t15 = 7:10; t60 = 11:13;
timescales = {t30s; t5; t15; t60};


% 0. initialize mean growth rate data, G
daily_G = [
    1.8391, 1.0965, 2.4439, 2.9497; % 30 sec, 2017-11-12
    2.1229, 1.2170, 2.4692, 2.6744; % 30 sec, 2017-11-14
    1.8373, 0.8430, 2.2755, 2.9907; % 30 sec, 2018-01-04
    
    1.6728, 1.4814, 2.6296, 2.8919; % 5 min, 2017-10-10
    1.6198, 1.3525, 2.3198, 2.9588; % 5 min, 2017-11-15
    1.3032, 0.9503, 2.1744,    NaN; % 5 min, 2018-01-11
    
    1.5064, 1.3282, 2.4434, 3.1206; % 15 min, 2017-11-13
    1.3302,    NaN, 2.4605,    NaN; % 15 min, 2018-01-12
    0.9695, 0.9870, 2.0733, 2.7382; % 15 min, 2018-01-16
    0.8595, 1.0341, 2.0362, 2.7467; % 15 min, 2018-01-17 

    1.2941, 0.7944, 2.3412, 2.8273; % 60 min, 2018-01-29 
    1.1211, 0.8357, 2.1422, 2.7684; % 60 min, 2018-01-31
    1.0272, 0.9130, 2.2123, 2.7572; % 60 min, 2018-02-01

    ];




% 1. calculate daily Jensens
daily_Jensens = (daily_G(:,low) + daily_G(:,high))/2;  



% 2. relative change calculations from G_ave and from G_jensens
daily_change_ave = (daily_G(:,fluc) - daily_G(:,ave))./daily_G(:,ave) * 100;   
daily_change_jensens = (daily_G(:,fluc) - daily_Jensens)./daily_Jensens * 100;   


%% Part Two. Perform calculations and plot

% 3. collect mean and standard error within a timescale
daily_means = zeros(4,2); % row is timescale: 30s, 5 min, 15 min, 60 min
daily_error = zeros(4,2); % column 1 = from ave; 
                          % column 2 = from Jensens;

for tt = 1:length(timescales)
    
    daily_means(tt,1) = nanmean(daily_change_ave(timescales{tt}));
    daily_means(tt,2) = nanmean(daily_change_jensens(timescales{tt}));
    
    sdev = nanstd(daily_change_ave(timescales{tt}));
    count = length(daily_change_ave(timescales{tt}));
    daily_error(tt,1) = sdev./sqrt(count);
    
    sdev_j = nanstd(daily_change_jensens(timescales{tt}));
    count_j = sum(~isnan(daily_change_jensens(timescales{tt})));
    daily_error(tt,2) = sdev_j./sqrt(count_j);
    
end



% 4. bar plot of percent change
spacing = [0.73, 0.91, 1.09, 1.27;
           1.73, 1.91, 2.09, 2.27];
       
    
       
figure(1)
bar(daily_means')
hold on
errorbar(spacing,daily_means',daily_error','.','Color',rgb('Black'))
ylabel('percent change')
xlabel('reference growth rate')
title('percent change in growth rate')





