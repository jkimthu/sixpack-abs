% plotMonod: biovolume production rate vs conc

% goal: plot monod curve using compiled LB data

% strategy:
%
%       0. initialize experiment data
%       1. for each experiment...load data
%               -  identify date
%               -  continue to exclude outlier: 2017-10-31
%               2. for each condition in experiment...
%                       3. isolate condition specific Mu_va and time data
%                       4. remove mu data with timestamps prior to and after stabilization
%                       5. remove zeros (always two at start and end of track) and negatives
%                       6. calculate: biovolume production rate = V(t) * mu(t) * ln(2)
%                       7. isolate data to stabilized regions of growth
%                       8. calculate average and s.e.m. of stablized data
%                       9. accumulate data for storage / plotting
%               10. store data from all conditions into measured data structure
%      11.  save stored data into data structure
%      12.  plot average biovolume production rate over time
%


% last updated: 2018 Mar 22
% commit: update comments to better reflect strategy


%% Add NEW experiment to bioProdRate data structure 

% 0. initialize experiment data
clear
clc

cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')

dataIndex = find(~cellfun(@isempty,storedMetaData));
bioProdRateData = cell(size(storedMetaData));
muData = cell(size(storedMetaData));

% initialize summary vectors for calculated data
experimentCount = length(dataIndex);

% 1. for each experiment, move to folder and load data

for e = 1:experimentCount
    
    % identify experiment by date
    index = dataIndex(e);
    date = storedMetaData{index}.date;
    
    % move directory to experiment data
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    
    % load data
    timescale = storedMetaData{index}.timescale;
    if ischar(timescale) == 0
        filename = strcat('lb-fluc-',date,'-window5-width1p4-1p7-jiggle-0p5.mat');
    elseif strcmp(date,'2017-09-26') == 1
        filename = 'lb-monod-2017-09-26-window5-va-jiggle-c12-0p1-c3456-0p5-bigger1p8.mat';
    elseif strcmp(date, '2017-11-09') == 1
        filename = 'lb-control-2017-11-09-window5-width1p4-jiggle-0p5.mat';
    end
    load(filename,'D5','M','M_va','T')
    
    % build experiment data matrix
    display(strcat('Experiment (', num2str(e),') of (', num2str(length(dataIndex)),')'))
    xy_start = 1;
    xy_end = length(D5);
    exptData = buildDM(D5,M,M_va,T,xy_start,xy_end,e);
    
    clear D5 M M_va T filename experimentFolder
   
    % 2. for each condition, calculate mean biovolume production rate per condition
    xys = storedMetaData{index}.xys;
    xy_dimensions = size(xys);
    totalConditions = xy_dimensions(1);
    
    for c = 1:totalConditions
        
        % 3. isolate all data from current condition
        conditionData = exptData(exptData(:,28) == c,:);
        
        % 4. isolate volume (Va), mu (mu_va) and time data from current condition
        volumes = conditionData(:,14);        % col 14 = calculated va_vals (cubic um)
        mus = conditionData(:,17);            % col 17 = calculated mu_va 
        timestamps = conditionData(:,2)/3600; % time in seconds converted to hours
        clear conditionData
        
        % 5. remove data for which mu = 0, as these were the edges of tracks that never get calculated
        trueMus = mus(mus > 0);
        trueVols = volumes(mus > 0);
        trueTimes = timestamps(mus > 0);
        clear volumes mus timestamps
        
        % 6. calculate: biovolume production rate = V(t) * mu(t) * ln(2)
        bioProdRate = trueVols .* trueMus * log(2); % log(2) in matlab = ln(2)
        clear trueVols  
        
        % 7. isolate data to stabilized regions of growth
        minTime = 3;  % hr
        maxTime = storedMetaData{index}.bubbletime(c);
        
        times_trim1 = trueTimes(trueTimes >= minTime);
        bioProdRate_trim1 = bioProdRate(trueTimes >= minTime);
        mu_trim1 = trueMus(trueTimes >= minTime);
        clear trueTimes trueMus
        
        if maxTime > 0
            bioProdRate_trim2 = bioProdRate_trim1(times_trim1 <= maxTime);
            mu_trim2 = mu_trim1(times_trim1 <= maxTime);
        else
            bioProdRate_trim2 = bioProdRate_trim1;
            mu_trim2 = mu_trim1;
        end
        clear times_trim1
        
        
        % 8. calculate average and s.e.m. of stabilized data
        mean_bioProdRate = mean(bioProdRate_trim2);
        count_BioProdRate = length(bioProdRate_trim2);
        std_BioProdRate = std(bioProdRate_trim2);
        sem_BioProdRate = std_BioProdRate./sqrt(count_BioProdRate);
        
        mean_mu = mean(mu_trim2);
        count_mu = length(mu_trim2);
        std_mu = std(mu_trim2);
        sem_mu = std_mu./sqrt(count_mu);
        
        
        % 9. accumulate data for storage / plotting
        compiledbioProdRate{c}.mean = mean_bioProdRate;
        compiledbioProdRate{c}.std = std_BioProdRate;
        compiledbioProdRate{c}.count = count_BioProdRate;
        compiledbioProdRate{c}.sem = sem_BioProdRate;
        
        compiledMu{c}.mean = mean_mu;
        compiledMu{c}.std = std_mu;
        compiledMu{c}.count = count_mu;
        compiledMu{c}.sem = sem_mu;
        

        clear mean_bioProdRate std_BioProdRate count_BioProdRate sem_BioProdRate
        clear bioProdRate bioProdRate_trim1 bioProdRate_trim2 maxTime 
        clear mean_mu std_mu count_mu sem_mu mu_trim1 mu_trim2
    
    end
    
    % 10. store data from all conditions into measured data structure        
    bioProdRateData{index} = compiledbioProdRate;
    muData{index} = compiledMu;
    
    clear compiledbioProdRate compiledMu
end


%% 11. Save new data into stored data structure
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
save('bioProdRateData.mat','bioProdRateData')
save('muData.mat','muData')

%% 12. plot average biovolume production rate over time
clc
clear

cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
load('bioProdRateData.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));
experimentCount = length(dataIndex);
latestExpt = '2018-01-31';

% initialize summary stats for fitting
counter = 0;
summaryMeans = zeros(1,(experimentCount-1)*3 + 6);
summaryConcentrations = zeros(1,(experimentCount-1)*3 + 6);


for e = 1:experimentCount
    
    % identify experiment by date
    index = dataIndex(e);
    date = storedMetaData{index}.date;
    
    % exclude outlier from analysis
    if strcmp(date, '2017-10-31') == 1
        disp(strcat(date,': excluded from analysis'))
        continue
    end
    disp(strcat(date, ': analyze!'))
    
%     % collect analyzed dates for legend
%     analyzedDates{e} = date;
    
    % load timescale
    timescale = storedMetaData{index}.timescale;
    
    % isolate biomass prod data for current experiment
    experimentData = bioProdRateData{index};
    
    % isolate concentration data for current experiment
    concentration = storedMetaData{index}.concentrations;
    
    % plot, labeled by experiment date
    figure(1)
    for c = 1:length(concentration)
        h(e) = errorbar(log(concentration(c)), experimentData{c}.mean, experimentData{c}.sem,'o','Color',[1 1 1]*e*.05,'MarkerSize',10);
        hold on
    end 
    legend(e)
    ylabel('biomass prodution rate (cubic um/hr)')
    xlabel('log fold LB dilution')
    
    % plot individual data, labeled by stable vs fluc
    figure(2)
    for c = 1:length(concentration)
        % if fluc experiment
        if ischar(timescale)
            errorbar(log(concentration(c)), experimentData{c}.mean, experimentData{c}.sem,'o','Color','k','MarkerSize',10);
            hold on

            % for stable conditions, accumulate data into summary vector
            counter = counter + 1;
            summaryMeans(counter) = experimentData{c}.mean;
            summaryConcentrations(counter) = concentration(c);

        elseif timescale == 30 && c == 1
            errorbar(log(concentration(c)), experimentData{c}.mean, experimentData{c}.sem,'o','Color',[0.25 0.25 0.9],'MarkerSize',10);
            hold on
        elseif timescale == 300 && c == 1
            errorbar(log(concentration(c)), experimentData{c}.mean, experimentData{c}.sem,'o','Color',[0 .7 .7],'MarkerSize',10);
            hold on
            legend('5 min')
        elseif timescale == 900 && c == 1
            errorbar(log(concentration(c)), experimentData{c}.mean, experimentData{c}.sem,'o','Color',[1 0.6 0],'MarkerSize',10);
            hold on
        elseif timescale == 3600 && c == 1
            errorbar(log(concentration(c)), experimentData{c}.mean, experimentData{c}.sem,'o','Color',[1 0.5 0.5],'MarkerSize',10);
            hold on
        else
            errorbar(log(concentration(c)), experimentData{c}.mean, experimentData{c}.sem,'o','Color','k','MarkerSize',10);
            hold on
            
            % for stable conditions, accumulate data into summary vector
            counter = counter + 1;
            summaryMeans(counter) = experimentData{c}.mean;
            summaryConcentrations(counter) = concentration(c);
        end
    end
    legend('30 sec','5 min','15 min','60 min','stable')
    ylabel('biomass prod rate (cubic um/hr)')
    xlabel('log fold LB dilution')
    title(strcat('up to :',latestExpt))
    
end

%% 13. calculate and plot monod fit for data


% initialize concentration and biovol production rate data for stable
% environments only
conc = summaryConcentrations;
bioVolProdRate = summaryMeans;

% lowest growth rate as y-intercept
C = min(summaryMeans);

% highest growth rate as Vmax
Vmax = max(summaryMeans);

% calculate fit using nonlinear regression
michaelisMenten = @(b,x)( (Vmax*x) ./ (b+x) + C);
beta0 = 0.01;
beta = nlinfit(conc,bioVolProdRate,michaelisMenten,beta0);

% generate fit data
Km = beta;
substrate = 0:0.0001:1;
for s = 1:10001
    V(s) = Vmax* ( substrate(s) / (Km + substrate(s)) ) + C;
end


figure(5)
plot(substrate,V,'r')
hold on
plot(summaryConcentrations, summaryMeans, 'o')
legend('fit','data')
xlabel('strength LB (fraction of full)')
ylabel('biovol production rate (um3/hr)')

figure(6)
plot(log(substrate),V,'r')
hold on
plot(log(summaryConcentrations),summaryMeans,'o')
legend('fit','data')
xlabel('log LB dilution')
ylabel('biovol production rate (um3/hr)')


%% Calculating mean growth at each stable concentration across experiments

highMeans = summaryMeans(summaryConcentrations == 0.02);
aveMeans = summaryMeans(summaryConcentrations == 0.0105);
lowMeans = summaryMeans(summaryConcentrations == 0.001);

figure(7)
errorbar(1, mean(lowMeans), std(lowMeans),'o')
hold on
errorbar(2, mean(aveMeans), std(aveMeans),'o')
hold on
errorbar(3, mean(highMeans),std(highMeans),'o')
xlabel('condition')
ylabel('mean biovolume prod rate with st dev')
legend('stable low','stable ave','stable high')
axis([0 4 0 16])
