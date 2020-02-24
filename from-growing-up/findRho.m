%% findRho

% Goal: solve for rho, the population growth rate
%
% 
% As written in Celani's Notes ('Can turbulence increase bacterial growth rates?'),  
% a population with a distribution of division times with always grow more
% quickly than a population in which individials divide at a uniform rate,
% given the same average division time.
% 
% This benefit of heterogeneity is due to Jensen's inequality. For those
% that grow slower than average, those that grow faster will more than make
% up for the loss - thanks to the exponential nature of growth.
%
% In this code, I attempt to solve for rho using division time data pulled
% from 6 experiments of E coli growing in steady-state growth in steady a
% stable glucose environment (0.5 uM).
%
% Please see distributeDoublingTime.m for the measured and fitted
% distributions of this data.
%
% 
% 
% Strategy:
%
%  I have a bunch of T values and the corresponding probabilities, f(T) values.
%  In other words, a division time T=30 min has some probability f(30min).
%  
%  These are my measured inputs for Celani's equality:
%
%           integral(zero and infinity) [ exp(-rho * T) * f(T)dt] = 1/2
%
%
%  Here, I will instead treat the integral as sum, using my data as discrete points.
%  So now, the integral has becom:
%
%           sum[ exp(-rho * T) * f(T) ]
%
%
%  Seeing as we don't yet know the value of rho, let's call this sum, S(rho).
%  Then...
%
%
%           0. initialize division time data
%           1. calculate PDF, or f(T) array over from T=0 to T=300
%           2. define function S(rho)
%           3. plot S(rho) over a range of rho values
%           4. find value of rho such that S(rho) = 0.5
% 
%
%
%
% last edit: jen, 2017 jun 16
%
% OK lez go!!

%%
clear

% 0. initialize division time data
newFolder = strcat('/Users/jen/Documents/StockerLab/Data analysis/PDF_divisionTime');
cd(newFolder);
clear newFolder;

% 0. load division time variables from distributeDoublingTime.m
load('dtVars');


% 1. calculate PDF, or f(T) array over from T=0 to T=300
T = linspace(0,300,301)';
f_T = normpdf(T, normalFit.mu, normalFit.sigma);


% 2. define and use function S(rho)
rho = linspace(0,0.01,101)';
S_rho = sumRho(rho,T,f_T);

% 3. plot S(rho) over a range of rho values
plot(rho,S_rho)
xlabel('Rho')
ylabel('S(Rho)')
    

% 4. find rho for which S(rho) = 0.5
y = ones(101,1)*0.5;

hold on
plot(rho,y,'r.')

%%

% 4. find value of rho* such that S(rho*) = 1/2
%from plot:
rho_star = 0.0056;


% 5. compute term for growth of a population with uniform growth rate

% compute expected value from distribution
byPairs = f_T.*T; % probability * value
e_T = sum(byPairs);

uni_Growth = exp(-rho_star * e_T);

%%
