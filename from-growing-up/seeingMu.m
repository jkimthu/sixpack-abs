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
%%
%
%  VISUAL CHECK: plot raw data and mu over time
%

% Load workspace from SlidingFits.m     (should be Year-Mon-Day-Mus-length.m)
load('2016-03-16-Mus-length.mat');

counter =0;
for n = 1:10:60
    counter = counter +1;
    m = 5;
    
    % Extracted mu
    Mu_track = M6{n}(m).Parameters(:,1);
    vectorLength = length(Mu_track);
    
    % Original length data (microns)
    Ltrack2 = D6{n}(m).MajAx(3:vectorLength+2);                                  % trimmed length trajectory to match Mu
    
    % Time data (hours)
    %dT = mean(mean(diff(T)));                                              % mean time between frames (seconds)
    %Ttrack = D6{n}(m).Frame(3:Num_mu+2);                                   % original frame # in trajectory
    timeTrack = T{n}/(60*60);
    
    figure(1)
    
    subplot(6,1,counter)
    plot(timeTrack(3:vectorLength+2),Ltrack2,'.',timeTrack(3:vectorLength+2),Mu_track*log(2),'r.');                          
    grid on;
    axis([0,11.3,-0.5,6])
    xlabel('Time (hours)')
    ylabel('Cell Length (um)')
    legend('Length','Mu');


    clear Mu_track Num_mu Ltrack2 Ttrack hr;

end
 

%% CHECK FOUR: plot average growth rate over time
%
%       - generates a single plot with all conditions 
%       - options to plot standard deviation or standard error
%       - saves average mu, standard dev, s.e.m., and number of tracks per bin per condition 
%       

% Initialize
clear;
load('2015-03-16-Mus-length.mat','D6','M6','T');
Mu_stats = {};



%    Condition One    %
Mu_cond1 = [];
Time_cond1 = [];

for n=31:40
    for m = 1:length(M6{n})
        
        %  assemble all instantaneous growth rates into a single vector
        Mu_cond1 = [Mu_cond1; M6{n}(m).Parameters(:,1)];
        
        %  assemble a corresponding timestamp vector
        vectorLength = length(M6{n}(m).Parameters(:,1));
        trackFrames = D6{n}(m).Frame(3:vectorLength+2);
        Time_cond1 = [Time_cond1; T{n}(trackFrames)];
        
    end
end

%  convert all timestamps from seconds to hours
Time_cond1 = Time_cond1/3600;

%  eliminate negative growth rates
%Mu_cond1(Mu_cond1<0)=NaN;

%  determine size of time bins 
Bins = ceil(Time_cond1*20);            % multiplying by 200 gives time bins of 0.005 hr

%  accumulate growth rates by bin, and calculate mean and std dev
Mu_Means = accumarray(Bins,Mu_cond1,[],@nanmean);
Mu_STDs = accumarray(Bins,Mu_cond1,[],@nanstd);


%   to calculate s.e.m.
%   1. count number of total tracks in each bin        
for j = 1:max(Bins)
    currentBin_count = find(Bins==j);        
    counter = 1;                   
    
    for i = 2:length(currentBin_count)
        if currentBin_count(i) == currentBin_count(i-1)+1;
            counter = counter;
        else
            counter = counter + 1;
        end
    end
    Mu_Counts(j) = counter;        
    clear i counter Kasten;
end

%   2. divide standard dev by square root of tracks per bin
Mu_sems = Mu_STDs./sqrt(Mu_Counts');

plot(Mu_Means,'k')
hold on
errorbar(Mu_Means,Mu_sems,'k')
hold on
%axis([0,11.3,0,.7])
axis([0,230,0,.4])
xlabel('Time (hours)')
ylabel('Elongation rate (1/hr)') 

% Saving stats
% Mu_stats(:,1) = {Mu_Means};
% Mu_stats(:,2) = {Mu_STDs};
% Mu_stats(:,3) = {Mu_sems};
% Mu_stats(:,4) = {Mu_Counts'};

clear vectorLength trackFrams Mu_Means Mu_STDs Mu_sems Bins hr dT Mu_Counts n m j;
clear Mu_cond1 Time_cond1;



%     Condition Two     %

Mu_cond2 = [];
Time_cond2 = [];

for n=41:50
    for m = 1:length(M6{n})
        
        %  assemble all instantaneous growth rates into a single vector
        Mu_cond2 = [Mu_cond2; M6{n}(m).Parameters(:,1)];
        
        %  assemble a corresponding timestamp vector
        vectorLength = length(M6{n}(m).Parameters(:,1));
        trackFrames = D6{n}(m).Frame(3:vectorLength+2);
        Time_cond2 = [Time_cond2; T{n}(trackFrames)];
        
    end
end

%  convert all timestamps from seconds to hours
Time_cond2 = Time_cond2/3600;

%  eliminate negative growth rates
%Mu_cond2(Mu_cond2<0)=NaN;

%  determine size of time bins 
Bins = ceil(Time_cond2*20);            % multiplying by 200 gives time bins of 0.005 hr

Mu_Means = accumarray(Bins,Mu_cond2,[],@nanmean);
Mu_STDs = accumarray(Bins,Mu_cond2,[],@nanstd);
       
for j = 1:max(Bins)
    currentBin_count = find(Bins==j);         
    counter = 1;                  
    
    for i = 2:length(currentBin_count)
        if currentBin_count(i) == currentBin_count(i-1)+1;
            counter = counter;
        else
            counter = counter + 1;
        end
    end
    
    Mu_Counts(j) = counter;                          
    clear i counter Kasten;
end

Mu_sems = Mu_STDs./sqrt(Mu_Counts');
                                    
plot(Mu_Means,'b')
hold on
errorbar(Mu_Means,Mu_sems,'b')
grid on
legend('condition 4', 'condition 5');

% Saving stats
% Mu_stats(:,5) = {Mu_Means};
% Mu_stats(:,6) = {Mu_STDs};
% Mu_stats(:,7) = {Mu_sems};
% Mu_stats(:,8) = {Mu_Counts'};

clear Num_mu Ttrack Mu_Means Mu_STDs Mu_sems Bins hr dT n m j n_track;
clear Mu_cond2 Time_cond2;


%save('2015-08-18-Mu-length.mat','D6','M6','Mu_stats','T');

%%
%
% Plot <growth rate> and std in bin 6 from all experiments
%
%

%  Data origin  %%
%
%                                  Monod   Condition                                   
%
%       1.   0 mM         Feb 13   ( 1 )       1
%       2.   0 mM         Feb 15   ( 2 )       1
%       3.   0 mM         Apr 20   ( 3 )       1
%       4.   0 mM         Jun 03   ( 6 )       1
%       
%       4.   16 nM        Apr 20   ( 3 )       2
%       5.   160 nM       Apr 20   ( 3 )       3
%       6.   160 nM       Feb 15   ( 2 )       2
%       7.   1.6 uM       Feb 15   ( 2 )       3
%       8.   16 uM        Apr 20   ( 3 )       4
%       9.   160 uM       Feb 15   ( 2 )       4
%       10.  1.6 mM       Feb 13   ( 1 )       2
%       11.  8 mM         Feb 13   ( 1 )       3
%       12.  16 mM        Feb 13   ( 1 )       4



% Initialize

load('Monod_compiled.mat');
Bin = 6;
Total_conditions = 20;
Monod = zeros(Total_conditions,4);

% in column 1 (conc.) entered '1' instead of '0' for log scale

% column 1 = glucose concentration (nM)
% column 2 = average Mu in Bin
% column 3 = standard deviation of Mus in Bin
% column 4 = standard error of the mean
%Monod(1,:) =  [1,        Stats_0213(Bin,1),  Stats_0213(Bin,2),  Stats_0213(Bin,2)/sqrt(Binpop_compiled{1}(Bin,1))];      
%Monod(2,:) =  [1,        Stats_0215(Bin,1),  Stats_0215(Bin,2),  Stats_0215(Bin,2)/sqrt(Binpop_compiled{2}(Bin,1))];      
%Monod(3,:) =  [1,        Stats_0420(Bin,1),  Stats_0420(Bin,2),  Stats_0420(Bin,2)/sqrt(Binpop_compiled{3}(Bin,1))];
Monod(4,:) =  [16,       Stats_0420(Bin,3),  Stats_0420(Bin,4),  Stats_0420(Bin,4)/sqrt(Binpop_compiled{3}(Bin,2))];
Monod(5,:) =  [160,      Stats_0420(Bin,5),  Stats_0420(Bin,6),  Stats_0420(Bin,6)/sqrt(Binpop_compiled{3}(Bin,3))];
Monod(6,:) =  [160,      Stats_0215(Bin,3),  Stats_0215(Bin,4),  Stats_0215(Bin,4)/sqrt(Binpop_compiled{2}(Bin,2))];
Monod(7,:) =  [1600,     Stats_0215(Bin,5),  Stats_0215(Bin,6),  Stats_0215(Bin,6)/sqrt(Binpop_compiled{2}(Bin,3))];
Monod(8,:) =  [16000,    Stats_0420(Bin,7),  Stats_0420(Bin,8),  Stats_0420(Bin,8)/sqrt(Binpop_compiled{3}(Bin,4))];
Monod(9,:) =  [160000,   Stats_0215(Bin,7),  Stats_0215(Bin,8),  Stats_0215(Bin,8)/sqrt(Binpop_compiled{2}(Bin,4))];
Monod(10,:) = [1600000,  Stats_0213(Bin,3),  Stats_0213(Bin,4),  Stats_0213(Bin,4)/sqrt(Binpop_compiled{1}(Bin,2))];
Monod(11,:) = [8000000,  Stats_0213(Bin,5),  Stats_0213(Bin,6),  Stats_0213(Bin,6)/sqrt(Binpop_compiled{1}(Bin,3))];
Monod(12,:) = [16000000, Stats_0213(Bin,7),  Stats_0213(Bin,8),  Stats_0213(Bin,8)/sqrt(Binpop_compiled{1}(Bin,4))];
%%
Monod(1,:) =  [1,        Stats_0213(12,1),   Stats_0213(12,2),   Stats_0213(12,2)/sqrt(Binpop_compiled{1}(12,1))];      
Monod(2,:) =  [1,        Stats_0215(12,1),   Stats_0215(12,2),   Stats_0215(12,2)/sqrt(Binpop_compiled{2}(12,1))];      
Monod(3,:) =  [1,        Stats_0420(18,1),   Stats_0420(18,2),   Stats_0420(18,2)/sqrt(Binpop_compiled{3}(18,1))]; 
Monod(4,:) =  [1,        Stats_0603{1}(10),  Stats_0603{2}(10),  Stats_0603{3}(10)];
Monod(5,:) =  [16,       Stats_0420(18,3),   Stats_0420(18,4),   Stats_0420(18,4)/sqrt(Binpop_compiled{3}(18,2))];
Monod(6,:) =  [16,       Stats_0506{5}(Bin), Stats_0506{6}(Bin), Stats_0506{7}(Bin)];
Monod(7,:) =  [16,       Stats_0603{5}(10),  Stats_0603{6}(10),  Stats_0603{7}(10)];

Monod(8,:) =  [80,       Stats_0506{9}(Bin), Stats_0506{10}(Bin), Stats_0506{11}(Bin)];
Monod(9,:) =  [160,      Stats_0506{13}(Bin), Stats_0506{14}(Bin), Stats_0506{15}(Bin)];

Monod(10,:) =  [160,      Stats_0420(18,5),   Stats_0420(18,6),   Stats_0420(18,6)/sqrt(Binpop_compiled{3}(18,3))];
Monod(11,:) =  [160,      Stats_0215(Bin,3),  Stats_0215(Bin,4),  Stats_0215(Bin,4)/sqrt(Binpop_compiled{2}(Bin,2))];
Monod(12,:) =  [1600,     Stats_0215(Bin,5),  Stats_0215(Bin,6),  Stats_0215(Bin,6)/sqrt(Binpop_compiled{2}(Bin,3))];
Monod(13,:) =  [16000,    Stats_0420(18,7),   Stats_0420(18,8),   Stats_0420(18,8)/sqrt(Binpop_compiled{3}(18,4))];
Monod(14,:) =  [16000,    Stats_0603{9}(10),  Stats_0603{10}(10), Stats_0603{11}(10)];
Monod(15,:) =  [16000,    Stats_0603{13}(10), Stats_0603{14}(10), Stats_0603{15}(10)];
Monod(16,:) =  [160000,   Stats_0215(Bin,7),  Stats_0215(Bin,8),  Stats_0215(Bin,8)/sqrt(Binpop_compiled{2}(Bin,4))];
Monod(17,:) = [1600000,  Stats_0213(Bin,3),  Stats_0213(Bin,4),  Stats_0213(Bin,4)/sqrt(Binpop_compiled{1}(Bin,2))];
Monod(18,:) = [8000000,  Stats_0213(Bin,5),  Stats_0213(Bin,6),  Stats_0213(Bin,6)/sqrt(Binpop_compiled{1}(Bin,3))];
Monod(19,:) = [16000000, Stats_0213(Bin,7),  Stats_0213(Bin,8),  Stats_0213(Bin,8)/sqrt(Binpop_compiled{1}(Bin,4))];
Monod(20,:) = [16000000, Stats_0506{1}(Bin), Stats_0506{2}(Bin), Stats_0506{3}(Bin)];

figure ()
errorbar(log10(Monod(:,1)),Monod(:,2),Monod(:,4),'.r')
grid on
hold on
plot(log10(Monod(:,1)),Monod(:,2),'.')
hold on
plot(log10(C),U)
axis([-.1,8,0,.6])
xlabel('Log10 glucose concentration (nM)')
ylabel('Elongation rate (1/hr)')
%
%%
%   Concept testing: Michaelis-Menten   %
%   
%   - Construct perfect function and observe effect of scale manipulations
%   - Test fit against selected points from this graph 
%   

% Formulation
%           
%       U  =  (Umax * C) / (Km + C)

%       where,  U   =  growth rate
%               C   =  glucose concentration
%               Km  =  saturation constant

% Calculate U, assuming Umax and Km
Umax = 0.4;                                                 % units = 1/hr
Km = 40;                                                    % units = nM
C = [1 4 8 12 16 32 64 80 160 800 1600 8000 16000 80000 160000 1600000 8000000];
U = zeros(1,length(C));

for i = 1:length(C)
    U(i) = Umax*C(i)/(Km+C(i));
end

plot(log(C),U)

%clear U Umax C Km;

