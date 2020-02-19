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
load('2016-08-05_b-Mus-length.mat');

counter =0;
for n = 1:10:40
    counter = counter +1;
    m = 1;
    
    % Extracted mu
    Mu_track = M6{n}(m).Parameters(:,1);
    vectorLength = length(Mu_track);
    
    % Original length data (microns)
    Ltrack2 = D6{n}(m).MajAx(3:vectorLength+2);                                  % trimmed length trajectory to match Mu
    
    % Time data (hours)
    %dT = mean(mean(diff(T)));                                              % mean time between frames (seconds)
    %Ttrack = D6{n}(m).Frame(3:Num_mu+2);                                   % original frame # in trajectory
    %timeTrack = T{n}/(60*60);
    timeTrack = T{n}/(3600);
    
    figure(1)
    
    subplot(4,1,counter)
    plot(timeTrack(3:vectorLength+2),Ltrack2,'.',timeTrack(3:vectorLength+2),Mu_track*log(2),'r.');                          
    grid on;
    axis([0,10,-0.5,6])
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
%clear;
load('2015-08-10-Mus-length.mat','D6','M6','T');

% defining conditions: col1 = first xy; col2 = final xy; col3 = time (hr) cutoff
conditions = [1 10 10; 11 20 10];% 21 30 4; 31 40 4];
%%

for i = 1:2 %number of conditions
    
    %    Condition One    %
    Mu_cond = [];
    Time_cond = [];
    
    for n = conditions(i,1):conditions(i,2)
        for m = 1:length(M6{n})
            
            %  assemble all instantaneous growth rates into a single vector
            Mu_cond = [Mu_cond; M6{n}(m).Parameters(:,1)];
            
            %  assemble a corresponding timestamp vector
            vectorLength = length(M6{n}(m).Parameters(:,1));
            trackFrames = D6{n}(m).Frame(3:vectorLength+2);
            %Time_cond = [Time_cond; T{n}(trackFrames)];
            Time_cond = [Time_cond; T(trackFrames,n)];
        end
    end
    
    %  convert all timestamps from seconds to hours
    Time_cond = Time_cond/3600;
    
    %  eliminate negative growth rates
    Mu_cond= Mu_cond(Mu_cond > 0);
    Time_cond = Time_cond(Mu_cond > 0);
    
    %  determine size of time bins
    BinsPerHour = 10;                              % multiplying by 10 gives bins of 0.1 hr
    Bins = ceil(Time_cond*BinsPerHour);            % multiplying by 200 gives time bins of 0.005 hr
    %plotUntil = floor(conditions(xy,3)*BinsPerHour);
    
    %  accumulate growth rates by bin, and calculate mean and std dev
    Mu_Means = accumarray(Bins,Mu_cond,[],@nanmean);
    Mu_STDs = accumarray(Bins,Mu_cond,[],@nanstd);
    
    
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
    
    %  eliminate bins without data
    bin_cond = linspace(1,length(Mu_sems),length(Mu_sems))';
    bin_cond2 = bin_cond(Mu_sems > 0);
    Mu_sems2= Mu_sems(Mu_sems > 0);
    Mu_Means2 = Mu_Means(Mu_sems > 0);
    
    errorbar(bin_cond2,Mu_Means2,Mu_sems2)
    %errorbar( Mu_Means(1:plotUntil),Mu_sems(1:plotUntil) )
    hold on
    grid on
    axis([0,200,0.3,.5])
    xlabel('Time')
    ylabel('Elongation rate (1/hr)')
    %forLegend = num2str(xy);
    %legend(forLegend)
    
    clear vectorLength trackFrams Mu_Means Mu_STDs Mu_sems Bins hr dT Mu_Counts n m j;
    clear Mu_cond Time_cond plotUntil;
    
end
legend('const', 'fluc')
%legend('fluc', 'low', 'ave', 'high');
%%
