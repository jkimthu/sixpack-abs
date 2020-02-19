%% widths

clear
clc


% Load workspace from SlidingFits.m     
load('t300_2017-01-18-Mus-length.mat');
conditions = [1 10; 11 20; 21 30; 31 40];


% average cell width for each condition
for i = 1:4 %number of conditions
    
    %    Condition One    %
    width_cond = [];
    length_cond = [];
    Time_cond = [];
    
    for n = conditions(i,1):conditions(i,2)
        for m = 1:length(D6{n})
            
            %  assemble all lengths and widths into a single vector
            width_cond = [width_cond; D6{n}(m).MinAx(:,1)];
            length_cond = [length_cond; D6{n}(m).MajAx(:,1)];
            
            %  assemble a corresponding timestamp vector
            vectorLength = length(D6{n}(m).MinAx(:,1));
            trackFrames = D6{n}(m).Frame;
            Time_cond = [Time_cond; T{n}(trackFrames)];
            
        end
    end
    clear n m vectorLength
    
    %  convert all timestamps from seconds to hours
    Time_cond = Time_cond/3600;
    
    %  eliminate negative growth rates
    %Mu_cond1(Mu_cond1<0)=NaN;
    
    %  determine size of time bins
    BinsPerHour = 60;                              % multiplying by 10 gives bins of 0.1 hr
    Bins = ceil(Time_cond*BinsPerHour);            % multiplying by 200 gives time bins of 0.005 hr
    %plotUntil = floor(conditions(xy,3)*BinsPerHour);
    
    %  accumulate widths by bin, and calculate mean and std dev
    width_Means = accumarray(Bins,width_cond,[],@nanmean);
    width_STDs = accumarray(Bins,width_cond,[],@nanstd);
    
    %  accumulate lengths by bin, and calculate mean and std dev
    length_Means = accumarray(Bins,length_cond,[],@nanmean);
    length_STDs = accumarray(Bins,length_cond,[],@nanstd);
    
    %   to calculate s.e.m.
    %   1. count number of total tracks in each bin
    for j = 1:max(Bins)                                 % for all bins
        currentBin_count = find(Bins==j);               % find which data belongs to which time bin
        width_Counts(j) = length(currentBin_count);     % count data points in each time bin
        
        clear k counter Kasten;
    end
    
    %   2. divide standard dev by square root of tracks per bin
    width_sems = width_STDs./sqrt(width_Counts');
    length_sems = length_STDs./sqrt(width_Counts');
    
    if i == 1
        color = 'c';
    end 
    
    if i == 2
        color = 'b';
    end
    
    if i == 3
        color = 'y';
    end
    
    if i == 4
        color = 'r';     
    end
    
    errorbar(width_Means,width_sems,color)
    hold on
    errorbar(length_Means,length_sems,color)
    grid on
    axis([0,550,1,5])
    xlabel('Time')
    ylabel('Elongation rate (1/hr)')
    
    
end
legend('fluc', 'low', 'ave', 'high');