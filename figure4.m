%% figures 4


%  Goal: plot a heated scatter of single-cell division size vs birth size.
%        consider only cell cycles born after 3 hrs


%  Strategy: 
%
%       a) initialize experimental data
%       b) identify complete cell cycles within each condition 
%       c) compile birth volumes, division volumes and time at birth
%       d) plot raw condition data and heatmap version, for each cycle
%       e) save both figures per experiment



%  Last edit: jen, 2019 April 30
%  Commit: first commit, full and heatmap scatter plots for each condition



%  OK let's go!

%% initialize

clc
clear

% 0. initialize complete meta data
%cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')

max_vbirth = 10;
max_vdiv = 20;

%%
% 1. for all experiments in dataset
exptArray = [2,3,4];

for e = 1:length(exptArray)
    
    % 1. collect experiment meta data
    index = exptArray(e);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    expType = storedMetaData{index}.experimentType;
    bubbletime = storedMetaData{index}.bubbletime;
    
    
    % 2. load measured data
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    if strcmp(date,'2017-11-12') == 1
        filename = strcat('lb-fluc-',date,'-width1p4-jiggle-0p5.mat');
    else
        filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
    end
    load(filename,'D5','T');
    
    
    % 3. compile experiment data matrix
    xy_start = min(min(storedMetaData{index}.xys));
    xy_end = max(max(storedMetaData{index}.xys));
    exptData = buildDM(D5, T, xy_start, xy_end,index,expType);
    clear D5 T xy_start xy_end e
    
    
    % 4. initialize colors for plotting
    palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
    environment = {'fluc','low','ave','high'};
    
    for condition = 1:length(environment)
        
        
        % 5. isolate condition specific data
        conditionData = exptData(exptData(:,21) == condition,:);  % col 21 = cond vals
        
        
        % 6. trim data to full cell cycles ONLY
        curveFinder = conditionData(:,5);                         % col 5 = curveFinder, ID of full cell cycles
        conditionData_fullOnly = conditionData(curveFinder > 0,:);
        clear curveFinder
  
        
        % 7. isolate volume
        timestamp_hr = conditionData_fullOnly(:,2)/3600;   % col 2  = raw timestamp (sec converted to hr)
        volume = conditionData_fullOnly(:,11);             % col 11  = volume (Va)
        
        
        % 8. identify unique cell cycles by ID number
        isDrop = conditionData_fullOnly(:,4);              % col 4 = isDrop; 1 = birth, 0 when not
        curveFinder = conditionData_fullOnly(:,5);         % col 5 = curveFinder, ID of full cell cycles
        unique_cc = curveFinder(isDrop == 1);
        birthTimes = timestamp_hr(isDrop == 1);
        
        
        % 9. remove birth times prior to 3 hr
        birthTimes_post3 = birthTimes(birthTimes > 3); 
        unique_cc_post3 = unique_cc(birthTimes > 3); 
        
        
        % 10. remove birth times post bubbles
        birthTimes_final = birthTimes_post3(birthTimes_post3 < bubbletime(condition));
        unique_cc_final = unique_cc_post3(birthTimes_post3 < bubbletime(condition));
        
        
        
        % 11. identify volume at birth and at division for each cell cycle
        clear unique_cc unique_cc_post3 birthTimes birthTimes_post3
        V_division = nan(length(unique_cc_final),1);
        V_birth = nan(length(unique_cc_final),1);
        cc_lengths = nan(length(unique_cc_final),1);
        
        for cc = 1:length(unique_cc_final)
            
            currentVolumes = volume(curveFinder == unique_cc_final(cc));
            cc_lengths(cc,1) = length(currentVolumes);
            V_division(cc,1) = currentVolumes(end);
            V_birth(cc,1) = currentVolumes(1);
            
        end
        clear cc
        
        
        
        % 12. trim cell cycles shorter than 6 timepoints
        V_division_6plus = V_division(cc_lengths > 5);
        V_birth_6plus = V_birth(cc_lengths > 5);
        
        
        
        % 13. trim outliers (those 3 std dev away from median) from final dataset
        %     a) id median and std of division size
        %     b) id median and std of birth size
        %     c) find indeces in both vectors that are within 3 std
        %     d) keep data from indeces in both vectors
        
        divSize_median = median(V_division_6plus); 
        divSize_std = std(V_division_6plus);
        
        birthSize_median = median(V_birth_6plus); 
        birthSize_std = std(V_birth_6plus);
        
        div_bigOutlier = find(V_division_6plus > (divSize_median+divSize_std*3));
        div_smallOutlier = find(V_division_6plus < (divSize_median-divSize_std*3));
        div_outliers = [div_bigOutlier; div_smallOutlier];
        
        birth_bigOutlier = find(V_birth_6plus > (birthSize_median+birthSize_std*3));
        birth_smallOutlier = find(V_birth_6plus < (birthSize_median-birthSize_std*3));
        birth_outliers = [birth_bigOutlier; birth_smallOutlier];

        V_division_binary = ones(length(V_division_6plus),1);
        V_division_binary(div_outliers) = 0;
        
        V_birth_binary = ones(length(V_birth_6plus),1);
        V_birth_binary(birth_outliers) = 0;
        
        V_summed = V_division_binary + V_birth_binary;
        
        V_division_final = V_division_6plus(V_summed == 2);
        V_birth_final = V_birth_6plus(V_summed == 2);
       
        clear birthSize_median birthSize_std divSize_median divSize_std
        clear div_bigOutlier div_smallOutlier birth_bigOutlier birth_smallOutlier
         
        
        
        
        % 14. calculate fit line and calculate correlation coeff
        
        %     i. fit linear regression
        p = polyfit(V_birth_final,V_division_final,1);
        x = V_birth_final;
        y = p(1)*x + p(2);
        
        
        %    ii. calculate correlation coefficient
        r = corrcoef(V_birth_final,V_division_final);
        R1(t) = r(1,2);
        
        
        
        % 14. plot
        color = rgb(palette(condition));
        
%         % for y=2x line
%         x = linspace(0,max_vbirth,10);
%         y = linspace(0,max_vdiv,10);
%         gray = rgb('Silver');
        
        
        % division size vs. birth size
        figure(4)
        subplot(2,2,condition)
        plot(V_birth_final,V_division_final,'o','Color',color)
        hold on
        plot(x,y,'Color',rgb('SlateGray'),'LineWidth',2)
        hold on
        txt = strcat('r=',num2str(r(1,2)));
        text(x(end),y(end),txt,'FontSize',14)
        legend(environment(condition))
        title(date)
        xlabel('birth size (cubic um)')
        ylabel('division size (cubic um)')
        axis([0 max_vbirth 0 max_vdiv])
        
        
        % 15. plot with colorbar indicating density
        binsPerMicron = 2;
        binSize = 1/binsPerMicron; % microns
        
        bin_birth = ceil(V_birth_final/binSize);
        bin_div = ceil(V_division_final/binSize);
        
        bin_mat = zeros(max_vdiv/binSize,max_vbirth/binSize);
        
        test = sub2ind(size(bin_mat),bin_div,bin_birth);
        
        a = unique(test);
        out = [a,histc(test(:),a)];
        
        idx = out(:,1);
        counts = out(:,2);
        
        
        bin_mat(out(:,1)) = out(:,2);
        [i,j] = ind2sub([40,20],idx);
        
        binned_vBirths = j.*binSize;
        binned_vDivs = i.*binSize;
        
        
     
        figure(5)
        subplot(2,2,condition)
        scatter(binned_vBirths,binned_vDivs,60,counts,'filled')
        colormap(parula) % see colormap documentation
        % colorMap = [linspace(0,color(1),max(counts))', linspace(0,color(2),max(counts))', linspace(0,color(3),max(counts))'];
        % colormap(colorMap);
        colorbar;
        plot(x,y,'Color',rgb('SlateGray'),'LineWidth',2)
        hold on
        txt = strcat('r=',num2str(r(1,2)));
        text(x(end),y(end),txt,'FontSize',14)
        axis([0 max_vbirth 0 max_vdiv])
        hold on
        plot(x,y,'Color',gray)
        title(condition)
        xlabel('birth size (cubic um)')
        ylabel('division size (cubic um)')

        


        
    end
    
    % 16. save plots in active folder
    %cd('/Users/jen/Documents/StockerLab/Data_analysis/currentPlots/')
    figure(4)
    plotName = strcat('figure4-div-v-birth-',date,'-noBinning');
    saveas(gcf,plotName,'epsc')
    close(gcf)
    
    figure(5)
    plotName = strcat('figure4-div-v-birth-',date,'-',num2str(binsPerMicron),'binsPerMicron');
    saveas(gcf,plotName,'epsc')
    close(gcf)
    clc
    
end

