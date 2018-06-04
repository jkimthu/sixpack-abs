% figure 32

% goal: collapse birth size normalized plots of added mass over binned birth size,
%       which mirror Fig 2a of Taheri-Araghi et al. Current Biology (2014)

%       quite like figure 28 & 31 save for the collapse. plots full and
%       constrained range of birth length. uses data structure of collected
%       cell cycle stats build in first part of figure 28.m script.

%       four output plots:
%           1. added vol/mean birth_vol vs. birth_vol/mean birth_vol,
%              broader x axis range
%           2. added vol/mean birth_vol vs. birth_vol/mean birth_vol,
%              Taheri x axis range
%           3. added length/mean birth_length vs. birth_length/mean birth_length,
%              broader x axis range
%           2. added length/mean birth_length vs. birth_length/mean birth_length,
%              Taheri x axis range


% strategy:
%
%       0. initialize data from structures built from figure 28.m
%          this consists of six structures:
%               - time, birth & division
%               - length, at birth & added
%               - volume, at birth & added

%           each with seven columns representing a different environmental condition
%           in this order: low; 30; 300; 900; 3600; ave; high

%           this data consists of data from cell cycles born after 3 hr and
%           divinding before bubble occurrence (if any)

%       0. initialize birth size ranges, as used in Taheri-Araghi et al (2015)
%       0. assign my environment conditions to their Taheri et al analogs,
%          based on similarity in inter-division time
%       1. for each environment...
%               2. isolate column and concatenate data
%               3. remove unphysical events, where added size is negative
%               4. remove non-single cell cycles, where division size > 3x birth size
%               5. normalize all values by mean of birth size
%               6. bin birth volumes by every 20th of mean birth size
%               7. accumulate added volume by birth bin 
%               8. determine x axis scale
%               9. plot all conditions overlaid
%      10. repeat for all!


% last update: jen, 2018 June 4

% commit: plot collapsed added mass over birth size, each normalized by the
%         population-average birth size


% OK let's go!

%%

clc
clear

% 0. initialize data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
load('addedSizeData.mat')

palette = {'Indigo','DarkTurquoise','SteelBlue','DeepSkyBlue','DodgerBlue','GoldenRod','FireBrick'};
environment_order = {'low','30','300','900','3600','ave','high'};


% 0. initialize birth size ranges, as used in Taheri-Araghi et al (2015)
%    from the article, I took the mean birth length per condition and
%    calculated the lower and upper range in birth length as a percentage
%    from the mean.

ranges = [24 38; 21 24; 21 27]; % as a percentage of the mean birth length
                                % left column is lower range
                                % rows are glycerol, glucose + 12aa, synthetic rich media
                                % with interdivision times: 51.3, 26.7, and 22.5 min
                                
% 0. assign my environment conditions to their Taheri et al analogs, based
%    on similarity in interdivision time
eRange = [1; 2; 2; 1; 1; 2; 3]; 
                                

% 1. for each environment... 
for eCol = 1:length(environment_order)
    
    % 2. isolate column and concatenate data
    withData = cellfun(@isempty,V_birth(:,eCol));
    
    % volume
    currCol_birthvols = V_birth(withData==0,eCol);
    currCol_addedvols = V_added(withData==0,eCol);
    concat_Vbirth = [];
    concat_Vadded = [];
    
    % length
    currCol_birthlengths = L_birth(withData==0,eCol);
    currCol_addedlengths = L_added(withData==0,eCol);
    concat_Lbirth = [];
    concat_Ladded = [];
    
    % times
    currCol_birthtimes = T_birth(withData==0,eCol);
    currCol_divtimes = T_division(withData==0,eCol);
    concat_Tbirth = [];
    concat_Tdiv = [];
    
    for replicate = 1:length(currCol_birthvols)
        
        % volume
        rep_Vbirths = currCol_birthvols{replicate};
        rep_Vadded = currCol_addedvols{replicate};
        concat_Vbirth = [concat_Vbirth; rep_Vbirths];
        concat_Vadded = [concat_Vadded; rep_Vadded];
        
        % length
        rep_Lbirths = currCol_birthlengths{replicate};
        rep_Ladded = currCol_addedlengths{replicate};
        concat_Lbirth = [concat_Lbirth; rep_Lbirths];
        concat_Ladded = [concat_Ladded; rep_Ladded];
        
        % times
        rep_Tbirths = currCol_birthtimes{replicate};
        rep_Tdivs = currCol_divtimes{replicate};
        concat_Tbirth = [concat_Tbirth; rep_Tbirths];
        concat_Tdiv = [concat_Tdiv; rep_Tdivs];
        
    end
    
    
    % 3. remove unphysical events, where added size is negative
    inconceivables_vol = concat_Vadded(concat_Vadded < 0);
    inconceivables_length = concat_Ladded(concat_Ladded < 0);
    
    conceivables_Vbirth = concat_Vbirth(concat_Vadded > 0);
    conceivables_Vadded = concat_Vadded(concat_Vadded > 0);
    
    conceivables_Lbirth = concat_Lbirth(concat_Ladded > 0);
    conceivables_Ladded = concat_Ladded(concat_Ladded > 0);
    
    
    % 4. remove non-single cell cycles, where division size > 3* birth size
    division_volume = conceivables_Vbirth + conceivables_Vadded;
    division_length = conceivables_Lbirth + conceivables_Ladded;
    
    vBirth_3x = conceivables_Vbirth * 3;
    lBirth_3x = conceivables_Lbirth * 3;
    
    singles_Vbirth = conceivables_Vbirth(division_volume <= vBirth_3x);
    singles_Vadded = conceivables_Vadded(division_volume <= vBirth_3x);
    
    singles_Lbirth = conceivables_Lbirth(division_length <= lBirth_3x);
    singles_Ladded = conceivables_Ladded(division_length <= lBirth_3x);
    
    
    % 5. normalize all values by mean of birth size
    mean_Vbirth = mean(singles_Vbirth);
    mean_Lbirth = mean(singles_Lbirth);
    
    normalized_Vbirth = singles_Vbirth./mean_Vbirth;
    normalized_Vadded = singles_Vadded./mean_Vbirth;
    
    normalized_Lbirth = singles_Lbirth./mean_Lbirth;
    normalized_Ladded = singles_Ladded./mean_Lbirth;
    
    
    % 6. bin birth volumes by every 0.1 um or cubic um OR 20th of mean birth size    
    fractionOfRaw = 10;
    fractionsOfMean = 20;
    
    % raw values
    birthBins_vol = floor(conceivables_Vbirth*fractionOfRaw);
    bins = 1:max(birthBins_vol);
    xbins_vol = bins'/ fractionOfRaw;

    birthBins_length = floor(conceivables_Lbirth*fractionOfRaw);
    bins = 1:max(birthBins_length);
    xbins_length = bins'/fractionOfRaw;

    % normalized by mean birth size
    birthBins_vol_normalized = floor(normalized_Vbirth*fractionsOfMean);
    bins_normalized = 1:max(birthBins_vol_normalized);
    xbins_vol_normalized = bins_normalized'/ fractionsOfMean;
    
    birthBins_length_normalized = floor(normalized_Lbirth*fractionsOfMean);
    bins_normalized = 1:max(birthBins_length_normalized);
    xbins_length_normalized = bins_normalized'/fractionsOfMean;
    
    
    % 7. accumulate added volume by birth bin 
    
    % raw values
    binned_addedVols = accumarray(birthBins_vol,conceivables_Vadded,[],@(x) {x});
    binned_v_means = cellfun(@mean,binned_addedVols);
    binned_v_stds = cellfun(@std,binned_addedVols);
    binned_v_counts = cellfun(@length,binned_addedVols);
    binned_v_sems = binned_v_stds./sqrt(binned_v_counts);
    
    binned_addedLengths = accumarray(birthBins_length,conceivables_Ladded,[],@(x) {x});
    binned_l_means = cellfun(@mean,binned_addedLengths);
    binned_l_stds = cellfun(@std,binned_addedLengths);
    binned_l_counts = cellfun(@length,binned_addedLengths);
    binned_l_sems = binned_l_stds./sqrt(binned_l_counts);
    
    % normalized by mean birth size
    binned_addedVols_normalized = accumarray(birthBins_vol_normalized,normalized_Vadded,[],@(x) {x});
    binned_v_norm_means = cellfun(@mean,binned_addedVols_normalized);
    binned_v_norm_stds = cellfun(@std,binned_addedVols_normalized);
    binned_v_norm_counts = cellfun(@length,binned_addedVols_normalized);
    binned_v_norm_sems = binned_v_norm_stds./sqrt(binned_v_norm_counts);

    binned_addedLengths_normalized = accumarray(birthBins_length_normalized,normalized_Ladded,[],@(x) {x});
    binned_l_norm_means = cellfun(@mean,binned_addedLengths_normalized);
    binned_l_norm_stds = cellfun(@std,binned_addedLengths_normalized);
    binned_l_norm_counts = cellfun(@length,binned_addedLengths_normalized);
    binned_l_norm_sems = binned_l_norm_stds./sqrt(binned_l_norm_counts);
    
    
    % 8. determine x axis scale
    percentz = ranges(eRange(eCol),:);
    lowerRange_raw_vol = mean_Vbirth - (percentz(1)/100)*mean_Vbirth;
    upperRange_raw_vol = mean_Vbirth + (percentz(2)/100)*mean_Vbirth;
    
    lowerRange_raw_l = mean_Lbirth - (percentz(1)/100)*mean_Lbirth;
    upperRange_raw_l = mean_Lbirth + (percentz(2)/100)*mean_Lbirth;
    
    lowerRange_norm = 1 - percentz(1)/100;
    upperRange_norm = 1 + percentz(2)/100;
    
    
    % 9. plot subplots for all conditions
    color = rgb(palette{eCol});
    condition = environment_order{eCol};


    % volume, collasped larger range
    figure(1)
    plot(xbins_vol_normalized,binned_v_norm_means,'o','Color',color)
    legend('low','30 sec','5 min','15 min','60 min','ave','high')
    hold on
    errorbar(xbins_vol_normalized,binned_v_norm_means,binned_v_norm_sems,'.','Color',color)
    axis([0 5 -.1 3])
    title('collapsed added vol vs birth vol: larger range')
    ylabel('added volume / mean birth volume')
    xlabel('birth volume / mean birth volume')
    
    
    % volume, collasped Taheri range
    figure(2)
    plot(xbins_vol_normalized,binned_v_norm_means,'o','Color',color)
    legend('low','30 sec','5 min','15 min','60 min','ave','high')
    hold on
    errorbar(xbins_vol_normalized,binned_v_norm_means,binned_v_norm_sems,'.','Color',color)
    axis([lowerRange_norm upperRange_norm -.1 2])
    title('collapsed added vol vs birth vol: Taheri range')
    ylabel('added volume / mean birth volume')
    xlabel('birth volume / mean birth volume')

    
    
    % length, collasped larger range
    figure(3)
    plot(xbins_length_normalized,binned_l_norm_means,'o','Color',color)
    legend('low','30 sec','5 min','15 min','60 min','ave','high')
    hold on
    errorbar(xbins_length_normalized,binned_l_norm_means,binned_l_norm_sems,'.','Color',color)
    axis([0 5 -1 3])
    title('collapsed added length vs birth length: larger range')
    ylabel('added length, mean + sem')
    xlabel('birth length, mean + sem')

    
    
    % length, collasped Taheri range
    figure(4)
    plot(xbins_length_normalized,binned_l_norm_means,'o','Color',color)
    legend('low','30 sec','5 min','15 min','60 min','ave','high')
    hold on
    errorbar(xbins_length_normalized,binned_l_norm_means,binned_l_norm_sems,'.','Color',color)
    axis([lowerRange_norm upperRange_norm -.1 2])
    title('collapsed added length vs birth length: Taheri range')
    ylabel('added length, mean + sem')
    xlabel('birth length, mean + sem')
    
end

