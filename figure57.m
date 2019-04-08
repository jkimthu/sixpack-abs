%% figure 57: calculate number of tracks per condition

%  Goals: loop through all experiments and total tracks per condition
%         save a data structure than holds eight values per experiment,
%         tracks in fluc, low, ave and high (rows 1 -4)
%         tracks total (bubble trimmed) and post-3h (columns 1-2)


%  Strategy:
%
%       0. initialize directory and meta data
%       1. create array of experiments of interest, for each:
%               2. initialize experiment meta data
%               3. load measured experiment data    
%               4. for single shift experiments, define which frames to ignore (noisy tracking)
%               5. for each condition in current experiment, build data matrix from specified condition
%                       6. isolate volume (Va), timestamp, drop, curve, and trackNum data
%                       7. calculate growth rate
%                       8. truncate data to non-erroneous (e.g. bubbles) timestamps
%                       9. isolate selected specific growth rate and timestamp
%                      10. if appropriate, assign NaN to all growth rates associated with frames to ignore
%                          else simply remove existing nans from analysis
%                      11. bin growth rate into time bins based on timestamp
%                      12. calculate mean, standard dev, counts, and standard error
%                      13. plot growth rate over time
%              14. save plot with experiment #, specific growth rate definition, and binning          
%      15. repeat for all experiments 


%  last updated: jen, 2019 April 8

%  commit: first commit, taking raw track values from assembled data matrix


% OK let's go!

%% initialize

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));


% 0. initialize new storate structure for track data
trackData = cell(size(storedMetaData));

%% Part 1. loop through experiments to count and record tracks

% 1. create array of experiments of interest, then loop through each:
exptArray = [2,3,4,5,6,7,9:12,13:15,21,22,26,27]; % use corresponding dataIndex values

for e = 1:length(exptArray)
    
    
    % 2. initialize experiment meta data
    index = exptArray(e);
    date = storedMetaData{index}.date;
    expType = storedMetaData{index}.experimentType;
    bubbletime = storedMetaData{index}.bubbletime;
    timescale = storedMetaData{index}.timescale;
    disp(strcat(date, ': analyze!'))
    
    
    
    % 3. load measured experiment data    
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    if strcmp(date,'2017-11-09') == 1
        filename = strcat('lb-control-',date,'-width1p4-jiggle-0p5.mat');
    elseif strcmp(date,'2017-09-26') == 1
        filename = 'lb-monod-2017-09-26-jiggle-c12-0p1-c3456-0p5-bigger1p8.mat';
    elseif strcmp(date,'2018-12-04') == 1
        filename = 'lb-monod-2018-12-04-c12-width1p7-c34-width1p4-jiggle-0p5.mat';
    elseif strcmp(expType,'origFluc') == 1
        filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
    elseif strcmp(expType,'fluc2stable') == 1
        filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
    else
        filename = strcat('lb-fluc-',date,'-width1p7-jiggle-0p5.mat');
        % single upshift and downshift data only uses larger width thresh
    end
    load(filename,'D5','T');
    
    
    
    % 4. for single shift experiments, define which frames to ignore (noisy tracking)
    %    these frames to be ignored only apply to condition 1.
    if strcmp(date,'2018-06-15') == 1
        ignoredFrames = [112,113,114];
    elseif strcmp(date,'2018-08-01') == 1
        ignoredFrames = [94,95];
    elseif strcmp(date,'2018-08-09') == 1
        ignoredFrames = [115,116,117];
    elseif strcmp(date,'2018-08-08') == 1
        ignoredFrames = [112,113,114];
    else
        ignoredFrames = [];
    end
    
    
    
    % 5. build data matrix from specified condition
    for condition = 1:length(bubbletime)
        
        
        xy_start = storedMetaData{index}.xys(condition,1);
        xy_end = storedMetaData{index}.xys(condition,end);
        conditionData = buildDM(D5, T, xy_start, xy_end,index,expType);
        
        
        
        % 6. isolate timestamp data
        timestamps_sec = getGrowthParameter(conditionData,'timestamp');   % ND2 file timestamp in seconds
        
        
        
        % 7. truncate data to non-erroneous (e.g. bubbles) timestamps
        maxTime = bubbletime(condition);
        timestamps_hr = timestamps_sec/3600; % time in seconds converted to hours
        
        if maxTime > 0
            conditionData_bubbleTrimmed = conditionData(timestamps_hr <= maxTime,:);
            timestamps_bubbleTrimmed = timestamps_hr(timestamps_hr <= maxTime,:);
        else
            conditionData_bubbleTrimmed = conditionData;
            timestamps_bubbleTrimmed = timestamps_hr;
        end
        
        
        
        % 8. remove data prior to 3h
        minTime = 3;
        post3Data = conditionData_bubbleTrimmed(timestamps_bubbleTrimmed >= minTime,:);
        clear timestamps_hr timestamps_sec timestamps_bubbleTrimmed


        
        % 9. isolate track number
        trackNum_total = getGrowthParameter(conditionData_bubbleTrimmed,'trackNum');          % track number, not ID from particle tracking
        trackNum_post3 = getGrowthParameter(post3Data,'trackNum'); 
        clear conditionData 
        
        
        
        % 10. count unique tracks in condition
        uniqs_total = unique(trackNum_total);
        uniqs_post3 = unique(trackNum_post3);
        clear trackNum_total trackNum_post3
        
        
        
        % 11. compile experiment data
        track_counts(condition,1) = length(uniqs_total); % row = condition
        track_counts(condition,2) = length(uniqs_post3); % column 1 = total tracks; column 2 = post 3h data
        clear uniqs_total uniqs_post3

        
    end
    
    % 12. store data in structure
    trackData{exptArray(e)} = track_counts;
    clear track_counts D5 T
    
end

save('trackData.mat','trackData')

%% Part 2. plot track counts to visualize full data

% goal: plot scatter of each condition count as total vs post 3 count


% 0. initialize track data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('trackData.mat')


% 0. initialize experiment array
exptArray = [2,3,4,5,6,7,9:12,13:15,21,22,26,27];


% 0. initialize colors for plotting
palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};


% 1. loop through experiments and plot a single point per row, colored by condition 
fluc_total = [];
fluc_post = [];
low_total = [];
low_post = [];
ave_total = [];
ave_post = [];
high_total = [];
high_post = [];

for e = 1:length(exptArray)
    
    currentData = trackData{exptArray(e)};
    
    for c = 1:4
        
        color = rgb(palette(c));
        rowData = currentData(c,:);
        
        % plot scatter
        figure(1)
        hold on
        if e > 13
            plot(rowData(1),rowData(2),'x','Color',color,'MarkerSize',7)
        else
            plot(rowData(1),rowData(2),'o','Color',color,'MarkerSize',7)
        end
        
        
        % concatenate condition data for distributions
        if c == 1
            fluc_total = [fluc_total; rowData(1)];
            fluc_post = [fluc_post; rowData(2)];
        elseif c == 2
            low_total = [low_total; rowData(1)];
            low_post = [low_post; rowData(2)];
        elseif c == 3
            ave_total = [ave_total; rowData(1)];
            ave_post = [ave_post; rowData(2)];
        elseif c ==4
            high_total = [high_total; rowData(1)];
            high_post = [high_post; rowData(2)];
        end
        
    end
    figure(1)
    xlabel('tracks after 3h')
    ylabel('total tracks')

end


