%% distributeDivisionTime

% Goal: this script plots the PDF of division times from a given condition


% Strategy:
%           0. for each experiment,
%                  1. initialize data matrix
%                  2. isolate condition of interest
%                  3. isolate time period of interest
%                  4. collect drop data and track ID
%                  5. count number of drops per track
%                        6. if at least two, there's a whole curve:
%                                i. align drop and time data
%                               ii. measure time between birth and division
%                              iii. save in div time array
%                           else, move onto next track
%                  7. build PDF from completed div time array
%              repeat for any other conditions of interest
%           8. compute mean and standard deviation, plot fitted gaussian on pdf


% last edit: jen, 2017 Jun 16
% commit: retired, 2018 Jan 25

% OK LEZ GO!

%%
% 0. initialize experiments for analysis

clear
experiments{1} = '2017-01-06';
experiments{2} = '2017-01-16';
experiments{3} = '2017-01-17';
experiments{4} = '2017-01-18';
experiments{5} = '2017-02-10';
experiments{6} = '2017-02-11';

%%
divisionTimes = [];
numCellCycles = 0; % counter

for i = 1:length(experiments)
    
    clearvars -except i experiments divisionTimes numCellCycles;
    
    % 1. initialize data matrix
    newFolder = strcat('/Users/jen/Documents/StockerLab/Data/',experiments{i},'  (t300)');
    cd(newFolder);
    load(['dm-t300-',experiments{i},'.mat']);
    
    load('meta.mat');
    meta = importdata('meta.mat');
    

    
    % 2. isolate condition of interest
    c = 3;
    matrix_dims = size(dataMatrix);
    trimmedData = dataMatrix(dataMatrix(:,matrix_dims(2)) == c,:);
    
    
    
    % 3. isolate time period of interest
    minTime = meta(c,3); % based on where growth rate seems to stabilize
    maxTime = meta(c,4);
    
    % remove data from timepoints earlier than minTime
    trimmedData2 = trimmedData(trimmedData(:,2) >= minTime,:);
    
    % remove data from timepoints later than maxTime
    trimmedData_final = trimmedData2(trimmedData2(:,2) <= maxTime,:);
    
    
    
    % 4. collect drop data and track ID
    drops = trimmedData_final(:,5);
    trackIDs = trimmedData_final(:,1);
    time = trimmedData_final(:,2);
    
    
    
    % 5. count number of drops per track
    
    % associate track ID to each drop
    dropIDs = trackIDs(drops == 1);
    
    % count drops per ID
    unique_dropIDs = unique(dropIDs);
    drop_counts_per_ID = histc(dropIDs,unique_dropIDs);
    
    % IDs with at least 2 drops
    atLeast2 = find(drop_counts_per_ID == 2);
    
    
    
    
    % 6. for each track with at least two drops
    for tr = 1:length(atLeast2)
        
        % isolate drop data
        tr_drops = drops(trackIDs == unique_dropIDs( atLeast2(tr) ));
        
        % get timestamp at all drops (when drop == 1)
        tr_time = time(trackIDs == unique_dropIDs( atLeast2(tr) ));
        drop_times = tr_time(tr_drops == 1);
        
        % calculate time duration between drops
        tr_divTimes = diff(drop_times);
        
        % concatenate array of all div times per condition
        divisionTimes = [divisionTimes; tr_divTimes]; %in hours
   
    end
    numCellCycles = length(divisionTimes)
    
    
    % 7. build PDF from completed div time array
    inMinutes = divisionTimes*60;
    
    figure(i)
    histogram(inMinutes,'Normalization','pdf','BinWidth',10);
    xlabel('Division Time (min)')
    ylabel('PDF')
    legend('0.5 uM glucose');
    
end

%%
% 8. compute mean and standard deviation for fitted plot on final PDF

clearvars -except i experiments divisionTimes inMinutes numCellCycles;

% trim values below 40 mins
inMinutes_trimmed = inMinutes(inMinutes >= 40);

% fit normal distribution to division time data
normalFit = fitdist(inMinutes_trimmed,'Normal');

% compute PDF from data, mean, and standard deviation
theo_divTimes = linspace(0,300,301);
computed_pdf = normpdf(theo_divTimes, normalFit.mu, normalFit.sigma);

% plot measured and computed distribution
figure()
histogram(inMinutes_trimmed,'Normalization','pdf','BinWidth',10);
xlabel('Division Time (min)')
ylabel('PDF')
hold on
plot(computed_pdf,'r','linewidth',2)
axis([0,300,0,0.015])
legend('measured','normal fit');


