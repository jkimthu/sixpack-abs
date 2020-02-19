%% distribute( dataz )

%  Goal: plot distributions of cell cycle duration and added mass,
%        normalized by population average
%
%  Goal: plot distribution of cell size at birth



%  Last edit: Jen Nguyen, April 26 2017


%  Section contents:
%  >> sections are separated based on input data formats
%
%       1. Cell cycle duration and added mass


% OK! Lez go!



%%  T W O.
%   plot distribution of birth size


% The intended input for these scripts is the following data matrix,
% saved with the naming convention of:

% dmMMDD-cond.mat

%      where,
%              dm  =  dataMatrix                  (see matrixBuilder.m)
%              MM  =  month of experimental date
%              DD  =  day of experimental date
%       condition  =  experimental condition      (fluc or const)
%



% Initialize data.
clear
load('dm-t900-2017-01-10.mat');
load('meta.mat');
meta = meta_2017jan10;


%%
%  Stragety:
%
%   For all for conditions:
%      0. designate time window of analysis
%      1. isolate data of interest (length and drop)
%      2. find length when drop == 1
%      3. plot!


for condition = 1:4
    
    % 0. designate time window of analysis
    
    firstTimepoint = meta(condition,3); % in hours
    lastTimepoint = meta(condition,4);
    
    % 1. isolate data from condition of interest 
    interestingData = find(dataMatrix(:,19) == condition);
    conditionData = dataMatrix(interestingData,:);
    
    % 1. isolate Length, Vc and Drop data
    data1 = conditionData(:,17); % vcAdded
    %data2 = conditionData(:,15); % mu_vc
    %veVals = conditionData(:,18); % veAdded
    durations = conditionData(:,8);
    %lengthVals = conditionData(:,3);
    drop = conditionData(:,5);
    %volVals = conditionData(:,13);
    timeStamps = conditionData(:,2);

    % 0. trim off timepoints earlier than first
    %data2 = data2(timeStamps >= firstTimepoint);
    drop = drop(timeStamps >= firstTimepoint);
    data1 = data1(timeStamps >= firstTimepoint);
    durations = durations(timeStamps >= firstTimepoint);
    lowTrimmed_timeStamps = timeStamps(timeStamps >= firstTimepoint);
    
    % 0. trim off timepoints later than last
    %data2 = data2(lowTrimmed_timeStamps <= lastTimepoint);
    drop = drop(lowTrimmed_timeStamps <= lastTimepoint);
    data1 = data1(lowTrimmed_timeStamps <= lastTimepoint);
    durations = durations(lowTrimmed_timeStamps <= lastTimepoint);
    finalTrimmed_timeStamps = lowTrimmed_timeStamps(lowTrimmed_timeStamps <= lastTimepoint);
    
    % 0. trim off timepoints without a drop
    %data2 = data2(drop == 1);
    data1 = data1(drop == 1);
    durations = durations(drop == 1);
    finalTrimmed_timeStamps = finalTrimmed_timeStamps(drop == 1);
    
    % 2. keep lengths when drop equals 1 (denotes birth)
    %compiled_data2{:,condition} = data2; %./durations;
    compiled_data1{:,condition} = data1./durations;
    compiled_time{:,condition} = finalTrimmed_timeStamps;
    
     % for cell cycle related measures, remove zeros from data sets.
     for i = 1:length(compiled_data1)
         compiled_data1{i}(compiled_data1{i} <= 0) = NaN;
         nanFilter = find(~isnan(compiled_data1{i}));
         compiled_data1{i} = compiled_data1{i}(nanFilter);
         compiled_time{i} = compiled_time{i}(nanFilter);
     end
     
%      for i = 1:length(compiled_data2)
%          compiled_data2{i}(compiled_data2{i} <= 0) = NaN;
%          nanFilter = find(~isnan(compiled_data2{i}));
%          compiled_data2{i} = compiled_data2{i}(nanFilter);
%      end
     clear i
     
      % for cell cycle related measures, remove super large outliers from data sets.
%      for i = 1:length(compiled_data1)
%          compiled_data1{i}(compiled_data1{i} >= 5) = NaN;
%          nanFilter = find(~isnan(compiled_data1{i}));
%          compiled_data1{i} = compiled_data1{i}(nanFilter);
%      end
     
%      for i = 1:length(compiled_data2)
%          compiled_data2{i}(compiled_data2{i} >= 5) = NaN;
%          nanFilter = find(~isnan(compiled_data2{i}));
%          compiled_data2{i} = compiled_data2{i}(nanFilter);
%      end
     
    % 3. plot
    figure(1)
    %distributionPlot(compiled_data2,'widthDiv',[2 1],'histOri','left','color','b','showMM',2)
    %distributionPlot(gca,compiled_data1,'widthDiv',[2 2],'histOri','right','color',[0,1,1],'showMM',2)
    subplot(4,1,condition)
    plot(compiled_time{condition}, compiled_data1{condition},'.')
    xlabel('Time')
    ylabel('V_delta / tau_D (1/hr)')
    hold on
    axis([2,8,-0.1,7])
    %histogram(birthLengths,'BinWidth',0.1)
    %hold on
    
    
    biomassProd_Means(condition) = mean(compiled_data1{condition});
    biomassProd_STDs(condition) = std(compiled_data1{condition});
    biomassProd_Counts(condition) = length(compiled_data1{condition});
    biomassProd_SEMs(condition) = biomassProd_STDs(condition)/sqrt( biomassProd_Counts(condition) );
    
    %mu_VE_sems = mu_VE_STDs./sqrt(Mu_Counts');
 
    figure(2)
    %plot(condition,biomassProd_SEMs(condition),condition,biomassProd_Means(condition),'o')
    errorbar(condition,biomassProd_Means(condition),biomassProd_SEMs(condition))
    hold on
    plot(condition,biomassProd_Means(condition),'o')
    axis([0,5,0,2])
    xlabel('Condition')
    ylabel('Mean biomass production rate (1/hr)')



end
%% testing distributionPlot.m

% functions from distributionPlot.m files, as shared by Jonas Dorn on File
% Exchange. 

data1 = randn(500,5);
data2 = bsxfun(@plus,randn(500,5),0:0.1:0.4);
figure
distributionPlot(data1,'widthDiv',[2 1],'histOri','left','color','b','showMM',4)
distributionPlot(gca,data2,'widthDiv',[2 2],'histOri','right','color','k','showMM',4)

% bsxfun(fun, A, B) : applies the function (fun) to arrays A and B
%                       ex. C = bsxfun(@minus, A, mean(A))
%                           subtracts the mean of A column-wise from each
%                           element in corresponding columns of A

% gca : helps return to previous plot

% widthDiv : [numberOfDivisions,currentDivision], allows comparison of
%              multiple distributions

% histOri : orientation of histogram ('center','left', or 'right'), with 'center' as default
%             'left' or 'right' only shows left or right half of violin plot

% color : uniform coloring of histograms. Supply either a color
%           string ('r'), or a truecolor vector ([1 0 0]). Use a
%           cell array of length nData to specify one color per
%           distribution. Default: 'k' 

% showMM : if 1, mean and median are shown as red crosses and
%                green squares, respectively. This is the default
%                2: only mean
%                3: only median
%                4: mean +/- standard error of the mean (no median)
%                5: mean +/- standard deviation (no median)
%                6: draw lines at the 25,50,75 percentiles (no mean)
%                0: plot neither mean nor median

