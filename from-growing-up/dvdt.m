% dvdt

% goal: given a data matrix for a specific condition,
%       calculate and report dV/dt values in a three column array:

%       col 1  =  dV/dt values
%       col 2  =  nutrient signal
%       col 3  =  difference in nutrient signal
%                      1 = upshift occurred between current timepoint and one before
%                     -1 = downshift occurred between current timepoint and one before
%                      0 = no change during between timepoint and one before


%       intended use: 
%
%                (1)  array provides info needed to determine growth rates
%                     specific to a LB dilution, without smoothing
%
%                (2)  for use in combination with data matrix, as indeces
%                     are kept consistent


% strategy: 
%       
%       0. initialize corrected timestamp and volume data
%       0. check that conditionData contains only one condition.
%          if more, report error
%       1. for each track, calculate dVdt reporting NaN at drops
%               2. isolate full curves and corresponding timestamps (lag corrected)
%               3. determine whether track contains any full curves
%               4. if no full curves, score for track = NaN
%               5. if full curves are present, for each curve
%               6. isolate volumes for current curve
%               7. calculate change in volume from previous timestep
%               8. translate timestamps into quarters of nutrient signal
%               9. from nutrient signal quarters, generate a binary nutrient signal where, 1 = high and 0 = low
%              10. determine whether each timestep was a step up or down
%              11. save calculated growth rate and meta data
%              12. concatenate scores for all tracks
%      13. output dvdt and nutrient meta data


% last updated: jen, 2018 Apr 30

% commit: edit to allow dV/dt data to be calculated for expts lacking
%         corrected time data (i.e. 2017-10-10)


% OK let's go!


function [dvdt_data] = dvdt(conditionData, timescale, date)

% 0. initialize corrected timestamp and volume data
correctedTime = conditionData(:,25);           % col 25 = timestamps, reflecing true time for all conditions in sec
rawTime = conditionData(:,2);                  % col 2  = raw timestamps
volume = conditionData(:,12);                  % col 12 = volumes (va)
curveID = conditionData(:,6);                  % col 6  = curve ID per track
trackNum = conditionData(:,22);                % col 22 = track number, not ID from particle tracking
condVals = conditionData(:,23);                % col 23 = condition number


% 0. check that conditionData contains only one condition
%    if more, report error
if length(unique(condVals)) > 1
    
    error('ERROR: dvdt function requires that conditionData contain only ONE condition')
    
end

% 1. calculate dvdt for each track, reporting NaN at drops
dvdt_data = [];

for tr = 1:max(trackNum)
    
    % 2. isolate full curves and corresponding timestamps (lag corrected)
    currentTrack_curves = curveID(trackNum == tr);
    currentTrack_correctedTimes = correctedTime(trackNum == tr);
    currentTrack_rawTimes = rawTime(trackNum == tr);
    currentTrack_volumes = volume(trackNum == tr);

    
    % 3. determine whether track contains any full curves
    if sum(currentTrack_curves) == 0
        
        % 4. if no full curves, score for track = NaN
        track_dv = NaN(length(currentTrack_curves),3);
        
    else
        
        % 5. if full curves are present, for each curve
        track_dv = NaN(length(currentTrack_curves),3);
        
        for cc = 1:max(currentTrack_curves)
            
            % 6. isolate volumes for current curve
            currentCurve_indeces = find(currentTrack_curves == cc);
            currentCurve_volumes = currentTrack_volumes(currentCurve_indeces);
            currentCurve_correctedTimes = currentTrack_correctedTimes(currentCurve_indeces);
            currentCurve_rawTimes = currentTrack_rawTimes(currentCurve_indeces);
            
            % 7. calculate change in volume from previous timestep
            dV = [NaN; diff(currentCurve_volumes)];
            
            if strcmp(date, '2017-10-10') == 1
                %disp(strcat(date,': calculate dV/dt with raw timestamp'))
                dt = [NaN; diff(currentCurve_rawTimes)];      % timestep in seconds
            else
                dt = [NaN; diff(currentCurve_correctedTimes)];      % timestep in seconds
            end
            dVdt = dV./dt * 3600;

            
            if unique(condVals) == 1
                
                % 8. translate timestamps into quarters of nutrient signal
                timeInPeriods = currentCurve_correctedTimes/timescale; % unit = sec/sec
                timeInPeriodFraction = timeInPeriods - floor(timeInPeriods);
                timeInQuarters = ceil(timeInPeriodFraction * 4);
                
                % 9. from nutrient signal quarters, generate a binary nutrient signal where, 1 = high and 0 = low
                binaryNutrientSignal = zeros(length(timeInQuarters),1);
                binaryNutrientSignal(timeInQuarters == 1) = 1;
                binaryNutrientSignal(timeInQuarters == 4) = 1;
                
                % 10. determine whether each timestep was a step up or down
                isShift = [NaN; diff(binaryNutrientSignal)];
                
                
            elseif unique(condVals) == 2
                
                binaryNutrientSignal = zeros(length(currentCurve_correctedTimes),1);
                isShift = [NaN; diff(binaryNutrientSignal)];
                
            elseif unique(condVals) == 3
                
                binaryNutrientSignal = ones(length(currentCurve_correctedTimes),1)*.5;
                isShift = [NaN; diff(binaryNutrientSignal)];
                
            elseif unique(condVals) == 4
                
                binaryNutrientSignal = ones(length(currentCurve_correctedTimes),1);
                isShift = [NaN; diff(binaryNutrientSignal)];
                
            end
            
            % 11. save calculated growth rate and meta data
            track_dv(currentCurve_indeces,1) = dVdt;
            track_dv(currentCurve_indeces,2) = binaryNutrientSignal;
            track_dv(currentCurve_indeces,3) = isShift;
            
        end
        
    end
    
    % 12. concatenate scores for all tracks
    dvdt_data = [dvdt_data; track_dv];
    
end

% 13. output dvdt and nutrient meta data
end