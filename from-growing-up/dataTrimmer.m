%% Extracting the useful subset from original data array D (ND2Proc.m)
%
%
%  SELECTION CRITERIA:
%
%   1) Tracks must be of reasonable size, SizeStrainer (4um)
%           - removes tracks that don't reach desired size  
%
%
%   2) Tracks that do not increase by more than JumpFrac (30% of size at previous timepoint) 
%           - removes tracks affected by settling of upstream cells
%
%           - finds tpt at which a jump in size occurs
%           - trims track at that tpt
%           - replaces original track with trimmed one
%
%
%   3) Tracks that are at least NumTpts (60 mins) long
%           - saturating conc double every ____________ mins
%
%           - finds length of each track
%           - removes tracks that are not of sufficient length
%
%
%   4) Tracks must undergo at least one doubling
%           - permits fitting of exponential to entire cell cycle
%           - permits identification of birth size and division size
%           - gives at least one doubling time
%
%           - finds ratio of max:min length in each track
%           - removes track if ratio is less than the GoldenRatio (1.8)
%           - note: this ideally catches doublings even if the cells are shrinking in size! 
%           
%

% retired: 2018 March 29. last use was over one year ago

%% ONE: throwing out particles that are too small to be growing cells
%
%           - catches particles that sneak past size selection in intitial track detection 
%           - keep threshold low to catch cells that divide but actually shrink in size
%

D2 = D;                                                                    % copies D for modification
SizeStrainer = 1.5;

for n = 1:length(D);                                                       % will give error for "undefined variable LData1" if some cells have no tracks, []    
    
    for i = 1:length(D2{n})
        % for all particles in series n
        LData1{i} = max(D2{n}(i).MajAx);                                   % gets largest MajAx value in track i   
    end

    LData1_double = cell2mat(LData1);                                      % covert to doubles for find function
    ByeTiny = find(LData1_double < SizeStrainer);                          % finds tracks that don't exceed __ um
    
    X = ['Removing ', num2str(length(ByeTiny)), ' small particles from D2(', num2str(n), ')...'];
    disp(X)
    
    toes = 0;
    for s = 1:length(ByeTiny)                                              % loops operates fine even with length(ByeTiny) = 0
        t = length(ByeTiny) - toes;                                        % remove structures based on row #
        D2{n}(ByeTiny(t)) = [];                                            % remove in reverse order
        toes = toes + 1;
        %disp(num2str(toes))
    end

    clear LData1 LData1_double i ByeTiny toes s t X;                       % important for loop
end 

clear SizeStrainer n;

%% TWO: trimming down tracks to remove unreasonable jumps in cell size
%          
%              - if too positive (cells shouldn't double within three minutes)
%              - in these cases, what causes these large jumps? check!
%              - negatives are OK because cells have to divide!

D3 = D2;  
JumpFrac = 0.3;                                                            % JumpFrac = threshold parameter
                                                                           % tracks that increase by a cell size fraction greater than JumpFrac will be eliminated from final dataset
for n = 1:length(D);                                                       
    counter = 0;                                                           
    D3{n} = rmfield(D3{n},'Conv');                                         % remove the 'Conv' field, as it is only one element
    
    for i = 1:length(D3{n})                                                % error "undefined variable" if some cells have no tracks, [], run groups separately as needed        
                                                                           
        % Rates{n,i} = diff(D5{n}(i).MajAx);                               % Rates{i} = instantaneous growth rate for track i
        Rates{i} = diff(D3{n}(i).MajAx);                                   % Rates = row of cells, where each cell is one track from that series
                                                                           
        Clumpers{n,i} = find(Rates{i} > JumpFrac);                         % Clumpers{n,i} creates an array where each row is a series, and each column is a track from that series.
                                                                           % Each cell lists the tpt when Rate jumps greater than fraction indicated (JumpFrac)                                                                                                                         
        tf = find(Clumpers{n,i});
        if tf > 0
            counter = counter + 1;
            %ClumpedTracks(n,counter) = i;                                 % ClumpedTracks(n,counter) = each row is a list of tracks(#) that contain increases greater than JumpFrac threshold
            ClumpedTrack(counter) = i;
            
            for z = 1:length(ClumpedTrack)
                Target = D3{1,n}(ClumpedTrack(z));                             % pulls out structure for target track z, will trim all variables (MajAx, X, Y, etc.)
                Tpt = Clumpers{n,ClumpedTrack(z)};                             % defines timepoint at point of jump
                TrimTarget = structfun(@(M) M(1:Tpt), Target, 'Uniform', 0);   % trims structure to desired timepoint
                D3{1,n}(ClumpedTrack(z)) = TrimTarget;                         % redefines track in data set as trimmed structure
            end
            
        else
            continue
        end
        
    end
    
    X = ['Trimming ', num2str(counter), ' jumps from D3(', num2str(n), ')...'];
    disp(X)
    
    if n < length(D)
        clear Rates ClumpedTrack z;                                        % erases info from current series, so that tracks don't roll into next iteration
    end
    
end
    
clear tf JumpFrac Target counter ClumpedTrack Rates Tpt TrimTarget i n z X;
clear Clumpers;

%% THREE: removing tracks that are insufficient in length (time)
%
%           - track trimming in Section Two can produce very short trajectories 
%           - remove these, as we only want to consider tracks with at least one doubling
%

D4 = D3;
Shortest = 20;                                                             % each timepoint = 3 mins; 30 tpts = 90 mins (1.5 hrs)

for n = 1:length(D);

    for i = 1:length(D4{n})
        % for all particles in series n
        LData2{i} = length(D4{n}(i).MajAx);                                % get number of timepoints in track i 
    end
    
    LData2_d = cell2mat(LData2);                                           % covert to doubles for find function
    ByeShortie = find(LData2_d < Shortest);                                % NumTpts must be set at start of section
    
    X = ['Removing ', num2str(length(ByeShortie)), ' short tracks from D4(', num2str(n), ')...'];
    disp(X)
    
    if isempty(ByeShortie) == 1
        continue
    end
    
    fingers = 0;
    for q = 1:length(ByeShortie)
        r = length(ByeShortie) - fingers;                                  % remove structures based on row #
        D4{n}(ByeShortie(r)) = [];                                         % remove in reverse order!
        fingers = fingers + 1;
        %disp(num2str(fingers));
    end
    clear  LData2 LData2_d i ByeShortie fingers q r X;
end

 clear Shortest n;
 
%% FOUR: removing tracks that do not double
%
%           - removes particles that are otherwise recognized as cells
%           

D5 = D4;
GoldenRatio = 1.3;                                                         % desired ratio between max and min lengths in a given track 

for n = 1:length(D);                                                       
    
    for i = 1:length(D5{n})                                                
        LData1{i} = arrayfun(@(Q) max(Q.MajAx), D5{n}(i));                 % pull out max MajAx value
        LData2{i} = arrayfun(@(Q) min(Q.MajAx), D5{n}(i));                 % pull out min MajAx value
        LData3{i} = LData1{i}/LData2{i};
    end
    
    LData3 = cell2mat(LData3);                                             % convert cells to double
    ToCut = find(LData3 < GoldenRatio);                                    % 'gt' doesn't work on cells
    
    X = ['Removing ', num2str(length(ToCut)), ' non-doublers from D5(', num2str(n), ')...']; 
    disp(X)                                                         
    
    %Check = length(DO) + length(ToCut)
                                                   
    counter = 0;
    for j = 1:length(ToCut)
        k = length(ToCut) - counter;                                       % remove structures based on row #
        D5{n}(ToCut(k)) = [];                                              % remove in reverse order to avoid changing smaller positions
        counter = counter + 1;
    end
    
    clear LData1 LData2 LData3 i j k X counter ToCut;        
    
end
clear GoldenRatio n;


%% Saving results

save('t300_2016-11-23-trimmed.mat', 'D', 'D2', 'D3', 'D4', 'D5', 'T')%, 'reader', 'ConversionFactor')


%% QUALITY CONTROL
%
%  Goal: manually examine length trajectories to ensure successful trimming
%  Approach: display trajectories for manual discard (if needed)
%
%
%  Step 1 - Plotting individual tracks from a series
%         - Choose to either accept or discard each track
%         

figure(1)

for n=31:40                                                                 % adjust n as needed!
    counter = 0;
    Delete = zeros(1,length(n));
    
    J = ['Click mouse to flag track for deletion'];
    disp(J);
    K = ['Hit keyboard to approve track'];
    disp(K);
    
    for m=1:length(D5{n})
        plot(T{n}(D5{n}(m).Frame(1:end))/3600,(D5{n}(m).MajAx),'Linewidth',2)
        axis([0,10,0,15])
        drawnow                                                             % displays current track m
        %pause()
        w = waitforbuttonpress;                                             % waits for user to approve or reject
        % mouse click = flags track for deletion (w = 0)
                                                                           % keyboard key = OK                      (w = 1) 
       if w == 0
           counter = counter + 1;
           Delete(counter) = m;                                            % saves track number for future removal
           X = ['Track ', num2str(m), ' of ', num2str(length(D5{n})), ' marked for deletion...'];
           disp(X);
       else
           disp('OK!');
       end
    end
    
    Rejects{n} = Delete;
    clear Delete m w X J K counter;
    save('2016-11-23-Rejects.mat', 'Rejects')                              % saves current Rejects after finishing each series
end                                                                        % both D5 and Rejects are saved for potential revisitation of removed data

clear n;

%% Quality control, continued...
%
%  Step 2 - Fine pass through Reject list to confirm deletion
%
%

figure(1)
for n=1:40                                                              % adjust n as needed!
    counter = 0;                                                           
    Confirmed = zeros(1,length(n));
    
    if Rejects{n} == 0
        continue
    end
    
    J = ['Click mouse to approve deletion'];
    disp(J);
    K = ['Hit keyboard to rescue track'];
    disp(K);
    
    
    for i=1:length(Rejects{n})                                             % number of tracks in Rejects piles
        m = Rejects{n}(i);

            plot(T{n}(D5{n}(m).Frame(1:end))/3600,(D5{n}(m).MajAx),'Linewidth',2)
            axis([0,11,0,15])
            %plot(T(D5{n}(m).Frame(1:end),n)/60,(D5{n}(m).MajAx),'color',[0,0,1]+(n-1)*[.05,.05,0],'Linewidth',2)
            %axis([0,1100,0,15])
            drawnow                                                             % displays current track m
            %pause()
            w = waitforbuttonpress;                                             % waits for user to approve or reject
            % mouse click = flags track for deletion (w = 0)
            % keyboard key = OK                      (w = 1)
            if w == 0
                counter = counter + 1;
                Confirmed(counter) = m;                                            % saves track number for future removal
                X = [ num2str(i), ' track of ', num2str(length(Rejects{n})), ' initial rejects confirmed for deletion...'];
                disp(X);
            else
                disp('Rescued!');
            end
            
        end
        Trash{n} = Confirmed;

    clear Confirmed m w X J K counter;
    save('2016-11-23-Rejects.mat', 'Rejects','Trash')
                                                                        % saves current Rejects after finishing each series
end                                                                        % both D5 and Rejects are saved for potential revisitation of removed data

clear n;

%% Quality control, continued... 
%
%  Step 3 - Removing flagged tracks from data set
%         - Delete based on track numbers saved in 'Rejects'
%
D6 = D5;
Trash = cellfun(@(x)x(logical(x)),Trash,'uni',false); 

for n = 1:length(D6);                                                       
    
    Remove = Trash(n);
    Remove = cell2mat(Remove);
               % replaces 0 with empty cells
    F = ['Removing ', num2str(length(Remove)), ' flagged tracks from D6(', num2str(n), ')...']; 
    disp(F)                                                         
                                                   
    counter = 0;
    for j = 1:length(Remove)
        k = length(Remove) - counter;                                      % remove structures based on row #
        D6{n}(Remove(k)) = [];                                             % remove in reverse order to avoid changing smaller positions
        counter = counter + 1;
    end
    clear j k F counter Remove;           
end
%clear n;
save('2016-11-21-trimmed.mat', 'D', 'D2', 'D3', 'D4', 'D5', 'D6', 'T')%, 'reader', 'ConversionFactor')