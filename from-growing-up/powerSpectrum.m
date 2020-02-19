%% Plot power spectrum from equilibrated regions of growth signals
%
%  This script:
%
%       1. plots section of fluctuating growth signal normalized to stable
%       2. uses that growth signal to generate a power spectrum
%  
%  We expect that any growth fluctuations resulting from fluctuating
%  nutrient signals should be distinct in frequency from noisy growth
%  fluctuations in stable environments
%

%%  SECTION ONE -- part A

%     Initializing growth rate signal     %
%

load('2016-07-30-Mus-length.mat');

% defining conditions: col1 = first xy; col2 = final xy; col3 = time (hr) cutoff
conditions = [1 10 6.5; 11 20 6.5; 21 30 4; 31 40 3.5];

% generate matrix of stats for each condition
% each row = a condition
% each column = a metric, or 'field'
muStats = [];


for i = 1:4 %number of conditions
    
    Mu_cond = [];
    Time_cond = [];
    
    for n = conditions(i,1):conditions(i,2)
        for m = 1:length(M6{n})
            
            %  assemble all instantaneous growth rates into a single vector
            Mu_cond = [Mu_cond; M6{n}(m).Parameters(:,1)];
            
            %  assemble a corresponding timestamp vector
            vectorLength = length(M6{n}(m).Parameters(:,1));
            trackFrames = D6{n}(m).Frame(3:vectorLength+2);
            Time_cond = [Time_cond; T{n}(trackFrames)];
            
        end
    end
    
    %  convert all timestamps from seconds to hours
    Time_cond = Time_cond/3600;
    
    %  eliminate negative growth rates
    %Mu_cond1(Mu_cond1<0)=NaN;
    

    BinsPerHour = 60;                 % multiplying by 10 gives bins of 0.1 hr
    Bins = ceil(Time_cond*BinsPerHour);            
    
    %  accumulate growth rates by bin, and calculate mean and std dev
    muMeans = accumarray(Bins,Mu_cond,[],@nanmean);
    muDevs = accumarray(Bins,Mu_cond,[],@nanstd);
    
    
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
        muCounts(j) = counter;
        clear i counter Kasten;
    end
    
    %   2. divide standard dev by square root of tracks per bin
    muSEMs = muDevs./sqrt(muCounts');
    
    s = struct('muMeans', muMeans,  'muDevs', muDevs, 'muCounts', muCounts, 'muSEMs', muSEMs);
    muStats = [muStats s];
    
    clear vectorLength trackFrames muMeans muDevs muSEMs Bins hr dT muCounts n m j;
    clear Mu_cond Time_cond;
    
end

clear BinsPerHour conditions currentBin_count s;


%%  SECTION ONE -- part B


%     Isolate and plot normalized region of interest within growth signal    %
%

frameRate = 1;                  % sampling frequency (1/mins)
L = 250;                        % number of points in growth rate data signal
time = T{:,1}/(60*60);          % time vector (hr), uses timestamp from 1st position

flucSignal = muStats(1).muMeans;               % growth rate signal for fluctuating environment
constSignal = muStats(3).muMeans;               %   "                 "  constant environment

%figure(1)
%plot(time(2:end-1),flucSignal)
%hold on
%plot(time(2:end-1),constSignal,'k')
%grid on



%%     Fast Fourier Transform     %
%       - method 1

fsTrimmed = flucSignal(150:250);
fsTrimmed = fsTrimmed-mean(fsTrimmed);
csTrimmed = constSignal(150:250);
csTrimmed = csTrimmed-mean(csTrimmed);
NFFT = 2^nextpow2(L);
FFT_f = fft(fsTrimmed,NFFT)/L;         % discrete Fourier transform of Mf
FFT_c = fft(csTrimmed,NFFT)/L;         % discrete Fourier transform of Mc

f = frameRate/2*linspace(0,1,NFFT/2+1);        % frequency vector

%figure(2)
plot(f,2*abs(FFT_f(1:NFFT/2+1)))        % single-sided amplitutde spectrum for fluctuating data
hold on
plot(f,2*abs(FFT_c(1:NFFT/2+1)),'k')    % single-sided amplitutde spectrum for constant data
grid on
axis([0,.2,0,.05])
xlabel('Frequency (1/min)')
ylabel('|FFT|(f)') 
hold on

%%
%     Plotting window of oscillation     %
%       

% raw data, post adaptation
figure(1)
plot(1:185, constSignal(1:185),'b')
hold on
plot(190:250, constSignal(190:250),'b')
hold on
plot(1:250,flucSignal(1:250),'k')
grid on
axis([1,250,-0.05,.53])
xlabel('Time (hr)')
ylabel('Average growth rate (1/hr)')

figure(2)
errorbar(1:185, muStats(3).muMeans(1:185), muStats(3).muSEMs(1:185),'b')
hold on
errorbar(190:250, muStats(3).muMeans(190:250), muStats(3).muSEMs(190:250),'b')
hold on
errorbar(1:250, muStats(1).muMeans(1:250),muStats(1).muSEMs(1:250),'k')
grid on
axis([1,250,-0.05,.53])
xlabel('Time (hr)')
ylabel('Average growth rate (1/hr)')

%%
% fluctuating signal, relative to mean rate in constant
normalizedFS_3 = flucSignal(180:390)/mean(constSignal(130:170)); 
normallizedFS_4 = flucSignal(180:390)-mean(constSignal(130:170));

figure(3)
plot(time(180:390),normalizedFS_3)
grid on
hold on
axis([3,6.5,0.6,1.4])
xlabel('Time (hr)')
ylabel('Relative growth rate (1/hr)') 


figure(4)
plot(time(180:390),normallizedFS_4)
grid on
hold on
axis([3,6.5,-.2,.2])
xlabel('Time (hr)')
ylabel('Relative growth rate (1/hr)') 


% driving signal (nutrient) and growth oscillations

%t2 = time(200:399);
%Ds = square(t2);
%Ds2 = .5*square(2*pi*15*t2);
%plot(t2,Ds,t2,Ds2)
%axis([5,7,-1.5,1.5])


