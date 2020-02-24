% showLengthTracks.m

% display length trajectories, just for funsies
%
%

load('2015-08-10-Mus-length.mat');

n = 1;
for m=1:length(D6{n})
    plot(D6{n}(m).MajAx,'Linewidth',2);
    hold on
    %axis([0,1000,0,15])
end
    