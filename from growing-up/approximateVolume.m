%% approximateVolume.m

%  Goal: approximate a cell as a cylinder and gather volume from our
%  measurements of length and width

%  Last update: jen, 2017 March 31

%% volume of a cylinder

%       V = pi * length * (width/2)^2

% 1. for each timepoint in each cell track,
%    caluculate V_cylinder using:

%        length = D6{n}(m).MajAx;
%        width  = D6{n}(m).MinAx;

% 2. append volume vector to D7{n}

%%

load('t300_2017-01-18-autoTrimmed.mat');
clear D D2 D3 D4 D5;

%%

for n = 1:length(D6)
    
    Vc = [];
    Ve = [];
    
    % initialize volume vectors
    for m = 1:length(D6{n})
        
        % gather length and width vs. time vectors
        length = D6{n}(m).MajAx;
        width = D6{n}(m).MinAx;
        
        % calculate instantaneous volume of cylinder
        v_cylinder = pi * length .* (width/2).^2;
        
        % calculate instantaneous colume of ellipse
        v_ellipse = 4/3 * pi * length/2 .* (width/2.^2);
     
        
        Vc{m,1} = v_cylinder;
        Ve{m,1} = v_ellipse;
        
    end
        
end

D7=D6;
