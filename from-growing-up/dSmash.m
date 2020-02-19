% smashing D's together

load('letstry-2017-06-12.mat');
D_smash = D;

%%
load('letstry-2017-06-12-largeDudes.mat');

for i = 31:60
    D_smash{i} = D{i};
end
%%

clear T
load('letstry-2017-06-12.mat');

save('letstry-2017-06-12-dSmash.mat','D_smash','T')