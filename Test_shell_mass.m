% 

dataRootPath = '/Users/CesarMB13/Google Drive/My Journal papers/In preparation/Shell_mass/Data/fits/';
dataOutPath = '/Users/CesarMB13/Google Drive/My Journal papers/In preparation/Shell_mass/Results/4c_Dvel_100/';

perc_train = 0.7;
perc_out = 0.3;
%alpha_range = 0.1:0.01:0.5;
alpha_range = [0.4,0.5];
%delta_range = 0:0.005:0.1;
delta_range = [0.4,0.5];

Ncross = 10;

Visualization = 0;

[ind_train, ind_test] = Crossvalidation(dataRootPath, dataOutPath, perc_train, perc_out, alpha_range, delta_range, Ncross, Visualization);

