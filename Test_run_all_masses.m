% 
clear
clc

%dataRootPath = '/Users/CesarMB13/Google Drive/My Journal papers/In preparation/Shell_mass/Data/fits/';
dataRootPath = '/N/dc2/projects/lifebid/code/ccaiafa/Shells/data'; %Karst path
%dataOutPath = '/Users/CesarMB13/Google Drive/My Journal papers/In preparation/Shell_mass/Results/4c_Dvel_100/';
dataOutPath = '/N/dc2/projects/lifebid/code/ccaiafa/Shells/Results'; %Karst path

alpha_range = 0.1:0.01:0.5;
%alpha_range = [0.15,0.18];
delta_range = 0.:0.001:0.1;
%delta_range = [0.15, 0.25];

Visualization = 0; % Visualization ONLY works in Cubo 1 (long range 80 - 130)

run_all_masses(dataRootPath, @load_shells_4c_Dvel_100, dataOutPath, alpha_range, delta_range, Visualization);
