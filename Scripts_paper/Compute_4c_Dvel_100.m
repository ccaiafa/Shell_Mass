% 
clear
clc

cd ..
code_path = pwd;
addpath(genpath(code_path));

cd ../..
root = pwd;
dataOutPath = strcat(root,'/code/Shell_Mass/Results');
dataRootPath = strcat(root,'/data');

OutputFile = 'Results_4c_Dvel_100';

alpha_range = 0.1:0.01:0.5;
%alpha_range = [0.15,0.18];
delta_range = 0.:0.001:0.1;
%delta_range = [0.15, 0.25];

Visualization = 0; % Visualization ONLY works in Cubo 1 (long range 80 - 130)
perc = 0.7; % 70% of channels

run_all_masses(dataRootPath, @load_shells_4c_Dvel_100, dataOutPath, OutputFile, alpha_range, delta_range, Visualization, perc);

rmpath(genpath(code_path));