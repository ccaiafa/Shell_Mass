% 
clear
clc
restoredefaultpath

clear
clc

cd ..
code_path = pwd;
addpath(genpath(code_path));

cd ../..
root = pwd;
dataOutPath = strcat(root,'/code/Shell_Mass/Results');
dataRootPath = strcat(root,'/data');

OutputFile = 'Results_3c_total';

dR = 0.50;

alpha_range = 0.2:0.01:1.0;
%alpha_range = [0.15,0.18];
delta_range = 0:0.01:0.75;
%delta_range = [0.15, 0.25];

Visualization = 0; % Visualization ONLY works in Cubo 1 (long range 80 - 130)
perc = 0.7; % 70% of channels

run_all_masses(dataRootPath, @load_shells_3c_total, dataOutPath, OutputFile, alpha_range, delta_range, Visualization, perc, dR);

rmpath(genpath(code_path));