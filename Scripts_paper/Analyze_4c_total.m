function [] = Analyze_4c_total()
restoredefaultpath
cd('/Users/CesarMB13/Google Drive/My Journal papers/In preparation/Shell_mass/code/Shell_Mass/Scripts_paper')
close all

data_type = '4c_total';

cd ..
code_path = pwd;
addpath(genpath(code_path));

cd ../..
root = pwd;
InputPath = strcat(root,'/code/Shell_Mass/Results');
dataRootPath = strcat(root,'/data');

dataOutPath = InputPath;

Analyze_masses(InputPath, data_type, dataRootPath, dataOutPath)

rmpath(genpath(code_path))

end