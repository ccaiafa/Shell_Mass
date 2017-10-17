function [] = Analyze_3c_total()
restoredefaultpath

data_type = '3c_total';

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