function [] = Analyze_4c_Dvel_100()
data_type = '4c_Dvel_100';

cd ..
code_path = pwd;
addpath(genpath(code_path));

cd ../..
root = pwd;
InputPath = strcat(root,'/Results');
dataRootPath = strcat(root,'/data');

dataOutPath = InputPath;

Analyze_masses(InputPath, data_type, dataRootPath, dataOutPath)

rmpath(genpath(code_path))

end