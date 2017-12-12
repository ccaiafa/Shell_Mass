function [] = Gen_fig_shell_4c_Dvel_100()
data_type = '4c_Dvel_100';

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