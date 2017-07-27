
function [ind_train, ind_test] = Crossvalidation(dataRootPath,dataOutPath, perc_train, perc_out, alpha_range, delta_range, Ncross, Visualization)

% Load full dataset
shell_all = load_shells_4c_Dvel_100();
N = size(shell_all,2);
N_train = round(perc_train*N);
N_test = N - N_train;

% Create and save Training and Testing datasets
ind_rnd = randperm(N);
ind_train = ind_rnd(1:N_train);
ind_test = ind_rnd(N_train+1:end);

shell_train = shell_all(ind_train);
shell_test = shell_all(ind_test);

save(fullfile(dataOutPath,sprintf('Train_Test_shells.mat')), 'shell_all','shell_train','shell_test','-v7.3')

% Tuning of hyperparameters alpha and delta by take p-out crossvalidation
[alpha_opt, delta_opt, mean_error_train, std_error_train, mean_error_test, std_error_test] = hyper_tuning(dataRootPath, shell_train, alpha_range, delta_range, perc_out, Ncross, Visualization);

end



