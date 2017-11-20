function [] = Analyze_masses_new(InputPath, data_type,dataRootPath1,dataOutPath1)


%load 'variables.mat'
%load 'variables2.mat' % threshold = 5
%load 'variables3.mat' % threshold = 25
%load 'variables4.mat' % threshold = 15 and L0 = 0.5b
%load 'variables5.mat' % threshold = 5 and L0 = 0.5b
%load 'variables6.mat' % threshold = 5 and L0 = 0.5b and Lmax = 1.25a
%load 'variables7.mat' % threshold = 3 and L0 = 0 and Lmax = 1.25a
%load 'variables8.mat' % threshold = 5 and L0 = 0.5b and Lmax = 1.25a compute AREA
%load 'variables9.mat' % threshold = 5 and L0 = 0.5b and Lmax = 1.25a compute AREA, 70% of channels
%load 'variables10.mat' % threshold = 5 and L0 = 0.5b and Lmax = 1.25a compute AREA, 70% of channels (CORRECTED INPUT)

%dataRootPath = '/Users/CesarMB13/Google Drive/My Journal papers/In preparation/Shell_mass/Data/fits/';
%dataRootPath = '/N/dc2/projects/lifebid/code/ccaiafa/Shells/data'; %Karst path
%dataOutPath = '/Users/CesarMB13/Google Drive/My Journal papers/In preparation/Shell_mass/Results/4c_Dvel_100/';


Results_file = fullfile(InputPath,strcat('Results_',data_type,'.mat'));
load(Results_file);

Area_vs_alpha = squeeze(Area(:,1,:));

dataRootPath = dataRootPath1;
dataOutPath = dataOutPath1;

shell_candidates = shell_all;
alpha = alpha_range;
delta = delta_range;

optimal_alpha = zeros(N,1);
ind_alpha = zeros(N,1);
error_alpha = zeros(N,1);

%% Search for the best alpha based on ellipse area
for n=1:N
    area_diff = abs((Area_vs_alpha(:,n) - pi*shell_candidates{n}.a*shell_candidates{n}.b)/(pi*shell_candidates{n}.a*shell_candidates{n}.b));
    [error_alpha(n),ind_alpha(n)] = min(area_diff);
    optimal_alpha(n) = alpha(ind_alpha(n));
end

%% Search for the best delta based on difference between shell mass and missing mass and store final estimated Mass and missing
optimal_delta = zeros(N,1);
ind_delta = zeros(N,1);
error_delta = zeros(N,1);
Mass_estimated = zeros(N,1);
Missing_estimated = zeros(N,1);
Mass_measured = zeros(N,1);
Missing_measured = zeros(N,1);
for n=1:N
    Mass_vs_delta = squeeze(Mass(ind_alpha(n),:,n));
    Missing_vs_delta = squeeze(Missing_Mass(ind_alpha(n),:,n));
    mass_diff = abs((Mass_vs_delta - Missing_vs_delta)./Mass_vs_delta);
    [error_delta(n),ind_delta(n)] = min(mass_diff);
    optimal_delta(n) = delta(ind_delta(n));
    Mass_estimated(n) = Mass(ind_alpha(n), ind_delta(n), n);
    Missing_estimated(n) = Missing_Mass(ind_alpha(n), ind_delta(n), n);    
    Mass_measured(n) = shell_candidates{n}.MassShell;
    Missing_measured(n) = shell_candidates{n}.MassMiss;
end

%% Calcular R2
Min_mean = mean(Missing_measured,1);
SSres_Min = sum((Missing_measured - Missing_estimated).^2,1);
SStot_Min = sum((Missing_estimated - Min_mean).^2,1);
R2_Min = 1 - SSres_Min/SStot_Min;

Mass_mean = mean(Mass_measured,1);
SSres_Mass = sum((Mass_measured - Mass_estimated).^2,1);
SStot_Mass = sum((Mass_measured - Mass_mean).^2,1);
R2_Mass = 1 - SSres_Mass/SStot_Mass;

SSres_datos = sum((Missing_measured - Mass_measured).^2,1);
SStot_datos = sum((Mass_measured - Min_mean).^2,1);
R2_datos = 1 - SSres_datos/SStot_datos;

rmse_Mass = sqrt(nanmean((Mass_measured - Mass_estimated).^2));
rmse_Miss = sqrt(nanmean((Missing_measured - Missing_estimated).^2));

disp(['Coeficiente de determinacion R2:', 'Missing=',num2str(R2_Min),', Mass=',num2str(R2_Mass),', Datos=',num2str(R2_datos)]);

%% Scatter plot de MIssing Mass
figure
hold on
scatter(Missing_measured, Missing_estimated,'DisplayName','Missing Mass') 
scatter(Mass_measured, Mass_estimated,'DisplayName','Shell Mass')
%scatter([masses(:,1);masses(:,2)],[masses(:,3);masses(:,4)],'DisplayName','Mass')
ylabel({'Automatically estimated'});
% Create xlabel
xlabel({'Estimated by hand'});
ylim([0 3000000]);
xlim([0 3000000]);
title(['rmse Missing=',num2str(rmse_Miss),', rmse Mass=', num2str(rmse_Mass)])

rango_masas = 0:1000:3000000;

plot(rango_masas, rango_masas,'k','DisplayName','Perfect estimation')
set(gca,'DataAspectRatio',[1499500 1499500 1],'XMinorTick','on','XScale','log','YMinorTick','on','YScale','log')
legend(gca,'show');

%% Print Masses
filename = fullfile(dataOutPath,strcat('Masses_',data_type,'.txt'));
fileID = fopen(filename,'w');
fprintf(fileID,'Name   \t Mass (hand)  \t Miss (hand) \t Mass (alg)  \t Miss (alg) \t dV\n');

for n=1:N
    fprintf(fileID,'%s   \t %12.1f \t %12.1f \t %12.1f \t %12.1f \t %12.1f \n', shell_all{n}.name, shell_candidates{n}.MassShell, shell_candidates{n}.MassMiss, Mass_estimated(n), Missing_estimated(n), shell_candidates{n}.dV); 
end
fclose(fileID);

%% Visualize structures
Visualization = 1;
perc = 0.7; % Velocity percentages
dir_name = fullfile(dataOutPath,strcat('Figs_',data_type));
mkdir(dir_name);
figure
%for n=1:N
for n = 1:N % 3cdataset
%for n=[31,53,37,59,11,20,56,46,2,49,62,25,9,29,44,47] %4c dataset
    disp(['Displaying shell ',num2str(n)])
    [~, cube, header] = select_cube(dataRootPath,shell_candidates{n});
    [ mass, missing_mass, diameter ] = compute_mass_V5( cube, header,  optimal_alpha(n), optimal_delta(n), shell_candidates{n}, Visualization, perc);
    file_fig_name = fullfile(dataOutPath,strcat('Figs_',data_type),strcat('Fig_shell_',num2str(n),'.pdf'));
    print(file_fig_name,'-dpdf')
    clf('reset')
end

end