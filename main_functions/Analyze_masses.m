function [] = Analyze_masses(InputPath, data_type,dataRootPath1,dataOutPath1)


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

dataRootPath = dataRootPath1;
dataOutPath = dataOutPath1;

shell_candidates = shell_all;
alpha = alpha_range;
delta = delta_range;

Top_Diam_Diff = 0.25; % Optimal rmse_miss=60,750,  rmse_mass=62,833
perc = 0.7; % Velocity percentages

Diff2 = zeros(size(Mass));
Diff_missing = zeros(size(Mass));

[Na,Nd,N] = size(Mass);
Error_local = zeros(N,1);
Error_auto_Mass = zeros(N,1);
Error_auto_Missing = zeros(N,1);
Error_global = zeros(N,1);
Global_est_mass = zeros(N,1);

Temp_img = zeros(N,1);
Temp_shell = zeros(N,1);
Temp_backg = zeros(N,1);
Temp_miss = zeros(N,1);

for n=1:N
    A = abs(Mass(:,:,n) - Missing_Mass(:,:,n)); % Mass exceeds x% of Missing Mass
    %A(A<0) = Inf;
    A(isnan(A)) = Inf;
    
    % Restrict tensor A to the cases where the estimated Diameter is close to the real
    % one
    %% Keep, for example, 10% Top most similar Diameter region.
%     B = abs((Diameter(:,:,n) - 2*shell_candidates{n}.a)/(2*shell_candidates{n}.a));
      B = abs((Area(:,:,n) - pi*shell_candidates{n}.a*shell_candidates{n}.b)/(pi*shell_candidates{n}.a*shell_candidates{n}.b));
     [~, index] = sort(B(:),'ascend');
     
     index = index(round(Top_Diam_Diff*length(index)):end);
     A(ind2sub(size(A),index)) = Inf;

%     A1 = abs((Area(:,:,n) - pi*shell_candidates{n}.a*shell_candidates{n}.b)/(pi*shell_candidates{n}.a*shell_candidates{n}.b));
%     A2 = abs((Mass(:,:,n) - Missing_Mass(:,:,n))./Mass(:,:,n)); % Mass exceeds x% of Missing Mass  
    
    % find minimum difference between mass and missing mass
    %Diff_missing(:,:,n) = A1;
    Diff_missing(:,:,n) = A;
    [val,ind] = min(reshape(Diff_missing(:,:,n),[Na*Nd,1]));
    [ia,id] = ind2sub([Na,Nd],ind);
    
    local_auto{n}.alpha = alpha(ia);
    local_auto{n}.delta = delta(id);
    local_auto{n}.ia = ia;
    local_auto{n}.id = id;
    local_auto{n}.Mass = Mass(ia,id,n);
    local_auto{n}.Missing = Missing_Mass(ia,id,n);  
    Error_auto_Missing(n) = abs(Missing_Mass(local_auto{n}.ia,local_auto{n}.id,n) - shell_candidates{n}.MassMiss)/shell_candidates{n}.MassMiss;
    Temp_img(n) = Temperatures{n}.img(ia,id);
    Temp_shell(n) = Temperatures{n}.shell(ia,id);
    Temp_backg(n) = Temperatures{n}.backg(ia,id);
    Temp_miss(n) = Temperatures{n}.miss(ia,id);
end

%% Print Temperatures
filename = fullfile(dataOutPath,strcat('Temperatures_',data_type,'.txt'));
fileID = fopen(filename,'w');
fprintf(fileID,'Name   \t Shell  \n',A);
masses = [];
for n=1:N
    %masses = [masses; shell_candidates{n}.MassMiss, shell_candidates{n}.MassShell, local_auto{n}.Missing, local_auto{n}.Mass];
    masses = [masses; shell_candidates{n}.MassMiss, shell_candidates{n}.MassShell, local_auto{n}.Missing, local_auto{n}.Mass];
    x{n} = [num2str(n),'=',shell_candidates{n}.name];
    fprintf(fileID,'%s   \t %12.8f \n', shell_all{n}.name, Temp_shell(n));
    
end
fclose(fileID);

%% Print Masses
filename = fullfile(dataOutPath,strcat('Masses_',data_type,'.txt'));
fileID = fopen(filename,'w');
fprintf(fileID,'Name   \t Mass (hand)  \t Miss (hand) \t Mass (alg)  \t Miss (alg) \t dV\n',A);

for n=1:N
    fprintf(fileID,'%s   \t %12.1f \t %12.1f \t %12.1f \t %12.1f \t %12.1f \n', shell_all{n}.name, shell_candidates{n}.MassShell, shell_candidates{n}.MassMiss, local_auto{n}.Mass, local_auto{n}.Missing, shell_candidates{n}.dV); 
end
fclose(fileID);

%% Calcular R2
%masses = log(masses);
Min_mean = mean(masses(:,1),1);
SSres_Min = sum((masses(:,1)-masses(:,3)).^2,1);
SStot_Min = sum((masses(:,3)-Min_mean).^2,1);
R2_Min = 1 - SSres_Min/SStot_Min;

Mass_mean = mean(masses(:,2),1);
SSres_Mass = sum((masses(:,2)-masses(:,4)).^2,1);
SStot_Mass = sum((masses(:,2)-Mass_mean).^2,1);
R2_Mass = 1 - SSres_Mass/SStot_Mass;

SSres_datos = sum((masses(:,1)-masses(:,2)).^2,1);
SStot_datos = sum((masses(:,2)-Min_mean).^2,1);
R2_datos = 1 - SSres_datos/SStot_datos;

disp(['Coeficiente de determinacion R2:', 'Missing=',num2str(R2_Min),', Mass=',num2str(R2_Mass),', Datos=',num2str(R2_datos)]);


figure
bar1 = bar(masses);
set(gca, 'XTick', 1:N, 'XTickLabel', x);
ax = gca; 
ax.XTickLabelRotation = 45;
set(bar1(1),'DisplayName','By Hand (MassMiss)');
set(bar1(2),'DisplayName','By Hand (MassShell)');
set(bar1(3),'DisplayName','Algorithm (MissingMass)');
set(bar1(4),'DisplayName','Algorithm (Mass)');

rmse_Mass = sqrt(nanmean((masses(:,2)-masses(:,4)).^2));
rmse_Miss = sqrt(nanmean((masses(:,1)-masses(:,3)).^2));

% Create legend
legend('show');
title(['Top ',num2str(100*Top_Diam_Diff),'%,  rmse Missing=',num2str(rmse_Miss),', rmse Mass=', num2str(rmse_Mass)])

%% Scatter plot de MIssing Mass
figure
hold on
scatter(masses(:,1),masses(:,3),'DisplayName','Missing Mass') 
scatter(masses(:,2),masses(:,4),'DisplayName','Shell Mass')
%scatter([masses(:,1);masses(:,2)],[masses(:,3);masses(:,4)],'DisplayName','Mass')
ylabel({'Automatically estimated'});
% Create xlabel
xlabel({'Estimated by hand'});
ylim([0 3000000]);
xlim([0 3000000]);

rango_masas = 0:1000:3000000;

plot(rango_masas, rango_masas,'k','DisplayName','Perfect estimation')
set(gca,'DataAspectRatio',[1499500 1499500 1],'XMinorTick','on','XScale','log','YMinorTick','on','YScale','log')
legend(gca,'show');
%plot(rango_masas, 0.6*rango_masas,'g')
%plot(rango_masas, 1.6*rango_masas,'g')

%% Compute Errors
diff_mass = 100*(masses(:,2) - masses(:,4))./masses(:,2);
diff_missing = 100*(masses(:,1) - masses(:,3))./masses(:,1);
figure
bar2 = bar([diff_mass, diff_missing]);
set(gca, 'XTick', 1:N, 'XTickLabel', x);
ax = gca; 
ax.XTickLabelRotation = 45;
set(bar2(1),'DisplayName','Error Mass');
set(bar2(2),'DisplayName','Error Missing');
legend('show');

%% Compute Differences
hand_diff = 100*(masses(:,2) - masses(:,1))./masses(:,2);
alg_diff = 100*(masses(:,4) - masses(:,3))./masses(:,3);
figure
bar3 = bar([hand_diff, alg_diff]);
set(gca, 'XTick', 1:N, 'XTickLabel', x);
ax = gca; 
ax.XTickLabelRotation = 45;
set(bar3(1),'DisplayName','Diff by Hand');
set(bar3(2),'DisplayName','Diff by Algorithm');
legend('show');


%% Visualize structures
%Visualization = 1;
%figure
%for n=15:N
%for n=1:N
%    [~, cube, header] = select_cube(dataRootPath,shell_candidates{n});
%    [ mass, missing_mass, diameter ] = compute_mass_V5( cube, header,  local_auto{n}.alpha, local_auto{n}.delta, shell_candidates{n}, Visualization, perc);
%    clf('reset')
%end
end