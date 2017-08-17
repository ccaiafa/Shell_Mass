clear


%load 'variables.mat'
%load 'variables2.mat' % threshold = 5
%load 'variables3.mat' % threshold = 25
%load 'variables4.mat' % threshold = 15 and L0 = 0.5b
%load 'variables5.mat' % threshold = 5 and L0 = 0.5b
%load 'variables6.mat' % threshold = 5 and L0 = 0.5b and Lmax = 1.25a
%load 'variables7.mat' % threshold = 3 and L0 = 0 and Lmax = 1.25a
%load 'variables8.mat' % threshold = 5 and L0 = 0.5b and Lmax = 1.25a compute AREA
load 'variables9.mat' % threshold = 5 and L0 = 0.5b and Lmax = 1.25a compute AREA, 70% of channels

dataRootPath = '/Users/CesarMB13/Google Drive/My Journal papers/In preparation/Shell_mass/Data/fits/';
%dataRootPath = '/N/dc2/projects/lifebid/code/ccaiafa/Shells/data'; %Karst path
dataOutPath = '/Users/CesarMB13/Google Drive/My Journal papers/In preparation/Shell_mass/Results/4c_Dvel_100/';

shell_candidates = shell_all;
alpha = alpha_range;
delta = delta_range;

Top_Diam_Diff = 0.25; % Optimal rmse_miss=60,750,  rmse_mass=62,833
perc = 0.7;

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
    Error_auto_Missing(n) = abs(Missing_Mass(local_auto{n}.ia,local_auto{n}.id,n) - shell_candidates{n}.MassMin)/shell_candidates{n}.MassMin;
    Temp_img(n) = Temperatures{n}.img(ia,id);
    Temp_shell(n) = Temperatures{n}.shell(ia,id);
    Temp_backg(n) = Temperatures{n}.backg(ia,id);
    Temp_miss(n) = Temperatures{n}.miss(ia,id);
end

filename = fullfile(dataOutPath,'Temperatures.txt');
fileID = fopen(filename,'w');
fprintf(fileID,'Name   \t Img   \t Shell   \t Backg   \t Miss \n',A);
masses = [];
for n=1:N
    %masses = [masses; shell_candidates{n}.MassMin, shell_candidates{n}.MassMax, local_auto{n}.Missing, local_auto{n}.Mass];
    masses = [masses; shell_candidates{n}.MassMin, shell_candidates{n}.MassMax, local_auto{n}.Missing, local_auto{n}.Mass];
    x{n} = [num2str(n),'=',shell_candidates{n}.name];
    fprintf(fileID,'%s   \t %12.8f \t %12.8f \t %12.8f \t %12.8f\n', shell_all{n}.name, Temp_img(n), Temp_shell(n), Temp_backg(n), Temp_miss(n));
    
end
fclose(fileID);
figure
bar1 = bar(masses);
set(gca, 'XTick', 1:N, 'XTickLabel', x);
ax = gca; 
ax.XTickLabelRotation = 45;
set(bar1(1),'DisplayName','By Hand (MassMin)');
set(bar1(2),'DisplayName','By Hand (MassMax)');
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
Visualization = 1;
figure
for n=59
%for n=1:N
    [~, cube, header] = select_cube(dataRootPath,shell_candidates{n});
    [ mass, missing_mass, diameter ] = compute_mass_V5( cube, header,  local_auto{n}.alpha, local_auto{n}.delta, shell_candidates{n}, Visualization, perc);
    clf('reset')
end

stop






% Figura de Errores
error = 0.7;
figure
axes('YScale','log');
errorbar(masses(:,1),error*masses(:,1),'r','Marker','o','LineStyle','none')
hold on
errorbar(masses(:,3),error*masses(:,3),'b','Marker','o','LineStyle','none')



figure
bar2 = bar([Error_auto_Mass , Error_auto_Missing]);
set(gca, 'XTick', 1:N, 'XTickLabel', x);
ax = gca; 
ax.XTickLabelRotation = 45;
set(bar2(1),'DisplayName','Relative Error Shell Mass');
set(bar2(2),'DisplayName','Relative Error Missing Mass');
% Create legend
legend('show');




dx=0.005;
dy=0.005;
figure
%scatter(global_opt.alpha, global_opt.delta,'o')
%h=text(global_opt.alpha+dx, global_opt.delta+dy, 'Global optimum');
%set(h,'rotation', 45)
hold on
for n=1:N
    s = local_opt{n}.Mass/50;
    %scatter(local_opt{n}.alpha, local_opt{n}.delta,s,'o')
    scatter(local_opt{n}.alpha,local_opt{n}.delta,'o')
    h=text(local_opt{n}.alpha+dx, local_opt{n}.delta+dy, shell_candidates{n}.name);
    set(h,'rotation', 45)
    
%     s = local_auto{n}.Mass/50;
%     scatter(local_auto{n}.alpha, local_auto{n}.delta,s,'s')
%     h=text(local_auto{n}.alpha+dx, local_auto{n}.delta+dy, shell_candidates{n}.name);
%     set(h,'rotation', 45)
end
xlim([min(alpha) max(alpha)]);
ylim([min(delta) max(delta)]);
xlabel('\alpha')
ylabel('\delta')
%xlim([0.1 0.33]);
%ylim([-0.02 0.15]);


 


for n=1:N
%     masa_real_mean = (shell_candidates{n}.MassMin + shell_candidates{n}.MassMax)/2;
%     masa_real_min = shell_candidates{n}.MassMin;
%     masa_real_max = shell_candidates{n}.MassMax;
    
    masa_real = shell_candidates{n}.MassMax;
    
    [val,ind] = max(Mass(:,:,n));
    [masa_est_max, ind2] = max(val);
    
    alpha_max = ind(ind2);
    delta_max = ind2;
    
    [val,ind] = min(Mass(:,:,n));
    [masa_est_min, ind2] = min(val);
    
    alpha_min = ind(ind2);
    delta_min = ind2;
    
    disp([shell_candidates{n}.name,' Masa Real=',num2str(masa_real)]);
    disp([' Masa Estimada: Max=',num2str(masa_est_max),'  :Min=',num2str(masa_est_min)])
    disp(['(\alpha,\delta) for Max =', '(',num2str(alpha_max),',',num2str(delta_max),')'])
    disp(['(\alpha,\delta) for Min =', '(',num2str(alpha_min),',',num2str(delta_min),')'])
    disp('  ')
    
%     figure 
%     imagesc(abs(Mass(:,:,n) - masa_real))
%     title(['Min rel error=',num2str(100*min(min(abs(Mass(:,:,n) - masa_real)))/masa_real),'%'])
%     pause(0.5)
end

save 'results.mat'