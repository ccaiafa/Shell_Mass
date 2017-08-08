clear

%load 'variables.mat'
%load 'variables2.mat' % threshold = 5
%load 'variables3.mat' % threshold = 25
load 'variables4.mat' % threshold = 15

shell_candidates = load_shells_4c_Dvel_100();
alpha = 0.1:0.01:0.5;
delta = 0.:0.001:0.1;

Top_Diam_Diff = 0.10; % Optimal rmse_miss=60,750,  rmse_mass=62,833

Diff2 = zeros(size(Mass));
Diff_missing = zeros(size(Mass));

[Na,Nd,N] = size(Mass);
Error_local = zeros(N,1);
Error_auto_Mass = zeros(N,1);
Error_auto_Missing = zeros(N,1);
Error_global = zeros(N,1);
Global_est_mass = zeros(N,1);

for n=1:N
    A = abs(Mass(:,:,n) - Missing_Mass(:,:,n)); % Mass exceeds x% of Missing Mass
    %A(A<0) = Inf;
    A(isnan(A)) = Inf;
    
    % Restrict tensor A to the cases where the estimated Diameter is close to the real
    % one
    %% Keep, for example, 10% Top most similar Diameter region.
    
    B = abs((Diameter(:,:,n) - 2*shell_candidates{n}.a)/(2*shell_candidates{n}.a));
    [~, index] = sort(B(:),'ascend');
    
    index = index(round(Top_Diam_Diff*length(index)):end);
    A(ind2sub(size(A),index)) = Inf;
    
    
    % find minimum difference between mass and missing mass
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
end

masses = [];
for n=1:N
    %masses = [masses; shell_candidates{n}.MassMin, shell_candidates{n}.MassMax, local_auto{n}.Missing, local_auto{n}.Mass];
    masses = [masses; shell_candidates{n}.MassMin, local_auto{n}.Missing, local_auto{n}.Mass];
    x{n} = shell_candidates{n}.name;
end
figure
bar1 = bar(masses);
set(gca, 'XTick', 1:N, 'XTickLabel', x);
ax = gca; 
ax.XTickLabelRotation = 45;
set(bar1(1),'DisplayName','By Hand (MassMin)');
set(bar1(2),'DisplayName','Algorithm (Mass)');
set(bar1(3),'DisplayName','Algorithm (MissingMass)');

rmse_Mass = sqrt(nanmean((masses(:,1)-masses(:,3)).^2));
rmse_Miss = sqrt(nanmean((masses(:,1)-masses(:,2)).^2));

% Create legend
legend('show');
title(['Top ',num2str(100*Top_Diam_Diff),'%,  rmse Missing=',num2str(rmse_Miss),', rmse Mass=', num2str(rmse_Mass)])


%% Scatter plot de MIssing Mass
figure
scatter(masses(:,1),masses(:,3),'DisplayName','Missing Mass')
hold on
%scatter(masses(:,2),masses(:,4),'DisplayName','Shell Mass')
% Create ylabel
ylabel({'Automatically estimated'});
% Create xlabel
xlabel({'Estimated by hand'});
ylim([0 300000]);
xlim([0 300000]);

rango_masas = 0:1000:300000;

plot(rango_masas, rango_masas,'k')
%plot(rango_masas, 0.6*rango_masas,'g')
%plot(rango_masas, 1.6*rango_masas,'g')





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