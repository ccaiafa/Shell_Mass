clear
%load 'variables-2.mat'
load 'variables-2.mat'
%load 'variables_sin_constraint.mat'
%close all
%
%Mass(:,1:10,:) = NaN;
%Missing_Mass(:,1:10,:) = NaN;
%Mass(1:10,:,:) = NaN;
%Missing_Mass(1:10,:,:) = NaN;

Diff2 = zeros(size(Mass));
Diff_missing = zeros(size(Mass));

[Na,Nd,N] = size(Mass);
Error_local = zeros(N,1);
Error_auto_Mass = zeros(N,1);
Error_auto_Missing = zeros(N,1);
Error_global = zeros(N,1);
Global_est_mass = zeros(N,1);

for n=1:N
%for n=11
    Diff2(:,:,n) = (Mass(:,:,n) - repmat(shell_candidates{n}.MassMax,[Na,Nd])).^2;
    [val,ind] = min(reshape(Diff2(:,:,n),[Na*Nd,1]));
    [ia,id] = ind2sub([Na,Nd],ind);
    local_opt{n}.alpha = alpha(ia);
    local_opt{n}.delta = delta(id);
    local_opt{n}.ia = ia;
    local_opt{n}.id = id;
    local_opt{n}.Mass = Mass(ia,id,n);
    local_opt{n}.Missing = Missing_Mass(ia,id,n);
    
    A = ((1-0.0)*Mass(:,:,n) - Missing_Mass(:,:,n)); % Mass exceeds x% of Missing Mass
    A(A<0) = Inf;
    A(isnan(A)) = Inf;
    %diam = Diameter(:,:,n);
    %A(abs(diam/2 - ones(size(diam))*shell_candidates{n}.a)/shell_candidates{n}.a >0.5) = Inf;
    
    % find minimum difference between mass and missing mass
    Diff_missing(:,:,n) = A;
    [val,ind] = min(reshape(Diff_missing(:,:,n),[Na*Nd,1]));
    [ia,id] = ind2sub([Na,Nd],ind);
    
    % find minimum difference in diameter
%     diam = Diameter(:,:,n);
%     diam_dif = abs(diam/2 - ones(size(diam))*shell_candidates{n}.a);
%     [val,ind] = min(reshape(diam_dif,[Na*Nd,1]));
%     [ia,id] = ind2sub([Na,Nd],ind);
    
    local_auto{n}.alpha = alpha(ia);
    local_auto{n}.delta = delta(id);
    local_auto{n}.ia = ia;
    local_auto{n}.id = id;
    local_auto{n}.Mass = Mass(ia,id,n);
    local_auto{n}.Missing = Missing_Mass(ia,id,n);
    %local_auto{n}.diameter = shell_candidates{n}.a;
        
end
cost_function = nansum(Diff2,3);
[val,ind] = nanmin(reshape(cost_function,[Na*Nd,1]));
[ia,id] = ind2sub([Na,Nd],ind);
global_opt.alpha = alpha(ia);
global_opt.delta = delta(id);
global_opt.ia = ia;
global_opt.id = id;

for n=1:N
    Global_est_mass(n) = Mass(ia,id,n);
    Error_global(n) = sqrt(Diff2(global_opt.ia,global_opt.id,n))/shell_candidates{n}.MassMax;
    Error_local(n) = sqrt(Diff2(local_opt{n}.ia,local_opt{n}.id,n))/shell_candidates{n}.MassMax;
    %Error_auto(n) = Diff_missing(local_auto{n}.ia,local_auto{n}.id,n)/Mass(local_auto{n}.ia,local_auto{n}.id,n);
    Error_auto_Mass(n) = abs(Mass(local_auto{n}.ia,local_auto{n}.id,n) - shell_candidates{n}.MassMax)/shell_candidates{n}.MassMax;
    Error_auto_Missing(n) = abs(Missing_Mass(local_auto{n}.ia,local_auto{n}.id,n) - shell_candidates{n}.MassMin)/shell_candidates{n}.MassMin;
end

masses = [];
for n=1:N
    %masses = [masses; shell_candidates{n}.MassMin, shell_candidates{n}.MassMax, local_auto{n}.Missing, local_auto{n}.Mass];
    masses = [masses; shell_candidates{n}.MassMin, local_auto{n}.Missing];
    x{n} = shell_candidates{n}.name;
end
bar1 = bar(masses);
set(gca, 'XTick', 1:N, 'XTickLabel', x);
ax = gca; 
ax.XTickLabelRotation = 45;
set(bar1(1),'DisplayName','By Hand');
set(bar1(2),'DisplayName','Algorithm');
%set(bar1(2),'DisplayName','MassMax');
%set(bar1(3),'DisplayName','Closest Estimated Mass');
%set(bar1(4),'DisplayName','Associated Missing Mass');
%set(bar1(3),'DisplayName','Automaticaly Estimated Missing Mass');
%set(bar1(4),'DisplayName','Automaticaly Estimated Mass');
% Create legend
legend('show');
title('Missing Mass')



% Scatter plot de MIssing Mass
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

