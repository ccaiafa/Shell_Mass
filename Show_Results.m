clear
clc
%close all


%shell_candidates = load_shells_v5();
%shell_candidates = shell_candidates(7:9);

%N = size(shell_candidates,2);


%alpha = 0.1:0.01:0.3;
%alpha = 0.16;
%Nalpha = length(alpha);
%delta = 0:0.005:0.1;
%delta = 0.07;

%Ndelta = length(delta);

%Mass = zeros(Nalpha,Ndelta,N); % HI-shell mass for every parameter (alpha,delta)
%Missing_Mass = zeros(Nalpha,Ndelta,N); % Missing mass for every parameter (alpha,delta)

load 'variables-2.mat'
load results

fig = figure('Position',[10,600,600,400]);

%for n=1:N
for n=23
    clf(fig)
    shell = shell_candidates{n};
    aux = shell.long;
    if (aux < 130 - shell.a) && ( aux > 80 + shell.a)
        disp(['shell ',num2str(n), ' en cubo 1']);
        disp(['loading data-cube...']);
        tic
        if shell.lat < 0 % Latitudes negativas
            fname = '/Volumes/HD-PCTU3/DATASET/Astronomy/HI/Cubo80-130/BL80-130.B10-50.FITS';
            %fname = '/Users/CesarMB13/Desktop/Cubo80-130/BL80-130.B10-50.FITS';
            %fname = 'Cubo80-130/BL80-130.B10-50.FITS';
        else             % Latitudes positivas
            fname = '/Volumes/HD-PCTU3/DATASET/Astronomy/HI/Cubo80-130/BL80-130.B50-10.FITS';
            %fname = '/Users/CesarMB13/Desktop/Cubo80-130/BL80-130.B50-10.FITS';
            %fname = 'Cubo80-130/BL80-130.B50-10.FITS';
        end
        cube = fitsread(fname);
        n_hdu = fitsinfo(fname);
        toc                
    elseif (aux < 170 - shell.a) && (aux >120 + shell.a)
        disp(['shell ',num2str(n), 'en cubo 2']);
    elseif (aux < 210 - shell.a) && (aux >160 + shell.a)
        disp(['shell ',num2str(n), 'en cubo 3']);
    elseif (aux < 250 - shell.a) && (aux >200 + shell.a) 
        disp(['shell ',num2str(n), 'en cubo 4']);
    elseif (aux < 290 - shell.a) && (aux >240 + shell.a)
        disp(['shell ',num2str(n), 'en cubo 5']);
    end
    
    % Comute mass associated to HI-shell n
    disp(['Computing HI-Shell Mass: ', shell.name, ' ',num2str(n),'/',num2str(N),'  ...']);
    
    
    alpha_aut = local_auto{n}.alpha;
    delta_aut = local_auto{n}.delta;
    ia = local_auto{n}.ia;
    id = local_auto{n}.id;
    
    %alpha = local_opt{n}.alpha;
    %delta = local_opt{n}.delta;
    
    [Mass_auto, Missing_Mass_aut] = compute_mass_V5(cube,n_hdu,alpha_aut,delta_aut,shell);
    
    Mass_alphas = Mass(:,id,n);
    Mass_deltas = Mass(ia,:,n);
    
    Missing_alphas = Missing_Mass(:,id,n);
    Missing_deltas = Missing_Mass(ia,:,n);
    
    figure
    subplot(2,2,1);
    %histogram(Mass_alphas,100)
    plot(alpha,Mass_alphas)
    title('Shell Mass vs alpha')
    
    subplot(2,2,2);
    %histogram(Mass_deltas,100)
    plot(delta,Mass_deltas)
    title('Shell Mass vs delta')
    
    subplot(2,2,3);
    %histogram(Missing_alphas,100)
    plot(alpha,Missing_alphas)
    title('Missing Mass vs alpha')
    
    subplot(2,2,4);
    %histogram(Missing_deltas,100)
    plot(delta,Missing_deltas)
    title('Missing Mass vs delta')
    
    pause(1)
    close
    
    
    
end


