clear
clc
close all
shell_candidates = load_shells_v5();
%shell_candidates = shell_candidates(1:10);

N = size(shell_candidates,2);


%alpha = 0.1:0.01:0.5;
alpha = 0.4;
Nalpha = length(alpha);
%delta = 0:0.005:0.1;
delta = 0.0;

Ndelta = length(delta);

Mass = zeros(Nalpha,Ndelta,N); % HI-shell mass for every parameter (alpha,delta)
Missing_Mass = zeros(Nalpha,Ndelta,N); % Missing mass for every parameter (alpha,delta)
Diameter = zeros(Nalpha,Ndelta,N);

fig = figure('Position',[10,600,600,400]);

%parfor n=1:N
for n=11
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
    [Mass(:,:,n), Missing_Mass(:,:,n), Diameter(:,:,n)] = compute_mass_V5(cube,n_hdu,alpha,delta,shell);
    
end

save 'variables_nueva_missing.mat'

