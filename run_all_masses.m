
function [] = run_all_masses(dataRootPath, fun, dataOutPath, alpha_range, delta_range, Visualization)

% Load full dataset
%shell_all = load_shells_4c_Dvel_100();
shell_all = fun();
N = size(shell_all,2);

Nalpha = length(alpha_range);
Ndelta = length(delta_range);

Mass = zeros(Nalpha,Ndelta,N); % HI-shell mass for every parameter (alpha,delta)
Missing_Mass = zeros(Nalpha,Ndelta,N); % Missing mass for every parameter (alpha,delta)
Area = zeros(Nalpha,Ndelta,N);

if Visualization
    fig = figure('Position',[10,600,600,400]);
end

%% Compute masses 
disp('Computing masses')
parfor s=1:N
%for s=59:N
    %fprintf('.')
    disp([num2str(s),'/',num2str(N)])
    shell = shell_all{s};
    [fname, cube, n_hdu] = select_cube(dataRootPath,shell);
    
    % Comute mass associated to HI-shell n
    %disp(['Computing HI-Shell Mass: ', shell.name, ' ',num2str(s),'/',num2str(N_subtrain),'  ...']);
    [Mass(:,:,s), Missing_Mass(:,:,s), Area(:,:,s)] = compute_mass_V5(cube,n_hdu,alpha_range,delta_range,shell,Visualization);
end

save('variables7.mat')
Diff = abs(Mass - Missing_Mass);


end

function [fname, cube, n_hdu] = select_cube(dataRootPath,shell)
aux = shell.long;
if (aux < 130 - shell.a) && ( aux > 80 + shell.a)
    %disp(['shell ',num2str(s), ' en cubo 1']);
    %disp(['loading data-cube...']);
    %tic
    if shell.lat < 0 % Latitudes negativas
        fname = deblank(ls(fullfile(dataRootPath,'BL80-130.B10-50.FITS')));
    else             % Latitudes positivas
        fname = deblank(ls(fullfile(dataRootPath,'BL80-130.B50-10.FITS')));
    end
    cube = fitsread(fname);
    n_hdu = fitsinfo(fname);
    %toc
elseif (aux < 170 - shell.a) && (aux >120 + shell.a)
    disp(['shell ', 'en cubo 2']);
    %disp(['loading data-cube...']);
    %tic
    if shell.lat < 0 % Latitudes negativas
        fname = deblank(ls(fullfile(dataRootPath,'BL170-120.B-10_50.FITS')));
    else             % Latitudes positivas
        fname = deblank(ls(fullfile(dataRootPath,'BL170-120.B50-10.FITS')));
    end
    cube = fitsread(fname);
    n_hdu = fitsinfo(fname);
    %toc
elseif (aux < 210 - shell.a) && (aux >160 + shell.a)
    %disp(['shell ',num2str(s), 'en cubo 3']);
    %disp(['loading data-cube...']);
    %tic
    if shell.lat < 0 % Latitudes negativas
        fname = deblank(ls(fullfile(dataRootPath,'BL160-210B-50.FITS')));
    else             % Latitudes positivas
        fname = deblank(ls(fullfile(dataRootPath,'BL160-210B+50.FITS')));
    end
    cube = fitsread(fname);
    n_hdu = fitsinfo(fname);
    %toc
elseif (aux < 250 - shell.a) && (aux >200 + shell.a)
    %disp(['shell ',num2str(s), 'en cubo 4']);
    %disp(['loading data-cube...']);
    %tic
    if shell.lat < 0 % Latitudes negativas
        fname = deblank(ls(fullfile(dataRootPath,'BL200-250.B-50.FITS')));
    else             % Latitudes positivas
        fname = deblank(ls(fullfile(dataRootPath,'BL200-250.B+50.FITS')));
    end
    cube = fitsread(fname);
    n_hdu = fitsinfo(fname);
    %toc
elseif (aux < 290 - shell.a) && (aux >240 + shell.a)
    %disp(['shell ',num2str(s), 'en cubo 5']);
    %disp(['loading data-cube...']);
    %tic
    if shell.lat < 0 % Latitudes negativas
        fname = deblank(ls(fullfile(dataRootPath,'BL240-290.B-50.FITS')));
    else             % Latitudes positivas
        fname = deblank(ls(fullfile(dataRootPath,'BL240-290.B+50.FITS')));
    end
    cube = fitsread(fname);
    n_hdu = fitsinfo(fname);
    %toc
end
end



