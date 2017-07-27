

function [alpha_opt, delta_opt, mean_error_train, std_error_train, mean_error_test, std_error_test] = hyper_tuning(dataRootPath, shell_train, alpha_range, delta_range, perc_out, Ncross, Visualization)

N = size(shell_train,2);
N_eval = round(perc_out*N);
N_subtrain = N-N_eval;

Nalpha = length(alpha_range);
Ndelta = length(delta_range);

% Crossvalidation loop
for n=1:Ncross
    % In each iteration we divide the training set into: subtrain_set and eval_set
    ind_rnd = randperm(N);
    ind_subtrain = ind_rnd(1:N_subtrain);
    ind_eval = ind_rnd(N_subtrain+1:end);

    shell_subtrain = shell_train(ind_subtrain);
    shell_eval = shell_train(ind_eval);

    Mass = zeros(Nalpha,Ndelta,N_subtrain); % HI-shell mass for every parameter (alpha,delta)
    Missing_Mass = zeros(Nalpha,Ndelta,N_subtrain); % Missing mass for every parameter (alpha,delta)
    Diameter = zeros(Nalpha,Ndelta,N_subtrain);
    
    if Visualization
        fig = figure('Position',[10,600,600,400]);
    end
    for s=1:N_subtrain
        shell = shell_subtrain{s};
        aux = shell.long;
        if (aux < 130 - shell.a) && ( aux > 80 + shell.a)
            disp(['shell ',num2str(s), ' en cubo 1']);
            disp(['loading data-cube...']);
            tic
            if shell.lat < 0 % Latitudes negativas
                fname = deblank(ls(fullfile(dataRootPath,'BL80-130.B10-50.FITS')));
            else             % Latitudes positivas
                fname = deblank(ls(fullfile(dataRootPath,'BL80-130.B50-10.FITS')));
            end
            cube = fitsread(fname);
            n_hdu = fitsinfo(fname);
            toc
        elseif (aux < 170 - shell.a) && (aux >120 + shell.a)
            disp(['shell ',num2str(s), 'en cubo 2']);
            disp(['loading data-cube...']);
            tic
            if shell.lat < 0 % Latitudes negativas
                disp('cubo 2 latitudes negativas no disponible')
                %fname = deblank(ls(fullfile(dataRootPath,'BL80-130.B10-50.FITS'));
            else             % Latitudes positivas
                fname = deblank(ls(fullfile(dataRootPath,'BL170-120.B-10_50.FITS')));
            end
            cube = fitsread(fname);
            n_hdu = fitsinfo(fname);
            toc
        elseif (aux < 210 - shell.a) && (aux >160 + shell.a)
            disp(['shell ',num2str(s), 'en cubo 3']);
            disp(['loading data-cube...']);
            tic
            if shell.lat < 0 % Latitudes negativas
                fname = deblank(ls(fullfile(dataRootPath,'BL160-210B-50.FITS')));
            else             % Latitudes positivas
                fname = deblank(ls(fullfile(dataRootPath,'BL160-210B+50.FITS')));
            end
            cube = fitsread(fname);
            n_hdu = fitsinfo(fname);
            toc            
        elseif (aux < 250 - shell.a) && (aux >200 + shell.a)
            disp(['shell ',num2str(s), 'en cubo 4']);
            disp(['loading data-cube...']);
            tic
            if shell.lat < 0 % Latitudes negativas
                fname = deblank(ls(fullfile(dataRootPath,'BL200-250.B-50.FITS')));
            else             % Latitudes positivas
                fname = deblank(ls(fullfile(dataRootPath,'BL200-250.B+50.FITS')));
            end
            cube = fitsread(fname);
            n_hdu = fitsinfo(fname);
            toc                 
        elseif (aux < 290 - shell.a) && (aux >240 + shell.a)
            disp(['shell ',num2str(s), 'en cubo 5']);
            disp(['loading data-cube...']);
            tic
            if shell.lat < 0 % Latitudes negativas
                fname = deblank(ls(fullfile(dataRootPath,'BL240-290.B-50.FITS')));
            else             % Latitudes positivas
                disp('cubo 5 latitudes positivas no disponible')
                %fname = deblank(ls(fullfile(dataRootPath,'BL200-250.B+50.FITS'));
            end
            cube = fitsread(fname);
            n_hdu = fitsinfo(fname);
            toc              
        end
        
        % Comute mass associated to HI-shell n
        disp(['Computing HI-Shell Mass: ', shell.name, ' ',num2str(s),'/',num2str(N_subtrain),'  ...']);
        [Mass(:,:,s), Missing_Mass(:,:,s), Diameter(:,:,s)] = compute_mass_V5(cube,n_hdu,alpha_range,delta_range,shell,Visualization);

    
    end

    
end





end