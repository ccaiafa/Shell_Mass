

function [alpha_opt_miss, delta_opt_miss, R2_miss_train, R2_miss_eval] = hyper_tuning(dataRootPath, shell_train, alpha_range, delta_range, perc_out, Ncross, Visualization)

N = size(shell_train,2);
N_eval = round(perc_out*N);
N_subtrain = N-N_eval;

Nalpha = length(alpha_range);
Ndelta = length(delta_range);


alpha_opt_miss = zeros(Ncross, 1);
delta_opt_miss = zeros(Ncross, 1);
R2_miss_train = zeros(Ncross, 1);
R2_miss_eval = zeros(Ncross, 1);

% Crossvalidation loop
for n=1:Ncross
    disp(' ')
    disp(['Crossvalidation iteration ',num2str(n)]);
    % In each iteration we divide the training set into: subtrain_set and eval_set
    ind_rnd = randperm(N);
    ind_subtrain = ind_rnd(1:N_subtrain);
    ind_eval = ind_rnd(N_subtrain+1:end);

    shell_subtrain = shell_train(ind_subtrain);
    shell_eval = shell_train(ind_eval);

    Mass = zeros(Nalpha,Ndelta,N_subtrain); % HI-shell mass for every parameter (alpha,delta)
    Missing_Mass = zeros(Nalpha,Ndelta,N_subtrain); % Missing mass for every parameter (alpha,delta)
    Diameter = zeros(Nalpha,Ndelta,N_subtrain);
    %Diff2_mass = zeros(Nalpha,Ndelta,N_subtrain);
    Diff2_miss = zeros(Nalpha,Ndelta,N_subtrain);  
    
    if Visualization
        fig = figure('Position',[10,600,600,400]);
    end

    %% Compute masses in the subtrain set
    fprintf('Training ')
    for s=1:N_subtrain
        fprintf('.')
        shell = shell_subtrain{s};
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
            %disp(['shell ',num2str(s), 'en cubo 2']);
            %disp(['loading data-cube...']);
            %tic
            if shell.lat < 0 % Latitudes negativas
                fname = deblank(ls(fullfile(dataRootPath,'BL170-120.B50-10.FITS')));
            else             % Latitudes positivas
                fname = deblank(ls(fullfile(dataRootPath,'BL170-120.B-10_50.FITS')));
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
        
        % Comute mass associated to HI-shell n
        %disp(['Computing HI-Shell Mass: ', shell.name, ' ',num2str(s),'/',num2str(N_subtrain),'  ...']);
        [Mass(:,:,s), Missing_Mass(:,:,s), Diameter(:,:,s)] = compute_mass_V5(cube,n_hdu,alpha_range,delta_range,shell,Visualization);
        %Diff2_mass(:,:,s) = abs(Mass(:,:,s) - shell.MassMin).^2; 
        Diff2_miss(:,:,s) = abs(Missing_Mass(:,:,s) - shell.MassMin).^2; 
    end
    % Find minimum error
    SSres_miss = nansum(Diff2_miss,3);
    
    SStot_miss = zeros(size(Missing_Mass));
    Mean_miss = nanmean(Missing_Mass,3);
    for s=1:N_subtrain
        SStot_miss(:,:,s) = (Missing_Mass(:,:,s) - Mean_miss).^2;
    end
    SStot_miss = nansum(SStot_miss,3);
    
    %cost_mass_function = nanmean(Diff2_mass,3);
    %cost_miss_function = nanmean(Diff2_miss,3);
    cost_miss_function = ones(size(SSres_miss)) - SSres_miss./SStot_miss;
    
%     [~,ind_mass] = nanmax(reshape(cost_mass_function,[Nalpha*Ndelta,1]));
%     [ia,id] = ind2sub([Nalpha,Ndelta],ind_mass);
%     alpha_opt_mass(n) = alpha_range(ia);
%     delta_opt_mass(n) = delta_range(id);
%     error_mass_train(n) = cost_mass_function(ia, id);
    
    [~,ind_miss] = nanmax(reshape(cost_miss_function,[Nalpha*Ndelta,1]));
    [ia,id] = ind2sub([Nalpha,Ndelta],ind_miss);
    alpha_opt_miss(n) = alpha_range(ia);
    delta_opt_miss(n) = delta_range(id);  
    R2_miss_train(n) = cost_miss_function(ia, id);
    disp(' ')
    disp(['(alpha, delta) = ', '(',num2str(alpha_opt_miss(n)),',',num2str(delta_opt_miss(n)),...
        '),  Training R2=',num2str(R2_miss_train(n))])
    
    
    %% Compute masses in the validation set
    %Diff2_mass_eval = zeros(N_eval,1);
    Diff2_miss_eval = zeros(N_eval,1); 
    %Mass_eval = zeros(N_eval,1);
    Missing_Mass_eval = zeros(N_eval,1);
    Diameter_eval = zeros(N_eval,1);
    
    disp(' '); fprintf('Evaluating ')
    for s=1:N_eval
        fprintf('.')
        shell = shell_eval{s};
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
            %disp(['shell ',num2str(s), 'en cubo 2']);
            %disp(['loading data-cube...']);
            %tic
            if shell.lat < 0 % Latitudes negativas
                fname = deblank(ls(fullfile(dataRootPath,'BL170-120.B50-10.FITS')));
            else             % Latitudes positivas
                fname = deblank(ls(fullfile(dataRootPath,'BL170-120.B-10_50.FITS')));
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
        
        % Comute mass associated to HI-shell n
        %disp(['Computing HI-Shell Mass: ', shell.name, ' ',num2str(s),'/',num2str(N_subtrain),'  ...']);
        [Mass_eval(s), Missing_Mass_eval(s), Diameter_eval(s)] = compute_mass_V5(cube,n_hdu,alpha_opt_miss(n),delta_opt_miss(n),shell,Visualization);
        %Diff2_mass_eval(s) = abs(Mass_eval(s) - shell.MassMin)^2;
        Diff2_miss_eval(s) = abs(Missing_Mass_eval(s) - shell.MassMin)^2;
    end
    
    SSres_miss_eval = nansum(Diff2_miss_eval);
    Mean_miss_eval = nanmean(Missing_Mass_eval);
    SStot_miss_eval = nansum((Missing_Mass_eval - Mean_miss_eval).^2);
    
    R2_miss_eval(n) = 1 - SSres_miss_eval/SStot_miss_eval;
    
    disp(' ')
    disp(['Evaluation R2=',num2str(R2_miss_eval(n))])
    
end





end