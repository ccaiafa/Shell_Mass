function [fname, cube, n_hdu] = select_cube(dataRootPath,shell)
aux = shell.long;
if (aux < 130 - shell.a) && ( aux > 80 + shell.a)
    disp(['shell ',' en cubo 1']);
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
    disp(['shell ', 'en cubo 3']);
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
    disp(['shell ', 'en cubo 4']);
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
    disp(['shell ', 'en cubo 5']);
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