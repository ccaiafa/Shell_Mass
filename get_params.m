%function [l0,dl,j0,b0,db,i0,v0,dv,k0,lmax,bmax] = get_params(Keywords)
function [l0,dl,j0,b0,db,i0,v0,dv,k0,lmax,bmax] = get_params(Keywords)
for n=1:size(Keywords,1)
    
    switch Keywords{n}
    case 'CRVAL1'
        l0 = Keywords{n,2};
    case 'CDELT1'
        dl = Keywords{n,2};
    case 'CRPIX1'
        j0 = Keywords{n,2};
    case 'CRVAL2'
        b0 = Keywords{n,2};
    case 'CDELT2'
        db = Keywords{n,2};
    case 'CRPIX2'
        i0 = Keywords{n,2};
    case 'CRVAL3'
        v0 = Keywords{n,2};
    case 'CDELT3'
        dv = Keywords{n,2};
    case 'CRPIX3'
        k0 = Keywords{n,2}; 
    end
end

%lmax = round(l0 - dl*j0); 
%bmax = round(b0 - db*i0); 

end