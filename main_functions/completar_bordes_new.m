function [x, scan] = completar_bordes_new(scan,scan_max,long0,lat0,L,N)

scan = fill_in_external(scan, scan_max, long0, lat0, L, N);


P = size(scan,2);
ang = 2*pi/P;
x = zeros(P,2);

delta = L/N;

for p=1:P
    if ~scan{p}.flag % if missing
        nant=p;
        cond = 0;
        while (~scan{nant}.flag)&&(~cond) % busca anterior
            nant=nant-1;
            if nant == 0
                nant = P;
            end
            if nant == p
                cond = 1;
            end
                    
        end
        nsig=p;
        cond = 0;
        while (~scan{nsig}.flag)&&(~cond) % busca siguiente
            nsig=nsig+1;
            if nsig == P+1;
                nsig = 1;
            end
            if nsig == p;
                cond =1;
            end
        end
        
        if ~cond
            rini = sqrt((scan{nant}.long-long0)^2 + (scan{nant}.lat-lat0)^2);
            angn = nant*2*pi/P; 
            rend = sqrt((scan{nsig}.long-long0)^2 + (scan{nsig}.lat-lat0)^2);
            
            %for n=nant+1:nsig-1
            n = nant+1;
            if nant>nsig
                Npoints = 1:(P -(nant-nsig)-1);
            elseif nant<nsig
                Npoints = 1:(nsig - nant -1);
            else
                Npoints = 1:P-1;
            end
            
            m = (rend-rini)/length(Npoints);
            
            for n1=Npoints
                if n == P+1
                    n = 1;
                end
                %n
                scan{n}.flag = 1;
                %r = rini + m*(n-nant);
                r = rini + m*n1;
                scan{n}.long = long0 - r*cos(angn);
                scan{n}.lat = lat0 - r*sin(angn);
                scan{n}.ind = round(r/delta +1);
                %scan{n}.val = scan{n}.profile(scan{p}.ind);
                if scan{n}.ind <= length(scan{n}.profile)
                    scan{n}.val = scan{n}.profile(scan{n}.ind);
                else
                    scan{n}.val = scan{n}.profile(end);
                end
                angn = angn+2*pi/P;
                n = n + 1;
            end
        end
    end
    ang = ang+2*pi/P;
end

p2=1;
for n = 1:P
    if (scan{n}.flag)
        x(p2,1) = scan{n}.long;
        x(p2,2) = scan{n}.lat;
        p2 = p2 + 1;
    end
end

x = x(1:p2-1,:);


end

function [scan] = fill_in_external(scan, scan_max, long0, lat0, L, N)
% Compute the medium gap
P = size(scan,2);
gap = [];
for p = 1:P
    if (scan{p}.flag && scan_max{p}.flag) % if found in both, max and outwall
        gap = [gap,scan{p}.ind - scan_max{p}.ind];
    end
end
mean_gap = round(mean(gap));

% Complete missing outwall points with the mean
d = L/(N-1);
for p = 1:P
    ang = p*(2*pi/P);
    if (~scan{p}.flag && scan_max{p}.flag)
        scan{p}.flag = 1;
        scan{p}.ind = scan_max{p}.ind + mean_gap; % add mean_gap
        scan{p}.long = long0 - (scan{p}.ind - 1)*d*cos(ang);
        scan{p}.lat = lat0 - (scan{p}.ind - 1)*d*sin(ang);
        scan{p}.val = scan{p}.profile(min(scan{p}.ind,length(scan{p}.profile)));
    end
end



end