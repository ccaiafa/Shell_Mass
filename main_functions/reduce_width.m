function [x, scan] = reduce_width(scan,scan_max,long0,lat0,L,N,dR,Reff)
max_width = dR*Reff;

scan = fill_in_external(scan, scan_max, long0, lat0, L, N);

P = size(scan,2);
ang = 2*pi/P;
delta = L/N;
x = zeros(P,2);

for p=1:P
    if scan{p}.flag % if not missing
        %if dist_grad(scan{p}, scan_max{p}) > 0.5
        d = (scan{p}.ind - scan_max{p}.ind)*delta;
        if d > max_width
            r = scan_max{p}.ind*delta + max_width ;
            scan{p}.long = long0 - r*cos(ang);
            scan{p}.lat = lat0 - r*sin(ang);
            scan{p}.ind = round(r/delta +1);
            if scan{p}.ind <= length(scan{p}.profile)
                scan{p}.val = scan{p}.profile(scan{p}.ind);
            else
                scan{p}.val = scan{p}.profile(end);
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