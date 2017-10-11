function [scan_in,scan_out,inwall,outwall] = buscar_paredes1(scan,inwall,outwall,delta,shell,L)
P = size(scan,2);
N = length(scan{1}.profile);

d = L/(N-1);

%ang = 2*pi/P;
for p=1:P
    ang = p*(2*pi/P);
    scan_in{p}.profile = scan{p}.profile;
    scan_out{p}.profile = scan{p}.profile;
    if scan{p}.flag
        %pared interna
        rangeint = 1:scan{p}.ind;
        scan_in{p}.ind = find_inwall(scan{p}.profile,rangeint,delta);
        if ~isempty(scan_in{p}.ind)
            scan_in{p}.flag = 1;
            scan_in{p}.long = shell.long - (scan_in{p}.ind - 1)*d*cos(ang);
            scan_in{p}.lat = shell.lat - (scan_in{p}.ind - 1)*d*sin(ang);
            scan_in{p}.val = scan{p}.profile(scan_in{p}.ind);
        else
            scan_in{p}.flag = 0;
        end
        
        % pared externa
        rangeext = scan{p}.ind:N;
        scan_out{p}.ind = find_outwall(scan{p}.profile,rangeext,delta);
        if ~isempty(scan_out{p}.ind)
            scan_out{p}.flag = 1;
            scan_out{p}.long = shell.long - (scan_out{p}.ind - 1)*d*cos(ang);
            scan_out{p}.lat = shell.lat - (scan_out{p}.ind - 1)*d*sin(ang);
            scan_out{p}.val = scan{p}.profile(scan_out{p}.ind);
        else
            scan_out{p}.flag = 0;
        end
    else
        scan_in{p}.flag = 0;
        scan_out{p}.flag = 0;
              
    end
    %ang = ang+2*pi/P;
end


p1=1;
p2=1;
for n = 1:size(scan,2)
    if (scan_in{n}.flag)
        inwall(p1,1) = scan_in{n}.long;
        inwall(p1,2) = scan_in{n}.lat;
        p1 = p1 + 1;
    end
    if (scan_out{n}.flag)
        outwall(p2,1) = scan_out{n}.long;
        outwall(p2,2) = scan_out{n}.lat;
        p2 = p2 + 1;
    end
end

if p1 > 1
    inwall = inwall(1:p1-1,:);
else
    inwall = [];
end
if p2 > 1
    outwall = outwall(1:p2-1,:);
else
    outwall = [];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Busca pared interna (se queda con la pared mas cercana al maximo)
function [i] = find_inwall(profile,rango,delta)
M = profile(rango(end)); % Maximo

%if profile(rango(1)) > M - delta;
if profile(rango(1)) > M - M*delta;
    i = []; % no existe pared
else
    i = rango(end);
    %while (profile(i) > M - delta)
    while (profile(i) > M - M*delta)
        i = i - 1;
    end
end
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Busca pared externa (se queda con la pared mas cercana al maximo)
function [i] = find_outwall(profile,rango,delta)
M = profile(rango(1)); % Maximo

i = rango(1);
%while (profile(i) > M - delta)&& (i < rango(end))
while (profile(i) > M - M*delta&& (i < rango(end)))
    i = i + 1;
end

if i == rango(end)
    i=[];
end
   
end