function [x, scan] = completar_bordes(scan,long0,lat0,L,N)
P = size(scan,2);
ang = 2*pi/P;
x = zeros(P,2);

delta = L/N;

% eliminar puntos aislados
% for n = 1:P
%     if n==1 
%         nant = P; 
%     else
%         nant = n-1; 
%     end
%     if n==P
%         nsig = 1;
%     else
%         nsig = n+1;
%     end
%     if (scan{n}.flag)&&(~scan{nant}.flag)&&(~scan{nsig}.flag)
%         scan{n}.flag = 0;
%     end
% end


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