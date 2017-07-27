function [scan,x,devant] = depurar_puntos_V5(scan,x,fig1,subplotnr,longmin,longmax,latmin,latmax,sub,shell,P,alpha,color_puntos)

deltapuntos = 1;
while (deltapuntos ~= 0) && (size(x,1))>5 %(~isempty(x))
    %if ~isempty(x)
        [dev(1), dev(2), xc, yc, tita,error,error1] = ellipse_fit_opt_V2(x(:,1), x(:,2));
    
        if (imag(xc) || imag(yc) || (error == Inf))
            scan = clear_flags(scan);
            x = [];
            deltapuntos = 0;
        else
    
        center = mean(x);
        
        cref = caxis;
        show_map(fig1,subplotnr,[longmin longmax],[latmin latmax], sub, cref)
    
        graficar_puntos(fig1,subplotnr,x,'o',color_puntos)
        
    
        [xep(:,1),xep(:,2)] = ellipse(dev(1),dev(2),tita*180/pi,center(1),center(2),color_puntos,P);

        pause(0.1);
        
        devant = dev;
        [scan, x, deltapuntos ] = depurar(scan, xep, alpha*dev(1)); 
        end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [scan, x, count] = depurar(scan, xep, tol)
P = size(scan,2);

count = 0;
for n = 1:P
    %if sqrt((scan{n}.long-x0)^2+(scan{n}.lat-y0)^2) > 2.5*dev1 || sqrt((scan{n}.long-x0)^2+(scan{n}.lat-y0)^2) < 0.5*dev2
    if (scan{n}.flag) && mindist(scan{n},xep) > tol
        scan{n}.flag = 0;
        count = count +1;
    end
end

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


m = mean(xep);

p=1;
for n = 1:size(scan,2)
    if (scan{n}.flag)
        x(p,1) = scan{n}.long;
        x(p,2) = scan{n}.lat;
        p = p + 1;
    end
end

if p == 1
    x = [];
end
    

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dmin] = mindist(point, xep)

dmin = Inf;
for p = 1:size(xep,1)
    d = sqrt((point.long - xep(p,1))^2 + (point.lat - xep(p,2))^2);
    if d < dmin
        pmin = p;
        dmin = d;
    end

end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = show_map(fig,subplotnr,rangelong, rangelat, sub, cref)

%clf(fig)
% show map
%axes1 = axes('Parent',fig,'Layer','top');
subplot(2,2,subplotnr);
imagesc(rangelong,rangelat,sub)
%image([longmin longmax],[latmin latmax],sub,'Parent',axes1,'CDataMapping','scaled')
set(gca,'YDir','normal')
set(gca,'XDir','reverse')
box(gca,'on');
hold(gca,'all')
caxis(cref)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = graficar_puntos(fig,subplotnr,x,mark,color)
figure(fig)
subplot(2,2,subplotnr);
hold on
scatter(x(:,1), x(:,2),12,color,mark);


end

