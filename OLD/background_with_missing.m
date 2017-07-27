function sub_back = background_with_missing(sub,scan,i0,j0,A,shell,L,P,Npts,Minimo)

%P = size(scan,2);
[I,J] = size(sub);
sub_back = zeros(I,J);
N = zeros(I,J);
step = 0.05;

delta = L/(Npts-1);

long0 = shell.long;
lat0 = shell.lat;

latmin = shell.lat - L;  Imin = round((latmin - A.lat0)/A.dlat + A.i0);
% if Imin<1
%     Imin =1;
% end
latmax = shell.lat + L;  %Imax = round((latmax - A.lat0)/A.dlat + A.i0);
longmin = shell.long + L; Jmin = round((longmin - A.long0)/A.dlong + A.j0);
% if Jmin<1
%     Jmin =1;
% end
longmax = shell.long - L; %Jmax = round((longmax - A.long0)/A.dlong + A.j0);

ang = 2*pi/P;
for p=1:round(P/2)-1
    pp=p;
    cond = 0;
    while (~scan{pp}.flag)&&(~cond)
        if pp==1
            pp = P;
        else
            pp = pp - 1;
        end
        if pp == p
            cond=1;
        end
    end
    
    if ~cond
        rho = scan{pp}.val;
        R = (scan{pp}.ind-1)*delta;
        pp_op = round(mod(pp+P/2,P));
        rho_op = scan{pp_op}.val;
        R_op = (scan{pp_op}.ind-1)*delta;
        
        m = (rho - rho_op)/(R+R_op);
        val_0 = rho_op + m*R_op;
        
        % move towards R
        i = i0;
        j = j0; 
        r = 0;
        while (i<=I)&&(i>=1)&&(j<=J)&&(j>=1)
            if r < R
                sub_back(i,j) = sub_back(i,j) + val_0 + m*r;
                %disp(num2str(val_0 + m*r));
                N(i,j) = N(i,j) + 1;
            else
                %sub_back(i,j) = sub(i,j);
                %plot(long,lat,'MarkerSize',12,'Marker','*','LineWidth',1,'LineStyle','none','Color',[1 1 1]);
            end

            iant = i;
            jant = j;

            while (i==iant)&&(j==jant)
                r = r + step;
                long = long0 - r*cos(ang);
                lat = lat0 - r*sin(ang);
                i = round((lat - A.lat0)/A.dlat + A.i0) - Imin + 1;
                j = round((long - A.long0)/A.dlong + A.j0) - Jmin + 1;
                %disp([num2str(i),' ',num2str(j)])
            end
            
        end
        
        % move towards R_op
        i = i0;
        j = j0; 
        r = 0;
        while (i<=I)&&(i>=1)&&(j<=J)&&(j>=1)
            if r < R_op           
                sub_back(i,j) = sub_back(i,j) + val_0 - m*r;
                %disp(num2str(val_0 - m*r));
                N(i,j) = N(i,j) + 1;
            else
                %sub_back(i,j) = sub(i,j);
                %plot(long,lat,'MarkerSize',12,'Marker','*','LineWidth',1,'LineStyle','none','Color',[1 1 1]);
            end

            iant = i;
            jant = j;

            while (i==iant)&&(j==jant)
                r = r + step;
                long = long0 - r*cos(ang+pi);
                lat = lat0 - r*sin(ang+pi);
                i = round((lat - A.lat0)/A.dlat + A.i0) - Imin + 1;
                j = round((long - A.long0)/A.dlong + A.j0) - Jmin + 1;
                %disp([num2str(i),' ',num2str(j)])
            end        
        end
    end
    
    
    ang = ang+2*pi/P;
end

for i=1:I
    for j=1:J
        if (N(i,j))
            sub_back(i,j) = sub_back(i,j)/N(i,j);
        else
            sub_back(i,j) = sub(i,j);
        end
    end
end



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function val = curve(m,rho,R,r)
% 
% val = m + (R - sqrt(R^2 - r^2))*(rho+m)/R;
% 
% end

function val = curve(m,rho,R,r)

val = rho -(rho-m)*sqrt(1-(r/R)^2);

end
