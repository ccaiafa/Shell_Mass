function [ mass, missing_mass, area , Temperatures] = compute_mass_V5( cube, header,  alpha, delta, shell, Visualization, perc)

%threshold = 3;
threshold = 5;
%threshold = 25;
%threshold = 15;
thresholddeltamax = 10;
%thresholddeltamax = 20;

%% Carga Cubo y sus parametros
I = header.PrimaryData.Size(1);
J = header.PrimaryData.Size(2);
K = header.PrimaryData.Size(3);

mass = zeros(length(alpha),length(delta));
missing_mass = zeros(length(alpha),length(delta));
area = zeros(length(alpha),length(delta));
Temperatures.img = zeros(length(alpha),length(delta));
Temperatures.shell = zeros(length(alpha),length(delta));
Temperatures.backg = zeros(length(alpha),length(delta));
Temperatures.miss = zeros(length(alpha),length(delta));

A.cube = cube;
clear cube;
%[A.long0,A.dlong,A.j0,A.lat0,A.dlat,A.i0,A.vel0,A.dvel,A.k0,A.longmax,A.latmax] = get_params(header.PrimaryData.Keywords);
[A.long0,A.dlong,A.j0,A.lat0,A.dlat,A.i0,A.vel0,A.dvel,A.k0] = get_params(header.PrimaryData.Keywords);

%% Generar suma de canales

L = 2*shell.a; % Tamano ventana
Imin = 0; Jmin=0; Imax = Inf; Jmax = Inf;
while (Imin < 1)||(Jmin < 1)||(Imax > size(A.cube,1))||(Jmax > size(A.cube,2))
    L = L*0.95;
    latmin = shell.lat - L;  Imin = round((latmin - A.lat0)/A.dlat + A.i0);
    latmax = shell.lat + L;  Imax = round((latmax - A.lat0)/A.dlat + A.i0);
    longmin = shell.long + L; Jmin = round((longmin - A.long0)/A.dlong + A.j0);
    longmax = shell.long - L; Jmax = round((longmax - A.long0)/A.dlong + A.j0);
end

%latmin = max(shell.lat - L, A.lat0);  Imin = round((latmin - A.lat0)/A.dlat + A.i0);
%latmax = min(shell.lat + L, A.latmax);  Imax = round((latmax - A.lat0)/A.dlat + A.i0);
%latmin = shell.lat - L;  Imin = round((latmin - A.lat0)/A.dlat + A.i0);
% if Imin<1
%     Imin =1;
% end
%latmax = shell.lat + L;  Imax = round((latmax - A.lat0)/A.dlat + A.i0);
% if Imax > size(A.cube,1)
%     Imax = size(A.cube,1);
% end
%longmin = shell.long + L; Jmin = round((longmin - A.long0)/A.dlong + A.j0);
% if Jmin<1
%     Jmin =1;
% end
%longmax = shell.long - L; Jmax = round((longmax - A.long0)/A.dlong + A.j0);
% if Jmax > size(A.cube,2)
%     Jmax = size(A.cube,2);
% end
dVred = shell.dV*perc;
vmin = shell.V0 - dVred/2; Kmin = round((vmin - A.vel0)/A.dvel + A.k0);
vmax = shell.V0 + dVred/2; Kmax = round((vmax - A.vel0)/A.dvel + A.k0);

Im = A.cube(Imin:Imax,Jmin:Jmax,Kmin:Kmax);
nI = Imax - Imin +1;
nJ = Jmax - Jmin + 1;
nK = Kmax - Kmin +1;

sub = nanmean(Im,3); % Promedio de todos los canales [V0-dV, V0+dV]
sub(isnan(sub)) = 0;

%% Generar perfiles radiales
P = 360; % numero de perfiles
N = 1000; % numero de puntos por perfil

%x = zeros(P,2); % maximos
%r = zeros(P,1); % radios

dif = L/N;

x0 = shell.long;
y0 = shell.lat;

i0 = round((y0 - A.lat0)/A.dlat + A.i0) - Imin + 1;
j0 =  round((x0 - A.long0)/A.dlong + A.j0) - Jmin + 1;

%fig2 = figure;
%axes2 = axes('Parent',fig2,'Layer','top');
%hold on
ang = 2*pi/P;
for p=1:P
        
    xend = x0 - L*cos(ang);
    yend = y0 - L*sin(ang);
    %plot(axes1,xend,yend,'MarkerSize',12,'Marker','.','LineWidth',1,'LineStyle','none','Color',[0 0 0]);
        
    iend = round((yend - A.lat0)/A.dlat + A.i0) - Imin + 1;
    jend = round((xend - A.long0)/A.dlong + A.j0) - Jmin + 1;
    scan_max{p}.profile = improfile(sub,[j0,jend],[i0,iend],N);
    
    %[scan_max{p}.val,scan_max{p}.ind] = max(scan_max{p}.profile);
    [scan_max{p}.val,scan_max{p}.ind] = local_max(scan_max{p}.profile,threshold,thresholddeltamax,0,L,shell);
    
    if scan_max{p}.ind 
        scan_max{p}.long = shell.long - (scan_max{p}.ind - 1)*dif*cos(ang);
        scan_max{p}.lat = shell.lat - (scan_max{p}.ind - 1)*dif*sin(ang);
        scan_max{p}.flag = 1;
        x(p,1) = scan_max{p}.long;
        x(p,2) = scan_max{p}.lat;
    else
        scan_max{p}.flag = 0;
        x(p,1) = NaN;
        x(p,2) = NaN;
    end
    %r(p) = sqrt((scan_max{p}.long - x0)^2 + (scan_max{p}.lat - y0)^2);
    
    %plot(axes1,scan_max{p}.long,scan_max{p}.lat,'MarkerSize',12,'Marker','*','LineWidth',1,'LineStyle','none','Color',[1 1 1]);
    %plot(axes1,shell.long + L*cos(ang),shell.lat+ L*sin(ang),'MarkerSize',12,'Marker','.','LineWidth',1,'LineStyle','none','Color',[0 0 0]);
    %plot profile
    
    %plot(axes2,scan_max{p}.profile);
    
    ang = ang+2*pi/P;
end

% delete rows with NaN
x(any(isnan(x),2),:)=[];

%fig1 = figure('Position',[10,500,1300,400]);
%hold on0
%axes1 = axes('Parent',fig1,'Layer','top');

if Visualization
    fig1 = gcf;
    subplot(2,2,1);
    %imagesc([longmin longmax],[latmin latmax],sub,'Parent',axes1)
    imagesc([longmin longmax],[latmin latmax],sub)
    cref = caxis;
    %image([longmin longmax],[latmin latmax],sub,'Parent',axes1,'CDataMapping','scaled')
    set(gca,'YDir','normal')
    set(gca,'XDir','reverse')
    box(gca,'on');
    hold(gca,'all');
    %plot(shell.long,shell.lat,'MarkerSize',12,'Marker','o','LineWidth',1,'LineStyle','none','Color',[1 0 0]);
else
    fig1 = [];
end

%fig2 = figure;

%fig3 = figure;

%fig4 = figure;
            

% Compute mass for every alpha, delta parameter value
inda =1;
scan_max0 = scan_max;
for a = alpha
    par_a = a;
    x_new = [];
    while size(x_new,1) < 2
        [scan_max_new,x_new,dev] = depurar_puntos_V5(scan_max,x,fig1,1,longmin,longmax,latmin,latmax,sub,shell,P,par_a,[1,0,0],Visualization);
        par_a = par_a*1.25;
    end
    
    
    %scan_max = scan_max_new;
    %x = x_new;
    
   
    if isempty(x_new)%|| (abs(2*dev(1)-shell.a)/shell.a > 0.5)||imag(dev(1))
        mass(inda,:) = NaN;
        missing_mass(inda,:) = NaN;
        area(inda,:) = NaN;
    else

        %% Buscar paredes
        im = round((shell.lat - A.lat0)/A.dlat + A.i0);
        jm = round((shell.long - A.long0)/A.dlong + A.j0);
        km = round((shell.V0 - A.vel0)/A.dvel + A.k0);
        Minimo = A.cube(im,jm,km);

        indd = 1;
        for d = delta
            inwall = zeros(P,2); % inner walls
            outwall = zeros(P,2); % outer walls
            [scan_in,scan_out,inwall,outwall] = buscar_paredes1(scan_max_new,inwall,outwall,d,shell,L);
            if size(outwall,1) < 1
                mass(inda,indd) = NaN;
                missing_mass(inda,indd) = NaN;
                area(inda,indd) = NaN;
            else
                area(inda,indd) = pi*dev(1)*dev(2);
                if Visualization
                    show_map(fig1,1,[longmin longmax],[latmin latmax], sub, cref)
                    graficar_puntos(fig1,1,x_new,'.',[1 0 0])
                end
                %if ~isempty(inwall)
                %    graficar_puntos(fig1,1,inwall,'x',[0 0 0])
                %end
                %graficar_puntos(fig1,1,outwall,'x',[0 0 1])
                % depurar bordes externos using 0.25*std (alpha=0.25)
                outwall_new = [];
                scan_out_new = [];
                d_par = 0.25;
                while size(outwall_new,2) < 2
                    [scan_out_new,outwall_new] = depurar_puntos_V5(scan_out,outwall,fig1,1,longmin,longmax,latmin,latmax,sub,shell,P,d_par,[0,0,0],0);
                    d_par = d_par*1.25;
                end
                scan_out = scan_out_new;
                
                % Constrain shell width to 0.5 degrees from maxima (Marcelo recommendation)
                [outwall, scan_out] = reduce_width(scan_out,scan_max_new,shell.long,shell.lat,L,N); 
                
                % Completar bordes externos
                [outwall, scan_out] = completar_bordes_new(scan_out,scan_max_new,shell.long,shell.lat,L,N); 
                if Visualization
                    %graficar_puntos(fig1,1,outwall,'x',[0 0 0])
                    plot(outwall(:,1),outwall(:,2),'Color',[1 1 1]) 
                    % mostrar elipse de catalogo
                    [~,~] = ellipse(shell.a,shell.b,shell.tita*180/pi,shell.X0,shell.Y0,[0,0,0],P,Visualization);
                    title(['SHELL ', shell.name,'MAP (alpha= ',num2str(a),'  delta=',num2str(d),')'])
                    set(gca,'fontsize', 9)
                end
                
                %
                %hold on
               
                %[xn,yn] = ellipse(shell.a,shell.b,shell.tita*180/pi,shell.X0,shell.Y0,[1 1 1],100);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 show_map(fig1,2,[longmin longmax],[latmin latmax], sub)
%                 graficar_puntos(fig1,2,x_new,'o',[1 0 0])
%                 if ~isempty(inwall)
%                     graficar_puntos(fig1,2,inwall,'x',[0 0 0])
%                 end
%                 %graficar_puntos(fig2,outwall,'x',[0 0 1])
%                 plot(outwall(:,1),outwall(:,2))

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Restar Local Background

                %sub_back = background(sub,scan_out,i0,j0,A,shell,L,P,N,Minimo);
                sub_back_flat = background_with_missing(sub,scan_out,i0,j0,A,shell,L,P,N,Minimo);

                if Visualization
                    %show_map(fig1,3,[longmin longmax],[latmin latmax], sub_back, cref)
                    show_map(fig1,3,[longmin longmax],[latmin latmax], sub_back_flat, cref)
                    hold on
                    plot(outwall(:,1),outwall(:,2),'Color',[1 1 1])
                    title(['BACKGROUND FLAT'])
                    set(gca,'fontsize', 9)
                end
                
                %sub_corrected = sub - sub_back;
                sub_corrected = sub - sub_back_flat;
                sub_corrected(sub_corrected < 0) = 0;

                if Visualization
                    show_map(fig1,2,[longmin longmax],[latmin latmax], sub_corrected, [])
                    hold on
                end
                
                %plot(outwall(:,1),outwall(:,2),'Color',[1 1 1])
                
                

                %% Calculo de masa
                C = 1.823e18;
                %dVT = (Kmax-Kmin+1)*A.dvel/1000; % en Km/s
                dVT = shell.dV/1000; % en Km/s
                %Omega = 4.759647e-6 ; % sr
                Omega = abs(A.dlong*A.dlat)*((pi/180)^2);
                constd = shell.D*3.085647758e21; % cm
                %constd = 2.8*3.085647758e21; % cm Distancia usada en el paper de Laura
                mHI = 1.67e-24;

                Masa = sub_corrected*C*dVT*Omega*(constd^2)*mHI; %gr
                Masa = sum(Masa(:)); %gr

                if Masa
                    Masa = Masa/1.989e33; % Masa solar
                else
                    Masa = NaN;
                end
                
                if Visualization
                    title([' Mass =',num2str(Masa,'%10.1e\n'),', ', ...
                        ' Ref=',num2str(shell.MassShell,'%10.1e\n'),...
                        ' Error=',num2str(100*(shell.MassShell-Masa)/(shell.MassShell+Masa)),'%'])
                    pause(0.01)
                    set(gca,'fontsize', 9)
                end

                mass(inda,indd) = Masa;
                Temperatures.img(inda,indd) = mean(sub(:)); 
                %Temperatures.shell(inda,indd) = mean(sub_corrected(:));
                
                sub_corrected(sub_corrected==0)=NaN;
                
                
                Temperatures.backg(inda,indd) = mean(sub_back_flat(:)); 
                Temperatures.shell(inda,indd) = nanmean(sub_corrected(:));
                %%
                %sub_back_with_missing = background_with_missing(sub,scan_out,i0,j0,A,shell,L,P,N,Minimo);
                
                %sub_missing = sub_back_with_missing - sub;
                sub_missing = sub_back_flat - sub;
                
                %sub_picos = -sub_missing(sub_missing<0); 
                sub_missing(sub_missing < 0) = 0;
                
                if Visualization
                    show_map(fig1,4,[longmin longmax],[latmin latmax], sub_missing, cref)
                    hold on
                    title(['MISSING'])
                    set(gca,'fontsize', 9)
                end
               
                
                %% Calculo de missing mass y masas picos
                C = 1.823e18;
                %dVT = (Kmax-Kmin+1)*A.dvel/1000; % en Km/s
                dVT = shell.dV/1000; % en Km/s
                %Omega = 4.759647e-6 ; % sr
                Omega = abs(A.dlong*A.dlat)*((pi/180)^2);
                constd = shell.D*3.085647758e21; % cm
                %constd = 2.8*3.085647758e21; % cm Distancia usada en el paper de Laura
                mHI = 1.67e-24;

                Masa_missing = sub_missing*C*dVT*Omega*(constd^2)*mHI; %gr
                Masa_missing = sum(Masa_missing(:)); %gr

                if Masa_missing
                    Masa_missing = Masa_missing/1.989e33; % Masa solar
                else
                    Masa_missing = NaN;
                end
                
                if Visualization
                    title(['Missing Mass= ',num2str(Masa_missing,'%10.1e\n'),...
                        ', Ref=',num2str(shell.MassMiss,'%10.1e\n'),...
                        ', Error=',num2str(100*(shell.MassMiss-Masa_missing)/(shell.MassMiss+Masa_missing)),'%'])
                    set(gca,'fontsize', 9)
                    pause(0.01)
         
                end

                missing_mass(inda,indd) = Masa_missing;
                Temperatures.miss(inda,indd) = mean(sub_missing(:)); 
                
                
                
            end

            indd = indd + 1;
            %disp(['alpha=',num2str(a),', delta=',num2str(d)])
        end
    end
    inda = inda + 1;
end

end

function  [val,ind] = local_max(profile,threshold,thresholddeltamax,thresholdmax,L,shell)
L0 = shell.b*0.5; % start searching at xx% of parameter b (ellipse semiaxis)
Lmax = shell.a*1.25; % stop searching at 125% of parameter a (ellipse semiaxis)
n2found = 0;
Npoints = length(profile);
dL = L/Npoints;
n0 = ceil(L0/dL);
nmax = min(ceil(Lmax/dL),Npoints);
%n0 = 1;
%nmax = Npoints;
for n2=n0:nmax
    if n2 > 1 
        if ((profile(n2) - profile(n2-1))/dL > threshold) && ~n2found;%missing
            n2found = n2;
        end
    end    
end

missingmax = 1;
minimum = profile(1);
if n2found; % if a bigslop was found at ang then we search for a maximum beyond that point
    n2 = n2found;
    while n2 < nmax && missingmax         
        ninf = n2 - 1;
        nsup = n2 + 1;                        
        while (ninf > 1) && abs(profile(n2) - profile(ninf))/dL < thresholddeltamax
            ninf = ninf -1;
        end
        while (nsup < nmax) && abs(profile(n2) - profile(nsup))/dL < thresholddeltamax
            nsup = nsup + 1;
        end

        nmean = round(ninf + (nsup - ninf)/2);
        if (profile(nmean) - minimum > thresholdmax) && ((profile(nmean) - profile(ninf))/dL >thresholddeltamax) && ((profile(nmean) - profile(nsup))/dL > thresholddeltamax) && (missingmax)         
            missingmax = 0;
        end
        n2 = n2 + 1;
    end
    
    if missingmax
        ind = 0;
        val = 0;
    else
        ind = n2-1;
        val = profile(n2);
    end
else
    ind = 0;
    val = 0;
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
if ~isempty(cref)
    caxis(cref)
end
colorbar()
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = graficar_puntos(fig,subplotnr,x,mark,color)
figure(fig)
subplot(2,2,subplotnr);
hold on
scatter(x(:,1), x(:,2),12,color,mark);


end


