a = [1.2170 7.4413 6.8185 1.8419 1.7020 2.0497]';
b = [0.1442 -2.4138 -2.1632 -0.09367 -0.01485 -0.05168]';
c = [-0.007552 0.3103 0.2887 0.02001 0.01522 0.01807]';
d = [0 -0.01222 -0.01162 0 0 0]';
tita_s = [40 275 280 280 280 280]';
tita_e = [250 620 625 500 500 405]';


Nang = 720;
step = 720/Nang;
theta = linspace(0,720,Nang);
thetaradians = deg2rad(theta);
lnr = a + b*thetaradians + c*(thetaradians.^2) + d*(thetaradians.^3);

figure
%hold on    
ind_s = round(tita_s/step);
ind_e = round(tita_e/step);
for n=1:size(lnr,1)
     polarplot(thetaradians(ind_s(n):ind_e(n)), lnr(n,ind_s(n):ind_e(n)))
     %plot(ind_s(n):ind_e(n), lnr(n,ind_s(n):ind_e(n)))
    hold on
end
