%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = graficar_puntos(fig,subplotnr,x,mark,color)
figure(fig)
subplot(1,3,subplotnr);
hold on
scatter(x(:,1), x(:,2),12,color,mark);


end