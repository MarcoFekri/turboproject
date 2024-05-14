function [Pt1,Tt1,alpha1,P1,M1,P0,T0]=igvsoff(Pt0,Tt0,md,A,alpha1_des)
y_loss=0.03;                         % offdesign loss coefficient for IGV
options = optimset('Display','off');
alpha1=alpha1_des;                 %For inlet guide vane the losses remains constant for M1<=0.5
[Mo]=fsolve(@(Mo) (md*sqrt(287*Tt0))/(Pt0*A*cosd(alpha1))-(sqrt(1.4)*Mo/((1+.2*Mo^2)^3)),.35,options);
P0=Pt0/((1+.2*Mo^2)^3.5);
T0=Tt0/(1+.2*Mo^2);
Tt1=Tt0;
% To calc the pressure substitute in loss equation for IGV
[x]=fsolve(@(x) [x(1)-x(2)*(1+.2*x(3)^2)^3.5; y_loss-(Pt0-x(1))/(x(1)-x(2));...
md*sqrt(287*Tt1)/(x(1)*A*cos(alpha1*pi/180))-sqrt(1.4)*x(3)*(1+.2*x(3)^2)^(-3)],[10^5 9*10^4 .5],options);
Pt1=x(1);
P1=x(2);
M1=x(3);
