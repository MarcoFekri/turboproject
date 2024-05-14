function [Pt3,Tt3,alpha3,max_eff,inced1,inced2]=offdesign(Pt1,Tt1,mdot,u1,A1,A2,alpha1)
options = optimset('Display','off');
%% Fitting of The Losses Curve
x_inc=[-.8 -.6 -.4 -.2 0 .2 .4 .6 .8];
y_loss=[.47 .3 .19 .1 .07 .1 .19 .3 .47];
w_loss=polyfit(x_inc,y_loss,2);
%% Fitting of The Deviation Curve
x_inc=[0 0.2  0.4 0.6 0.8];
y_dev=[0 0.05 0.2 0.4 0.7];
d_dev=polyfit(x_inc,y_dev,2);
%% At design
delta_des = 12.76; % [deg] eps=beta2-beta2d
eps_des = 23.1926;   % [deg] eps=beta1-beta2 
i_star=0.2266;      % [deg] same incidence angle for stator and rotor due to symmetry 

%% Across Rotor
beta1d=41.7; % [deg]
beta2d=16.7; % [deg]
[M1]=fsolve(@(M1) (mdot*sqrt(287*Tt1))/(Pt1*A1*cos(alpha1*pi/180))-(sqrt(1.4)*M1*(1+0.2*M1^2)^-3),.5,options);
T1=Tt1/(1+.2*M1^2);
P1=Pt1/((1+.2*M1^2)^3.5);
c1=M1*sqrt(1.4*287*T1);
cx1=c1*cos(alpha1*pi/180);
cth1=c1*sind(alpha1);
wth1=u1-cth1;
beta1=atand(wth1/cx1);       %[deg]
i_r=beta1-beta1d;
inc_r=(i_r-i_star)/(eps_des);
max_eff=0;
% In Order To Give Indication That The Incidence Exceed The Staling limits
% (I Assume The Limits Are That of The Range Of The Given Curves) 
if (inc_r > .8) || (inc_r < -.8)
    max_eff=1;
end
if inc_r <=0
    dev_r=0;
else
dev_r=polyval(d_dev,inc_r);
end
deltar=dev_r*eps_des+delta_des; 
beta2=(beta2d+deltar);          
wr=polyval(w_loss,inc_r);
w1=cx1/cos(beta1*pi/180);
Ttr1=T1+(w1^2/(2*1004.5));
Ptr1=P1*((Ttr1/T1)^3.5);
Ptr2=Ptr1*(1-wr)+P1*wr;
Ttr2=Ttr1;                     %constant across roter
[Mrel2]=fsolve(@(Mrel2) (mdot*sqrt(287*Ttr2))/(Ptr2*A2*cos(beta2*pi/180))-(sqrt(1.4)*Mrel2*((1+.2*Mrel2^2)^-3)),0.4,options);
T2=Ttr2/(1+.2*Mrel2^2);
P2=Ptr2/((1+.2*Mrel2^2)^3.5);
w2=Mrel2*sqrt(1.4*287*T2);
cx2= w2*cos(beta2*pi/180);
wth2=w2*sin(beta2*pi/180);
cth2=u1-wth2;
c2=sqrt(cx2^2+cth2^2);
M2=c2/sqrt(1.4*287*T2);
Tt2=T2*(1+.2*M2^2);
Pt2=P2*((1+.2*M2^2)^3.5);
%% Across Stator:
alpha2d=50.6;
alpha3d=14.67;
Tt3=Tt2;                      %Since No Heat added and no work across the rotor
alpha2=atan(cth2/cx2)*(180/pi);       % [deg]
i_s=alpha2-alpha2d;
inc_s=(i_s-i_star)/(eps_des);
if (inc_s > .8) || (inc_s < -.8)
    max_eff=1;
end
if inc_s <=0
    dev_s=0;
else
dev_s=polyval(d_dev,inc_s);
end
deltas=dev_s*eps_des+delta_des;%[deg]
alpha3=alpha3d+deltas; % [deg]
ws=polyval(w_loss,inc_s);
Pt3=Pt2*(1-ws)+ws*P2;
inced1=(beta1-beta1d);   %[deg] Incidence to rotor
inced2=(alpha2-alpha2d); %[deg] Incidence to stator