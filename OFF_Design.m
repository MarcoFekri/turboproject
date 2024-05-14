clc
clear all 
close all
%% from design problem
U=250;
m_dot_des=50;
Pt0=10^5;
Tt0=288;
A=[ 0.2801 0.3085 0.2801 0.3085];
alpha1_des=27.3508;
n_stage=1;
pi_c_des=3;
i_star=-1.703;
delta_star=12.76;
epsi_star=21.64;
%% calculations
N_rel=0.5:0.1:1.4;
mdot_rel=0.1:0.005:1.4;
for count1=1:length(N_rel)
    U_off=U*N_rel(count1);
    Limt1=0;
    for count2=1:length(mdot_rel)
        mdot=mdot_rel(count2)*m_dot_des;
        [Pt1,Tt1,alpha1,P1,M1,P0,T0]=igvsoff(Pt0,Tt0,mdot,A(1),alpha1_des);
        for n=1:n_stage
            [Pt3,Tt3,alpha3,max_eff,ince_off(2*n),ince_off(2*n+1)]=offdesign(Pt1,Tt1,mdot,U_off,A(2*n),A(2*n+1),alpha1);
            Pt1=Pt3;
            Tt1=Tt3;
            alpha1=alpha3;
        end
        mfpout_off=mdot*sqrt(Tt3)/Pt3;
        pic_off=Pt3/Pt0;
        toic_off=Tt3/Tt0;
        etac_off=(pic_off^(.4/1.4)-1)/(toic_off-1);
        % Applying The Condition That Efficiency is Between 0 & 1 %
        if (etac_off > 1) || (etac_off < 0) || (max_eff == 1)
        else
            pic_rel(count1,count2)=pic_off/pi_c_des;
            etac_rel(count1,count2)=etac_off;
            % To Get The surge points %
            if pic_rel(count1,count2) > Limt1
                picmax(count1)=pic_rel(count1,count2);
                mdmin(count1)=mdot_rel(count2);
                beg(count1)=count2;
                Limt1=pic_rel(count1,count2);
            end
            mdrr(count1,count2)=mdot_rel(count2);
            en(count1)=count2;
        end
        
    end
end

%% plots

for count1=1:length(N_rel)
    plot(mdrr(count1,beg(count1):en(count1)),pic_rel(count1,beg(count1):en(count1)),'color',rand(1,3),'LineWidth',3)
    hold on
end
surge=polyfit(mdmin,picmax,5);
plot(mdot_rel,polyval(surge,mdot_rel),'b-.','LineWidth',2)
axis([0.2 0.6 0.2 0.6])
legend('At Part Speed = 50%','At Part Speed = 60%','At Part Speed = 70%','At Part Speed = 80%','At Part Speed = 90%','At Part Speed = 100%','At Part Speed = 110%','At Part Speed = 120%','At Part Speed = 130%','At Part Speed = 140%','The Surge Line','The Operating Line')
title('Compressor Map','FontWeight','bold','color','b')
xlabel('MFP)_{rel}','FontWeight','bold','color','r')
ylabel('\pi_{c})_{rel}','FontWeight','bold','color','r')
grid on

figure(2)
for count1=1:length(N_rel)
    xi = mdrr(count1,beg(count1):en(count1));
    zi = etac_rel(count1,beg(count1):en(count1));
    plot(xi,zi,'color',rand(1,3),'LineWidth',3)
    hold on
end
legend('At Part Speed = 50%','At Part Speed = 60%','At Part Speed = 70%','At Part Speed = 80%','At Part Speed = 90%','At Part Speed = 100%','At Part Speed = 110%','At Part Speed = 120%','At Part Speed = 130%','At Part Speed = 140%','The Surge Line')
title('Efficiency Variation','FontWeight','bold','color','r')
xlabel('MFP)_{rel}','FontWeight','bold','color','b')
ylabel('\eta_{c})_{rel}','FontWeight','bold','color','b')
grid on
