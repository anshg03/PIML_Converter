%% buck_control_design.m
% This file is for the controller design of the buck converter stage
% calculation of input impedance


clear all;
close all;
%%
P=bodeoptions;
P.FreqUnits='Hz';
P.PhaseWrapping='on';
P.Grid='on';
w=logspace(1,5,10000);
s=tf('s');
%[]
%operating point
% R0=8;
% Vdc0=48;
% Vin0=24;
% Uout0=24;
% iLf0=Vdc0/R0;
% Iload0=Vdc0/R0;
% uRef0=Vdc0;
% D0=Vin0/Vdc0;
% Imax=100;
P0= 186;
Vdc0=48;
Vin0=24;
Iload0=P0/Vdc0;
uRef0=Vin0;
Uout0=Vin0;
D0=Vin0/Vdc0;
iLf0=Iload0/(1-D0);
Imax=2*iLf0;
% filter paramters
fs=20e3;              % switching frequency
fsamp=20e3;          % sampling frequency
Tsamp=1/fsamp;
Lf=1/fs*1/(0.4*iLf0)*(1-0.5)*Vdc0;     
%Cf=1/(Lf*(0.1*2*pi*fs)^2);
Cf=(Iload0 * (1-D0))/(fs * 0.8);

Cdc=163e-6;
rCdc=0.02;
Ldc=0;
rLdc=0;

% total plant delay
%Td=0.5*1/fs+0.5*1/fsamp;   % delay of double update mode + fast sampling
Td=0.75/fs;
[num, den]=pade(Td,2);  % approximation for dead-time
Gt=tf(num, den);

% current controller design
Gi=-1/(Lf*s)*Gt;

wBWi=2*pi*fs*0.1;
Tni=5*1/(wBWi);
kpi=wBWi*Lf;
Ri=kpi*(1+s*Tni)/(s*Tni);

Giol=-Ri*Gi;
Gicl=feedback(Giol,1);

figure();
bode(Gi,w,P);
grid on;
hold on;
bode(Ri,w,P);
bode(Giol,w,P);
legend('G_{iL}','R_{i}','R_{i}*G_{iL}');
title('compensated open-loop tf of converter current');

% voltage controller design
Gu=Vdc0/iLf0*(1-s*Lf*iLf0/(D0*Vdc0))/(1+s*Cdc*Vdc0/(D0*iLf0))*Gicl;

wBWu=0.2*wBWi;
kpu=wBWu*Cdc/D0;
Tnu=10*1/(wBWu);
Ru=kpu*(1+s*Tnu)/(s*Tnu);

Guol=Ru*Gu;

figure();
bode(Gu,w,P);
grid on;
hold on;
bode(Ru,w,P);
bode(Guol,w,P);
legend('G_{uC}','R_{u}','R_{u}*G_{uC}');
title('compensated open-loop tf of converter output voltage');


kmff=0;kiff=0;
Zout=-(((Lf*s*uRef0 + Gt*Ri*(iLf0*kiff*Lf*s + uRef0 - kiff*uRef0))*Vdc0^2)/(uRef0*((-uRef0)*(Gt*iLf0*(Ri + kmff*Lf*s) + uRef0 - Gt*kmff*uRef0) + ...
     Gt*Ri*Ru*(iLf0*Lf*s - uRef0)*Vdc0 - Cdc*s*(Gt*Ri + Lf*s)*Vdc0^2)));
kmff=1;kuff=0;
Zoutmff=-(((Lf*s*uRef0 + Gt*Ri*(iLf0*kiff*Lf*s + uRef0 - kiff*uRef0))*Vdc0^2)/(uRef0*((-uRef0)*(Gt*iLf0*(Ri + kmff*Lf*s) + uRef0 - Gt*kmff*uRef0) + ...
     Gt*Ri*Ru*(iLf0*Lf*s - uRef0)*Vdc0 - Cdc*s*(Gt*Ri + Lf*s)*Vdc0^2)));
kmff=0;kiff=1;
Zoutiff=-(((Lf*s*uRef0 + Gt*Ri*(iLf0*kiff*Lf*s + uRef0 - kiff*uRef0))*Vdc0^2)/(uRef0*((-uRef0)*(Gt*iLf0*(Ri + kmff*Lf*s) + uRef0 - Gt*kmff*uRef0) + ...
     Gt*Ri*Ru*(iLf0*Lf*s - uRef0)*Vdc0 - Cdc*s*(Gt*Ri + Lf*s)*Vdc0^2)));
kmff=1;kiff=1;
Zoutmffiff=-(((Lf*s*uRef0 + Gt*Ri*(iLf0*kiff*Lf*s + uRef0 - kiff*uRef0))*Vdc0^2)/(uRef0*((-uRef0)*(Gt*iLf0*(Ri + kmff*Lf*s) + uRef0 - Gt*kmff*uRef0) + ...
     Gt*Ri*Ru*(iLf0*Lf*s - uRef0)*Vdc0 - Cdc*s*(Gt*Ri + Lf*s)*Vdc0^2)));

figure();
bode(Zout,w,P);
grid on;
hold on;
% bode(Zoutmff,w,P);
bode(Zoutiff,w,P);
% bode(Zoutmffiff,w,P);
legend('Z_{out}','Z_{out,iff}');
title('Output impedance of converter');


((1 + Cdc*rCdc*s)*(Lf*s*uRef0 + Gt*Ri*(iLf0*kiff*Lf*s + uRef0 - kiff*uRef0))*Vdc0^2)/ ...
(uRef0*((1 + Cdc*rCdc*s)*uRef0*(Gt*iLf0*(Ri + kmff*Lf*s) + uRef0 - Gt*kmff*uRef0) + Gt*Ri*Ru*(1 + Cdc*rCdc*s)*((-iLf0)*Lf*s + uRef0)*Vdc0 + Cdc*s*(Gt*Ri + Lf*s)*Vdc0^2));

% 
%% effect of ESR DC capacitance
% kmff=0;kiff=0;
% Zout=((1 + Cdc*rCdc*s)*(Lf*s*uRef0 + Gt*Ri*(iLf0*kiff*Lf*s + uRef0 - kiff*uRef0))*Vdc0^2)/ ...
% (uRef0*((1 + Cdc*rCdc*s)*uRef0*(Gt*iLf0*(Ri + kmff*Lf*s) + uRef0 - Gt*kmff*uRef0) + Gt*Ri*Ru*(1 + Cdc*rCdc*s)*((-iLf0)*Lf*s + uRef0)*Vdc0 + Cdc*s*(Gt*Ri + Lf*s)*Vdc0^2));
% kmff=1;kiff=0;
% Zoutmff=((1 + Cdc*rCdc*s)*(Lf*s*uRef0 + Gt*Ri*(iLf0*kiff*Lf*s + uRef0 - kiff*uRef0))*Vdc0^2)/ ...
% (uRef0*((1 + Cdc*rCdc*s)*uRef0*(Gt*iLf0*(Ri + kmff*Lf*s) + uRef0 - Gt*kmff*uRef0) + Gt*Ri*Ru*(1 + Cdc*rCdc*s)*((-iLf0)*Lf*s + uRef0)*Vdc0 + Cdc*s*(Gt*Ri + Lf*s)*Vdc0^2));
% kmff=0;kiff=1;
% Zoutiff=((1 + Cdc*rCdc*s)*(Lf*s*uRef0 + Gt*Ri*(iLf0*kiff*Lf*s + uRef0 - kiff*uRef0))*Vdc0^2)/ ...
% (uRef0*((1 + Cdc*rCdc*s)*uRef0*(Gt*iLf0*(Ri + kmff*Lf*s) + uRef0 - Gt*kmff*uRef0) + Gt*Ri*Ru*(1 + Cdc*rCdc*s)*((-iLf0)*Lf*s + uRef0)*Vdc0 + Cdc*s*(Gt*Ri + Lf*s)*Vdc0^2));
% kmff=1;kiff=1;
% Zoutmffiff=((1 + Cdc*rCdc*s)*(Lf*s*uRef0 + Gt*Ri*(iLf0*kiff*Lf*s + uRef0 - kiff*uRef0))*Vdc0^2)/ ...
% (uRef0*((1 + Cdc*rCdc*s)*uRef0*(Gt*iLf0*(Ri + kmff*Lf*s) + uRef0 - Gt*kmff*uRef0) + Gt*Ri*Ru*(1 + Cdc*rCdc*s)*((-iLf0)*Lf*s + uRef0)*Vdc0 + Cdc*s*(Gt*Ri + Lf*s)*Vdc0^2));
% 
% figure();
% bode(Zout,w,P);
% grid on;
% hold on;
% bode(Zoutiff,w,P);
% 
% legend('Z_{out}','Z_{out,iff}');
% title('Output impedance incl. Cdc ESR=0.2');
% % 
% % 
% 
% % cascaded droop and active damping
% 
% kd=5/10; % dV/dI
% Glp=1/(1+s/(2*pi*100));
% Rad=0;
% 
% 
% kmff=0;kiff=0;
% ZoutDroop=((1 + Cdc*rCdc*s)*Vdc0*((-Gt)*(Glp*kd + Rad)*Ri*Ru*(iLf0*Lf*s - uRef0)*uRef0 + Lf*s*uRef0*Vdc0 + Gt*Ri*(iLf0*kiff*Lf*s + uRef0 - kiff*uRef0)*Vdc0))/ ...
%   (uRef0*((1 + Cdc*rCdc*s)*uRef0*(Gt*iLf0*(Ri + kmff*Lf*s) + uRef0 - Gt*kmff*uRef0) + Gt*Ri*Ru*(1 + Cdc*rCdc*s)*((-iLf0)*Lf*s + uRef0)*Vdc0 + Cdc*s*(Gt*Ri + Lf*s)*Vdc0^2));
% 
% Thp=1/(2*pi*100);
% Tlp=1/(2*pi*400);
% Rad=20/10*s*Thp/(1+s*Thp)*1/(1+s*Tlp);
% 
% kmff=0;kiff=0;
% ZoutDroopAd=((1 + Cdc*rCdc*s)*Vdc0*((-Gt)*(Glp*kd + Rad)*Ri*Ru*(iLf0*Lf*s - uRef0)*uRef0 + Lf*s*uRef0*Vdc0 + Gt*Ri*(iLf0*kiff*Lf*s + uRef0 - kiff*uRef0)*Vdc0))/ ...
%   (uRef0*((1 + Cdc*rCdc*s)*uRef0*(Gt*iLf0*(Ri + kmff*Lf*s) + uRef0 - Gt*kmff*uRef0) + Gt*Ri*Ru*(1 + Cdc*rCdc*s)*((-iLf0)*Lf*s + uRef0)*Vdc0 + Cdc*s*(Gt*Ri + Lf*s)*Vdc0^2));
% 
% Zoutcasc=ZoutDroopAd;
% 
% figure();
% bode(Zout,w,P);
% grid on;
% hold on;
% bode(ZoutDroop,w,P);
% bode(ZoutDroopAd,w,P);
% legend('Z_{out}','Z_{out,droop}','Z_{out,droop,ad}');
% title('Output impedance cascaded droop');
% 
% 
% % current droop and active damping
% rCdc=0.2
% kd=20/10; % dV/dI
% Glp=1;
% Gad=0;
% 
% kmff=0;kiff=0;
% Zoutweak=(kd*(Gt*Ri + Lf*s)*(1 + Cdc*rCdc*s)*uRef0*Vdc0^2)/(kd*(1 + Cdc*rCdc*s)*uRef0^2*(Gt*iLf0*(Ri + kmff*Lf*s) + uRef0 - Gt*kmff*uRef0) + ...
%    ((-Gt)*iLf0*(Glp + Gad*kd)*Lf*Ri*s*(1 + Cdc*rCdc*s) + (Gt*(Glp + Gad*kd)*Ri + Cdc*Gt*(kd + Glp*rCdc + Gad*kd*rCdc)*Ri*s + Cdc*kd*Lf*s^2)*uRef0)*Vdc0^2);
% 
% kd=5/10; % dV/dI
% Glp=1;
% Gad=0;
% 
% kmff=0;kiff=0;
% Zoutstrong=(kd*(Gt*Ri + Lf*s)*(1 + Cdc*rCdc*s)*uRef0*Vdc0^2)/(kd*(1 + Cdc*rCdc*s)*uRef0^2*(Gt*iLf0*(Ri + kmff*Lf*s) + uRef0 - Gt*kmff*uRef0) + ...
%    ((-Gt)*iLf0*(Glp + Gad*kd)*Lf*Ri*s*(1 + Cdc*rCdc*s) + (Gt*(Glp + Gad*kd)*Ri + Cdc*Gt*(kd + Glp*rCdc + Gad*kd*rCdc)*Ri*s + Cdc*kd*Lf*s^2)*uRef0)*Vdc0^2);
% 
% Glp=1/(1+s/(2*pi*100));
% Gad=0;
% 
% kmff=0;kiff=0;
% ZoutDroop=(kd*(Gt*Ri + Lf*s)*(1 + Cdc*rCdc*s)*uRef0*Vdc0^2)/(kd*(1 + Cdc*rCdc*s)*uRef0^2*(Gt*iLf0*(Ri + kmff*Lf*s) + uRef0 - Gt*kmff*uRef0) + ...
%    ((-Gt)*iLf0*(Glp + Gad*kd)*Lf*Ri*s*(1 + Cdc*rCdc*s) + (Gt*(Glp + Gad*kd)*Ri + Cdc*Gt*(kd + Glp*rCdc + Gad*kd*rCdc)*Ri*s + Cdc*kd*Lf*s^2)*uRef0)*Vdc0^2);
% 
% Thp=1/(2*pi*100);
% Tlp=1/(2*pi*10000);
% Gad=0.5*s*Thp/(1+s*Thp)*1/(1+s*Tlp);
% 
% kmff=0;kiff=0;
% ZoutDroopAd=(kd*(Gt*Ri + Lf*s)*(1 + Cdc*rCdc*s)*uRef0*Vdc0^2)/(kd*(1 + Cdc*rCdc*s)*uRef0^2*(Gt*iLf0*(Ri + kmff*Lf*s) + uRef0 - Gt*kmff*uRef0) + ...
%    ((-Gt)*iLf0*(Glp + Gad*kd)*Lf*Ri*s*(1 + Cdc*rCdc*s) + (Gt*(Glp + Gad*kd)*Ri + Cdc*Gt*(kd + Glp*rCdc + Gad*kd*rCdc)*Ri*s + Cdc*kd*Lf*s^2)*uRef0)*Vdc0^2);
% 
% 
% figure();
% bode(Zoutweak,w,P);
% grid on;
% hold on;
% bode(Zoutstrong,w,P);
% bode(ZoutDroop,w,P);
% bode(ZoutDroopAd,w,P);
% legend('Z_{out,kd2.0}','Z_{out,kd0.5}','Z_{out,LP}','Z_{out,LP+ad}');
% title('Output impedance current droop');
% 
% figure();
% bode(Zoutcasc,w,P);
% grid on;
% hold on;
% bode(ZoutDroopAd,w,P);
% legend('Z_{casc,ad}','Z_{direct,ad}');
% title('Comparison cascaded / current droop');
% 






