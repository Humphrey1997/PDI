%新的参数设置，统一用105GHz作为频率归一化

Bs = 3.34;
Te = 650;
Ti = 650;
ne = 3.09e19;
f0 = 105e9;
F0 = 105e9;

ep0 = 8.8542e-12;
c = 3e8;
me = 9.10956e-31;
e = 1.6e-19;
mratio = 2*1836; % D plasma
mi = me*mratio;

global wpe wpi wce wci w0;
W0 = 2*pi*F0;
w0 = f0/F0;
wpe = sqrt(ne*e^2/ep0/me)/W0;
wpi = sqrt(ne*e^2/ep0/mi)/W0;
wce = e*Bs/me/W0;
wci = e*Bs/mi/W0;
wLH = sqrt((wce^2*wci^2+wpi^2*wce^2)/(wpe^2+wce^2));

global ve vi
ve = sqrt(2*Te*e/me)/c;
vi = sqrt(2*Ti*e/mi)/c;


global nsp n1p n2p n0p
global nsz n1z n2z n0z
global des de1 de2 de0
n0p = 4.64;
n0z = 1e-5;
nsp = 4.85;
nsz = 1e-6;
de0 = 0;
de1 = 180/180*pi;%+180
des = de1-asin(n0p/nsp*sin(de1));
n1p = sqrt(nsp^2+n0p^2-2*nsp*n0p*cos(des));
n2p = sqrt(nsp^2+n0p^2+2*nsp*n0p*cos(des));
de2 = asin(nsp/n2p*sin(des));
n1z = nsz-n0z;
n2z = nsz+n0z;

global beta0 e0x e0y e0z;
E0 = 3e6;
beta0 = E0*(e/me/(ve*c)/W0);


d0 = D0(n0p);
[Vec,Dig] = eig(d0);
e0x = Vec(1,1);
e0y = Vec(2,1);
e0z = Vec(3,1);

vci = wci*1e3;