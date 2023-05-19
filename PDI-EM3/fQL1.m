function res = fQL1(ws,n,r)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
global wpe wce ve;
global nsp n1p n0p;
global nsz n1z n0z;
global e0x e0y e0z;
global des de1 de0;
global beta0 w0;

wc = wce;
vt = ve;
wp = wpe;

w1 = ws-w0;
de = n*de1-r*des;
rn = r - n;

e0xc = conj(e0x);
e0yc = conj(e0y);
e0zc = conj(e0z);
e0p = e0x + 1i*e0y;
e0s = e0x - 1i*e0y;
e0pc = conj(e0p);
e0sc = conj(e0s);

as = nsp*vt/wc;
a1 = n1p*vt/wc;
a0 = n0p*vt/wc;
asx = as*cos(des);
asy = as*sin(des);
a1x = a1*cos(de1);
a1y = a1*sin(de1);

Eis = exp(1i*des);
Eisc = conj(Eis);
Ei1 = exp(1i*de1);
Ei1c = conj(Ei1);

cs = nsz*vt;
c1 = n1z*vt;
c0 = n0z*vt;
w1n = w1-n*wc;
wsr = ws-r*wc;
w0rn = w0-rn*wc;
zeta1n = w1n/c1;
zetasr = wsr/cs;
zeta0rn = w0rn/c0;

ZA012=g2(0,zeta1n,1,zetasr,2)/c1^1/cs^2;
ZA011=g2(0,zeta1n,1,zetasr,1)/c1^1/cs^1;
ZA112=g2(1,zeta1n,1,zetasr,2)/c1^1/cs^2;
ZA212=g2(2,zeta1n,1,zetasr,2)/c1^1/cs^2;
ZA111=g2(1,zeta1n,1,zetasr,1)/c1^1/cs^1;
ZA211=g2(2,zeta1n,1,zetasr,1)/c1^1/cs^1;
ZA311=g2(3,zeta1n,1,zetasr,1)/c1^1/cs^1;
ZAmx1=ZA012;ZAmx2=ZA011;ZAmx3=ZA111;
ZAmy1=ZA012;ZAmy2=ZA011;ZAmy3=ZA111;
ZAmz1=ZA112;ZAmz2=ZA011;ZAmz3=ZA111;ZAmz4=ZA211;
ZAzx1=ZA112;ZAzx2=ZA111;ZAzx3=ZA211;
ZAzy1=ZA112;ZAzy2=ZA111;ZAzy3=ZA211;
ZAzz1=ZA212;ZAzz2=ZA111;ZAzz3=ZA211;ZAzz4=ZA311;


ZB012=g2(0,zeta1n,1,zeta0rn,2)/c1^1/c0^2;
ZB011=g2(0,zeta1n,1,zeta0rn,1)/c1^1/c0^1;
ZB112=g2(1,zeta1n,1,zeta0rn,2)/c1^1/c0^2;
ZB212=g2(2,zeta1n,1,zeta0rn,2)/c1^1/c0^2;
ZB111=g2(1,zeta1n,1,zeta0rn,1)/c1^1/c0^1;
ZB211=g2(2,zeta1n,1,zeta0rn,1)/c1^1/c0^1;
ZB311=g2(3,zeta1n,1,zeta0rn,1)/c1^1/c0^1;
ZBmx1=ZB012;ZBmx2=ZB112;ZBmx3=ZB011;ZBmx4=ZB111;ZBmx5=ZB211;
ZBmy1=ZB012;ZBmy2=ZB112;ZBmy3=ZB011;ZBmy4=ZB111;ZBmy5=ZB211;
ZBmz1=ZB012;ZBmz2=ZB112;ZBmz3=ZB011;ZBmz4=ZB111;ZBmz5=ZB211;
ZBzx1=ZB112;ZBzx2=ZB212;ZBzx3=ZB111;ZBzx4=ZB211;ZBzx5=ZB311;
ZBzy1=ZB112;ZBzy2=ZB212;ZBzy3=ZB111;ZBzy4=ZB211;ZBzy5=ZB311;
ZBzz1=ZB112;ZBzz2=ZB212;ZBzz3=ZB111;ZBzz4=ZB211;ZBzz5=ZB311;

    function y = fxyz(x)
        jsr0 = besselj(r,as*x);
        jsr1 = besselj(r-1,as*x);
        jsr2 = besselj(r+1,as*x);
        jsr3 = besselj(r-2,as*x);
        jsr4 = besselj(r+2,as*x);
        
        j0rn0 = besselj(rn,a0*x);
        j0rn1 = besselj(rn-1,a0*x);
        j0rn2 = besselj(rn+1,a0*x);
        j0rn3 = besselj(rn-2,a0*x);
        j0rn4 = besselj(rn+2,a0*x);
        
        j1n0 = besselj(n,a1*x);
        j1n1 = besselj(n-1,a1*x);
        j1n2 = besselj(n+1,a1*x);
        
        ej1 = Eis*jsr1+Eisc*jsr2;
        ej2 = Eis*jsr1-Eisc*jsr2;
        ej3 = e0sc*j0rn1+e0pc*j0rn2;
        ej4 = e0sc*j0rn1-e0pc*j0rn2;
        ej5 = Eis*(jsr0-jsr3)+Eisc*(jsr0-jsr4);
        ej6 = Eis*(jsr0-jsr3)-Eisc*(jsr0-jsr4);
        ej7 = e0sc*(j0rn0-j0rn3)+e0pc*(j0rn0-j0rn4);
        ej8 = e0sc*(j0rn0-j0rn3)-e0pc*(j0rn0-j0rn4);
        
        YAmx1 = (1i*(-1/2)).*cs.*ej1.*w0.^(-1).*x.^3.*(2.*e0zc.*j0rn0.*w0rn+c0.*ej3.*x);
        YAmx2 = (1/4).*w0.^(-1).*x.^2.*(1i.*as.*ej3.*ej6.*w0.*x+2.*ej1.*((-1i).*(ej3+ej4.*r).*w0+2.*j0rn0.*(asy.*e0xc.*w0-e0yc.*(asx.*w0+a0.*r.*wc)).*x+((2i).*ej3.*w0+a0.*e0yc.*(1i.*asy.*(j0rn1-j0rn2)+asx.*(j0rn1+j0rn2)).*wc).*x.^2));
        YAmx3 = (1/4).*w0.^(-1).*x.*((-4i).*e0zc.*ej1.*j0rn0.*rn.*wc+(2i).*(c0.*ej1.*(ej3+ej4.*r)+e0zc.*(a0.*ej1.*(-j0rn1+j0rn2).*r+as.*ej6.*j0rn0.*rn).*wc).*x+((-1i).*as.*c0.*ej3.*ej6+4.*ej1.*j0rn0.*(-asy.*c0.*e0xc+asx.*c0.*e0yc+(2i).*e0zc.*w0+a0.*asy.*e0zc.*wc)).*x.^2);
        YBmx1 = (1i*(1/2)).*c0.*cs.*ej1.*ej3.*ws.^(-1).*x.^4;
        YBmx2 = 1i.*c0.*cs.*e0zc.*ej1.*j0rn0.*ws.^(-1).*x.^3;
        YBmx3 = (1/4).*ws.^(-1).*x.^2.*((2i).*ej3.*(ej1+ej2.*rn).*ws+((4i).*cs.*e0zc.*ej1.*j0rn0+4.*asy.*ej3.*jsr0.*rn.*wc+(-1i).*a0.*ej1.*ej8.*ws).*x+(-2).*ej1.*ej3.*(a0.*asy.*wc+(2i).*ws).*x.^2);
        YBmx4 = (1i*(1/4)).*ws.^(-1).*x.*(cs.*x.*((-2).*ej3.*(ej1+ej2.*rn)+a0.*ej1.*ej8.*x)+2.*e0zc.*(2.*ej2.*j0rn0.*rn.*ws+(2i).*asy.*j0rn0.*wc.*x.*((-2).*jsr0.*rn+a0.*ej1.*x)+ej1.*ws.*x.*(a0.*(j0rn1-j0rn2)+(-4).*j0rn0.*x)));
        YBmx5 = (1i*(-1/2)).*cs.*e0zc.*ws.^(-1).*x.*(2.*ej2.*j0rn0.*rn+a0.*ej1.*(j0rn1-j0rn2).*x);
        
        YAmy1 = (-1/2).*cs.*ej2.*w0.^(-1).*x.^3.*(2.*e0zc.*j0rn0.*w0rn+c0.*ej3.*x);
        YAmy2 = (1/4).*w0.^(-1).*x.^2.*(as.*ej3.*ej5.*w0.*x+2.*ej2.*(-(ej3+ej4.*r).*w0+(2i).*j0rn0.*(-asy.*e0xc.*w0+asx.*e0yc.*w0+a0.*e0yc.*r.*wc).*x+(2.*ej3.*w0+a0.*e0yc.*(asy.*(j0rn1-j0rn2)+(-1i).*asx.*(j0rn1+j0rn2)).*wc).*x.^2));
        YAmy3 = (1/4).*w0.^(-1).*x.*((-4).*e0zc.*ej2.*j0rn0.*rn.*wc+2.*(c0.*ej2.*(ej3+ej4.*r)+e0zc.*(a0.*ej2.*(-j0rn1+j0rn2).*r+as.*ej5.*j0rn0.*rn).*wc).*x-(as.*c0.*ej3.*ej5+(4i).*ej2.*j0rn0.*(-asy.*c0.*e0xc+asx.*c0.*e0yc+(2i).*e0zc.*w0+a0.*asy.*e0zc.*wc)).*x.^2);
        YBmy1 = (1/2).*c0.*cs.*ej2.*ej3.*ws.^(-1).*x.^4;
        YBmy2 = c0.*cs.*e0zc.*ej2.*j0rn0.*ws.^(-1).*x.^3;
        YBmy3 = (1/4).*ws.^(-1).*x.^2.*(2.*ej3.*(ej2+ej1.*rn).*ws+(4.*cs.*e0zc.*ej2.*j0rn0-a0.*ej2.*ej8.*ws+(-4).*ej3.*jsr0.*(asx.*rn.*wc+a0.*ws)).*x+2.*ej3.*(a0.*asx.*ej1.*wc+(-2).*ej2.*ws).*x.^2);
        YBmy4 = (1/4).*ws.^(-1).*x.*(4.*e0zc.*ej1.*j0rn0.*rn.*ws+(-2).*(cs.*ej3.*(ej2+ej1.*rn)+4.*asx.*e0zc.*j0rn0.*jsr0.*rn.*wc+a0.*e0zc.*(-ej2.*j0rn1+ej2.*j0rn2+4.*j0rn0.*jsr0).*ws).*x+(a0.*(cs.*ej2.*ej8+4.*cs.*ej3.*jsr0+4.*asx.*e0zc.*ej1.*j0rn0.*wc)+(-8).*e0zc.*ej2.*j0rn0.*ws).*x.^2);
        YBmy5 = (1/2).*cs.*e0zc.*ws.^(-1).*x.*((-2).*ej1.*j0rn0.*rn+a0.*(ej2.*(-j0rn1+j0rn2)+4.*j0rn0.*jsr0).*x);
        
        YAmz1 = (-1i).*cs.*jsr0.*w0.^(-1).*x.^2.*(2.*e0zc.*j0rn0.*w0rn+c0.*ej3.*x);
        YAmz2 = (-1i).*jsr0.*w0.^(-1).*x.^2.*(2.*e0zc.*j0rn0.*w0rn+c0.*ej3.*x);
        YAmz3 = (1/2).*w0.^(-1).*x.*((-2i).*ej4.*jsr0.*r.*w0+(4.*asy.*e0xc.*j0rn0.*jsr0.*w0+(-4).*asx.*e0yc.*j0rn0.*jsr0.*w0+(-1i).*as.*ej3.*(jsr1-jsr2).*w0+(-4).*a0.*e0yc.*j0rn0.*jsr0.*r.*wc).*x+2.*jsr0.*((2i).*ej3.*w0+a0.*e0yc.*(1i.*asy.*(j0rn1-j0rn2)+asx.*(j0rn1+j0rn2)).*wc).*x.^2);
        YAmz4 = (1/2).*w0.^(-1).*x.*((2i).*(c0.*ej4.*jsr0.*r+e0zc.*(a0.*(-j0rn1+j0rn2).*jsr0.*r+as.*j0rn0.*(-jsr1+jsr2).*rn).*wc)+(4.*asx.*c0.*e0yc.*j0rn0.*jsr0+1i.*(as.*c0.*ej3.*(jsr1-jsr2)+8.*e0zc.*j0rn0.*jsr0.*w0)+4.*asy.*j0rn0.*jsr0.*(-c0.*e0xc+a0.*e0zc.*wc)).*x);
        YBmz1 = (-1i).*c0.*ej3.*jsr0.*(-wsr).*ws.^(-1).*x.^3;
        YBmz2 = (-2i).*c0.*e0zc.*j0rn0.*jsr0.*(-wsr).*ws.^(-1).*x.^2;
        YBmz3 = (-2i).*e0zc.*j0rn0.*jsr0.*(-wsr).*ws.^(-1).*x.^2;
        YBmz4 = (1i*(-1/2)).*ws.^(-1).*x.*(a0.*ej8.*jsr0.*r.*wc.*x+as.*ej3.*(-jsr1+jsr2).*rn.*wc.*x+(-2).*ej3.*jsr0.*(r.*wc+(1i.*a0.*asy.*wc+(-2).*ws).*x.^2));
        YBmz5 = 1i.*e0zc.*ws.^(-1).*x.*(as.*j0rn0.*(jsr1-jsr2).*rn.*wc+(-4).*j0rn0.*jsr0.*ws.*x+a0.*jsr0.*wc.*(j0rn1.*r-j0rn2.*r+(2i).*asy.*j0rn0.*x));
        
        YAzx1 = (1i*(-1/2)).*cs.*ej1.*w0.^(-1).*x.^2.*(2.*e0zc.*j0rn0.*w0rn+c0.*ej3.*x);
        YAzx2 = (1/4).*w0.^(-1).*x.*(1i.*as.*ej3.*ej6.*w0.*x+2.*ej1.*((-1i).*(ej3+ej4.*r).*w0+2.*j0rn0.*(asy.*e0xc.*w0-e0yc.*(asx.*w0+a0.*r.*wc)).*x+((2i).*ej3.*w0+a0.*e0yc.*(1i.*asy.*(j0rn1-j0rn2)+asx.*(j0rn1+j0rn2)).*wc).*x.^2));
        YAzx3 = (1/4).*w0.^(-1).*((-4i).*e0zc.*ej1.*j0rn0.*rn.*wc+(2i).*(c0.*ej1.*(ej3+ej4.*r)+e0zc.*(a0.*ej1.*(-j0rn1+j0rn2).*r+as.*ej6.*j0rn0.*rn).*wc).*x+((-1i).*as.*c0.*ej3.*ej6+4.*ej1.*j0rn0.*(-asy.*c0.*e0xc+asx.*c0.*e0yc+(2i).*e0zc.*w0+a0.*asy.*e0zc.*wc)).*x.^2);
        YBzx1 = (1i*(1/2)).*c0.*cs.*ej1.*ej3.*ws.^(-1).*x.^3;
        YBzx2 = 1i.*c0.*cs.*e0zc.*ej1.*j0rn0.*ws.^(-1).*x.^2;
        YBzx3 = (1/4).*ws.^(-1).*x.*((2i).*ej3.*(ej1+ej2.*rn).*ws+((4i).*cs.*e0zc.*ej1.*j0rn0+4.*asy.*ej3.*jsr0.*rn.*wc+(-1i).*a0.*ej1.*ej8.*ws).*x+(-2).*ej1.*ej3.*(a0.*asy.*wc+(2i).*ws).*x.^2);
        YBzx4 = (1i*(1/4)).*ws.^(-1).*(cs.*x.*((-2).*ej3.*(ej1+ej2.*rn)+a0.*ej1.*ej8.*x)+2.*e0zc.*(2.*ej2.*j0rn0.*rn.*ws+(2i).*asy.*j0rn0.*wc.*x.*((-2).*jsr0.*rn+a0.*ej1.*x)+ej1.*ws.*x.*(a0.*(j0rn1-j0rn2)+(-4).*j0rn0.*x)));
        YBzx5 = (1i*(-1/2)).*cs.*e0zc.*ws.^(-1).*(2.*ej2.*j0rn0.*rn+a0.*ej1.*(j0rn1-j0rn2).*x);
        
        YAzy1 = (-1/2).*cs.*ej2.*w0.^(-1).*x.^2.*(2.*e0zc.*j0rn0.*w0rn+c0.*ej3.*x);
        YAzy2 = (1/4).*w0.^(-1).*x.*(as.*ej3.*ej5.*w0.*x+2.*ej2.*(-(ej3+ej4.*r).*w0+(2i).*j0rn0.*(-asy.*e0xc.*w0+asx.*e0yc.*w0+a0.*e0yc.*r.*wc).*x+(2.*ej3.*w0+a0.*e0yc.*(asy.*(j0rn1-j0rn2)+(-1i).*asx.*(j0rn1+j0rn2)).*wc).*x.^2));
        YAzy3 = (-1/4).*w0.^(-1).*(4.*e0zc.*ej2.*j0rn0.*rn.*wc+(-2).*(c0.*ej2.*(ej3+ej4.*r)+e0zc.*(a0.*ej2.*(-j0rn1+j0rn2).*r+as.*ej5.*j0rn0.*rn).*wc).*x+(as.*c0.*ej3.*ej5+(4i).*ej2.*j0rn0.*(-asy.*c0.*e0xc+asx.*c0.*e0yc+(2i).*e0zc.*w0+a0.*asy.*e0zc.*wc)).*x.^2);
        YBzy1 = (1/2).*c0.*cs.*ej2.*ej3.*ws.^(-1).*x.^3;
        YBzy2 = c0.*cs.*e0zc.*ej2.*j0rn0.*ws.^(-1).*x.^2;
        YBzy3 = (1/4).*ws.^(-1).*x.*(2.*ej3.*(ej2+ej1.*rn).*ws+(4.*cs.*e0zc.*ej2.*j0rn0-a0.*ej2.*ej8.*ws+(-4).*ej3.*jsr0.*(asx.*rn.*wc+a0.*ws)).*x+2.*ej3.*(a0.*asx.*ej1.*wc+(-2).*ej2.*ws).*x.^2);
        YBzy4 = (1/4).*ws.^(-1).*(4.*e0zc.*ej1.*j0rn0.*rn.*ws+(-2).*(cs.*ej3.*(ej2+ej1.*rn)+4.*asx.*e0zc.*j0rn0.*jsr0.*rn.*wc+a0.*e0zc.*(-ej2.*j0rn1+ej2.*j0rn2+4.*j0rn0.*jsr0).*ws).*x+(a0.*(cs.*ej2.*ej8+4.*cs.*ej3.*jsr0+4.*asx.*e0zc.*ej1.*j0rn0.*wc)+(-8).*e0zc.*ej2.*j0rn0.*ws).*x.^2);
        YBzy5 = (1/2).*cs.*e0zc.*ws.^(-1).*((-2).*ej1.*j0rn0.*rn+a0.*(ej2.*(-j0rn1+j0rn2)+4.*j0rn0.*jsr0).*x);
        
        YAzz1 = (-1i).*cs.*jsr0.*w0.^(-1).*x.*(2.*e0zc.*j0rn0.*w0rn+c0.*ej3.*x);
        YAzz2 = (-1i).*jsr0.*w0.^(-1).*x.*(2.*e0zc.*j0rn0.*w0rn+c0.*ej3.*x);
        YAzz3 = (1/2).*w0.^(-1).*((-2i).*ej4.*jsr0.*r.*w0+(4.*asy.*e0xc.*j0rn0.*jsr0.*w0+(-4).*asx.*e0yc.*j0rn0.*jsr0.*w0+(-1i).*as.*ej3.*(jsr1-jsr2).*w0+(-4).*a0.*e0yc.*j0rn0.*jsr0.*r.*wc).*x+2.*jsr0.*((2i).*ej3.*w0+a0.*e0yc.*(1i.*asy.*(j0rn1-j0rn2)+asx.*(j0rn1+j0rn2)).*wc).*x.^2);
        YAzz4 = (1/2).*w0.^(-1).*((2i).*(c0.*ej4.*jsr0.*r+e0zc.*(a0.*(-j0rn1+j0rn2).*jsr0.*r+as.*j0rn0.*(-jsr1+jsr2).*rn).*wc)+(4.*asx.*c0.*e0yc.*j0rn0.*jsr0+1i.*(as.*c0.*ej3.*(jsr1-jsr2)+8.*e0zc.*j0rn0.*jsr0.*w0)+4.*asy.*j0rn0.*jsr0.*(-c0.*e0xc+a0.*e0zc.*wc)).*x);
        YBzz1 = (-1i).*c0.*ej3.*jsr0.*(-wsr).*ws.^(-1).*x.^2;
        YBzz2 = (-2i).*c0.*e0zc.*j0rn0.*jsr0.*(-wsr).*ws.^(-1).*x;
        YBzz3 = (-2i).*e0zc.*j0rn0.*jsr0.*(-wsr).*ws.^(-1).*x;
        YBzz4 = (1i*(1/2)).*ws.^(-1).*(-a0.*ej8.*jsr0.*r.*wc.*x+as.*ej3.*(jsr1-jsr2).*rn.*wc.*x+2.*ej3.*jsr0.*(r.*wc+(1i.*a0.*asy.*wc+(-2).*ws).*x.^2));
        YBzz5 = 1i.*e0zc.*ws.^(-1).*(as.*j0rn0.*(jsr1-jsr2).*rn.*wc+(-4).*j0rn0.*jsr0.*ws.*x+a0.*jsr0.*wc.*(j0rn1.*r-j0rn2.*r+(2i).*asy.*j0rn0.*x));
        
        cmx = exp(-x.^2).*(YAmx1*ZAmx1+YAmx2*ZAmx2+YAmx3*ZAmx3+YBmx1*ZBmx1+YBmx2*ZBmx2+YBmx3*ZBmx3+YBmx4*ZBmx4+YBmx5*ZBmx5);
        cmy = exp(-x.^2).*(YAmy1*ZAmy1+YAmy2*ZAmy2+YAmy3*ZAmy3+YBmy1*ZBmy1+YBmy2*ZBmy2+YBmy3*ZBmy3+YBmy4*ZBmy4+YBmy5*ZBmy5);
        cmz = exp(-x.^2).*(YAmz1*ZAmz1+YAmz2*ZAmz2+YAmz3*ZAmz3+YAmz4*ZAmz4+YBmz1*ZBmz1+YBmz2*ZBmz2+YBmz3*ZBmz3+YBmz4*ZBmz4+YBmz5*ZBmz5);
        czx = exp(-x.^2).*(YAzx1*ZAzx1+YAzx2*ZAzx2+YAzx3*ZAzx3+YBzx1*ZBzx1+YBzx2*ZBzx2+YBzx3*ZBzx3+YBzx4*ZBzx4+YBzx5*ZBzx5);
        czy = exp(-x.^2).*(YAzy1*ZAzy1+YAzy2*ZAzy2+YAzy3*ZAzy3+YBzy1*ZBzy1+YBzy2*ZBzy2+YBzy3*ZBzy3+YBzy4*ZBzy4+YBzy5*ZBzy5);
        czz = exp(-x.^2).*(YAzz1*ZAzz1+YAzz2*ZAzz2+YAzz3*ZAzz3+YAzz4*ZAzz4+YBzz1*ZBzz1+YBzz2*ZBzz2+YBzz3*ZBzz3+YBzz4*ZBzz4+YBzz5*ZBzz5);
        
        cxx = cmx.*(Ei1c*j1n1+Ei1*j1n2)/2;
        cyx = cmx.*1i.*(Ei1c*j1n1-Ei1*j1n2)/2;
        cxy = cmy.*(Ei1c*j1n1+Ei1*j1n2)/2;
        cyy = cmy.*1i.*(Ei1c*j1n1-Ei1*j1n2)/2;
        cxz = cmz.*(Ei1c*j1n1+Ei1*j1n2)/2;
        cyz = cmz.*1i.*(Ei1c*j1n1-Ei1*j1n2)/2;
        czx = czx*j1n0;
        czy = czy*j1n0;
        czz = czz*j1n0;
        
        y = [cxx,cxy,cxz;cyx,cyy,cyz;czx,czy,czz];
    end

mat = integral(@(x) fxyz(x),0, inf,'ArrayValued',true);
res = mat*beta0*wp^2*(w0/w1)*exp(1i*de);

end

