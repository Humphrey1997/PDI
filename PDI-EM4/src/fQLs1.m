function res = fQLs1(ws,n,r)
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
de = n*des-r*de1;
nr = n-r;

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
wsn = ws-n*wc;
w1r = w1-r*wc;
w0nr = w0-nr*wc;
zetasn = wsn/cs;
zeta1r = w1r/c1;
zeta0nr = w0nr/c0;

ZA012=g2(0,zetasn,1,zeta1r,2)/cs^1/c1^2;
ZA011=g2(0,zetasn,1,zeta1r,1)/cs^1/c1^1;
ZA112=g2(1,zetasn,1,zeta1r,2)/cs^1/c1^2;
ZA212=g2(2,zetasn,1,zeta1r,2)/cs^1/c1^2;
ZA111=g2(1,zetasn,1,zeta1r,1)/cs^1/c1^1;
ZA211=g2(2,zetasn,1,zeta1r,1)/cs^1/c1^1;
ZA311=g2(3,zetasn,1,zeta1r,1)/cs^1/c1^1;

ZAmx1=ZA012;ZAmx2=ZA011;ZAmx3=ZA111;
ZAmy1=ZA012;ZAmy2=ZA011;ZAmy3=ZA111;
ZAmz1=ZA112;ZAmz2=ZA011;ZAmz3=ZA111;ZAmz4=ZA211;
ZAzx1=ZA112;ZAzx2=ZA111;ZAzx3=ZA211;
ZAzy1=ZA112;ZAzy2=ZA111;ZAzy3=ZA211;
ZAzz1=ZA212;ZAzz2=ZA111;ZAzz3=ZA211;ZAzz4=ZA311;

ZB012=g2(0,zetasn,1,zeta0nr,2)/cs^1/c0^2;
ZB011=g2(0,zetasn,1,zeta0nr,1)/cs^1/c0^1;
ZB112=g2(1,zetasn,1,zeta0nr,2)/cs^1/c0^2;
ZB212=g2(2,zetasn,1,zeta0nr,2)/cs^1/c0^2;
ZB111=g2(1,zetasn,1,zeta0nr,1)/cs^1/c0^1;
ZB211=g2(2,zetasn,1,zeta0nr,1)/cs^1/c0^1;
ZB311=g2(3,zetasn,1,zeta0nr,1)/cs^1/c0^1;

ZBmx1=ZB012;ZBmx2=ZB112;ZBmx3=ZB011;ZBmx4=ZB111;ZBmx5=ZB211;
ZBmy1=ZB012;ZBmy2=ZB112;ZBmy3=ZB011;ZBmy4=ZB111;ZBmy5=ZB211;
ZBmz1=ZB012;ZBmz2=ZB112;ZBmz3=ZB011;ZBmz4=ZB111;ZBmz5=ZB211;
ZBzx1=ZB112;ZBzx2=ZB212;ZBzx3=ZB111;ZBzx4=ZB211;ZBzx5=ZB311;
ZBzy1=ZB112;ZBzy2=ZB212;ZBzy3=ZB111;ZBzy4=ZB211;ZBzy5=ZB311;
ZBzz1=ZB112;ZBzz2=ZB212;ZBzz3=ZB111;ZBzz4=ZB211;ZBzz5=ZB311;

    function y = fxyz(x)
        j1r0 = besselj(r,a1*x);
        j1r1 = besselj(r-1,a1*x);
        j1r2 = besselj(r+1,a1*x);
        j1r3 = besselj(r-2,a1*x);
        j1r4 = besselj(r+2,a1*x);
        
        j0nr0 = besselj(nr,a0*x);
        j0nr1 = besselj(nr-1,a0*x);
        j0nr2 = besselj(nr+1,a0*x);
        j0nr3 = besselj(nr-2,a0*x);
        j0nr4 = besselj(nr+2,a0*x);
        
        jsn0 = besselj(n,as*x);
        jsn1 = besselj(n-1,as*x);
        jsn2 = besselj(n+1,as*x);
        
        ej1 = Ei1*j1r1+Ei1c*j1r2;
        ej2 = Ei1*j1r1-Ei1c*j1r2;
        ej3 = e0s*j0nr1+e0p*j0nr2;
        ej4 = e0s*j0nr1-e0p*j0nr2;
        ej5 = Ei1*(j1r0-j1r3)+Ei1c*(j1r0-j1r4);
        ej6 = Ei1*(j1r0-j1r3)-Ei1c*(j1r0-j1r4);
        ej7 = e0s*(j0nr0-j0nr3)+e0p*(j0nr0-j0nr4);
        ej8 = e0s*(j0nr0-j0nr3)-e0p*(j0nr0-j0nr4);
        
        YAmx1 = (1i*(-1/2)).*c1.*ej1.*w0.^(-1).*x.^3.*(2.*e0z.*j0nr0.*w0nr+c0.*ej3.*x);
        YAmx2 = (1/4).*w0.^(-1).*x.^2.*(1i.*a1.*ej3.*ej6.*w0.*x+2.*ej1.*(1i.*ej4.*r.*w0+(-2).*j0nr0.*(-a1y.*e0x.*w0+a1x.*e0y.*w0+a0.*e0y.*r.*wc).*x+a0.*e0y.*((-1i).*a1y.*(j0nr1-j0nr2)+a1x.*(j0nr1+j0nr2)).*wc.*x.^2+1i.*ej3.*w0.*(-1+2.*x.^2)));
        YAmx3 = (1/4).*w0.^(-1).*x.*((-4i).*e0z.*ej1.*j0nr0.*nr.*wc+(2i).*(c0.*ej1.*(ej3-ej4.*r)+e0z.*(a1.*ej6.*j0nr0.*nr+a0.*ej1.*(j0nr1-j0nr2).*r).*wc).*x+((-1i).*a1.*c0.*ej3.*ej6+4.*ej1.*j0nr0.*(-a1y.*c0.*e0x+a1x.*c0.*e0y+(2i).*e0z.*w0+a0.*a1y.*e0z.*wc)).*x.^2);
        YBmx1 = (1i*(-1/2)).*c0.*c1.*ej1.*ej3.*w1.^(-1).*x.^4;
        YBmx2 = (-1i).*c0.*c1.*e0z.*ej1.*j0nr0.*w1.^(-1).*x.^3;
        YBmx3 = (1i*(1/4)).*w1.^(-1).*x.^2.*(2.*ej3.*nr.*(ej2.*w1+(-2i).*a1y.*j1r0.*wc.*x)+ej1.*((-2).*ej3.*w1+(-4).*c1.*e0z.*j0nr0.*x+a0.*ej8.*w1.*x+2.*ej3.*(2.*w1+1i.*a0.*a1y.*wc).*x.^2));
        YBmx4 = (1i*(1/4)).*w1.^(-1).*x.*(4.*e0z.*ej2.*j0nr0.*nr.*w1+2.*(c1.*ej3.*(ej1-ej2.*nr)+a0.*e0z.*ej1.*(-j0nr1+j0nr2).*w1+(-4i).*a1y.*e0z.*j0nr0.*j1r0.*nr.*wc).*x+ej1.*(-a0.*c1.*ej8+8.*e0z.*j0nr0.*w1+(4i).*a0.*a1y.*e0z.*j0nr0.*wc).*x.^2);
        YBmx5 = (1i*(-1/2)).*c1.*e0z.*w1.^(-1).*x.*(2.*ej2.*j0nr0.*nr+a0.*ej1.*(-j0nr1+j0nr2).*x);
        
        YAmy1 = (-1/2).*c1.*ej2.*w0.^(-1).*x.^3.*(2.*e0z.*j0nr0.*w0nr+c0.*ej3.*x);
        YAmy2 = (1/4).*w0.^(-1).*x.^2.*(a1.*ej3.*ej5.*w0.*x+2.*ej2.*(ej4.*r.*w0+(2i).*j0nr0.*(-a1y.*e0x.*w0+a1x.*e0y.*w0+a0.*e0y.*r.*wc).*x+a0.*e0y.*(a1y.*(-j0nr1+j0nr2)+(-1i).*a1x.*(j0nr1+j0nr2)).*wc.*x.^2+ej3.*w0.*(-1+2.*x.^2)));
        YAmy3 = (1/4).*w0.^(-1).*x.*((-4).*e0z.*ej2.*j0nr0.*nr.*wc+2.*(c0.*ej2.*(ej3-ej4.*r)+e0z.*(a1.*ej5.*j0nr0.*nr+a0.*ej2.*(j0nr1-j0nr2).*r).*wc).*x-(a1.*c0.*ej3.*ej5+(4i).*ej2.*j0nr0.*(-a1y.*c0.*e0x+a1x.*c0.*e0y+(2i).*e0z.*w0+a0.*a1y.*e0z.*wc)).*x.^2);
        YBmy1 = (-1/2).*c0.*c1.*ej2.*ej3.*w1.^(-1).*x.^4;
        YBmy2 = (-1).*c0.*c1.*e0z.*ej2.*j0nr0.*w1.^(-1).*x.^3;
        YBmy3 = (1/4).*w1.^(-1).*x.^2.*(2.*ej3.*(ej1.*nr.*w1+(-2).*j1r0.*(a0.*w1+a1x.*nr.*wc).*x+a0.*a1x.*ej1.*wc.*x.^2)+ej2.*((-4).*c1.*e0z.*j0nr0.*x+a0.*ej8.*w1.*x+2.*ej3.*w1.*(-1+2.*x.^2)));
        YBmy4 = (1/4).*w1.^(-1).*x.*(c1.*x.*(2.*ej2.*ej3+(-2).*ej1.*ej3.*nr-a0.*ej2.*ej8.*x+4.*a0.*ej3.*j1r0.*x)+2.*e0z.*x.*(a0.*(ej2.*(-j0nr1+j0nr2)+(-4).*j0nr0.*j1r0).*w1+(-4).*a1x.*j0nr0.*j1r0.*nr.*wc+4.*ej2.*j0nr0.*w1.*x)+4.*e0z.*ej1.*j0nr0.*(nr.*w1+a0.*a1x.*wc.*x.^2));
        YBmy5 = (1/2).*c1.*e0z.*w1.^(-1).*x.*((-2).*ej1.*j0nr0.*nr+a0.*(ej2.*(j0nr1-j0nr2)+4.*j0nr0.*j1r0).*x);
        
        YAmz1 = (-1i).*c1.*j1r0.*w0.^(-1).*x.^2.*(2.*e0z.*j0nr0.*w0nr+c0.*ej3.*x);
        YAmz2 = (-1i).*j1r0.*w0.^(-1).*x.^2.*(2.*e0z.*j0nr0.*w0nr+c0.*ej3.*x);
        YAmz3 = (1/2).*w0.^(-1).*x.*((2i).*ej4.*j1r0.*r.*w0+x.*((-4).*a0.*e0y.*j0nr0.*j1r0.*r.*wc+2.*a0.*e0y.*((-1i).*a1y.*(j0nr1-j0nr2)+a1x.*(j0nr1+j0nr2)).*j1r0.*wc.*x+w0.*(4.*a1y.*e0x.*j0nr0.*j1r0+(-4).*a1x.*e0y.*j0nr0.*j1r0+(-1i).*a1.*ej3.*(j1r1-j1r2)+(4i).*ej3.*j1r0.*x)));
        YAmz4 = (1/2).*w0.^(-1).*x.*((2i).*e0z.*(a1.*j0nr0.*(-j1r1+j1r2).*nr+a0.*(j0nr1-j0nr2).*j1r0.*r).*wc+4.*e0z.*j0nr0.*j1r0.*((2i).*w0+a0.*a1y.*wc).*x+c0.*((-2i).*ej4.*j1r0.*r+4.*(-a1y.*e0x+a1x.*e0y).*j0nr0.*j1r0.*x+1i.*a1.*ej3.*(j1r1-j1r2).*x));
        YBmz1 = (-1i).*c0.*ej3.*j1r0.*w1.^(-1).*w1r.*x.^3;
        YBmz2 = (-2i).*c0.*e0z.*j0nr0.*j1r0.*w1.^(-1).*w1r.*x.^2;
        YBmz3 = (-2i).*e0z.*j0nr0.*j1r0.*w1.^(-1).*w1r.*x.^2;
        YBmz4 = (1i*(1/2)).*w1.^(-1).*x.*(a0.*ej8.*j1r0.*r.*wc.*x+ej3.*((-2).*j1r0.*r.*wc+a1.*(j1r1-j1r2).*nr.*wc.*x+2.*j1r0.*(2.*w1+1i.*a0.*a1y.*wc).*x.^2));
        YBmz5 = 1i.*e0z.*w1.^(-1).*x.*(a1.*j0nr0.*(j1r1-j1r2).*nr.*wc+4.*j0nr0.*j1r0.*w1.*x+a0.*j1r0.*wc.*(-j0nr1.*r+j0nr2.*r+(2i).*a1y.*j0nr0.*x));
        
        YAzx1 = (1i*(-1/2)).*c1.*ej1.*w0.^(-1).*x.^2.*(2.*e0z.*j0nr0.*w0nr+c0.*ej3.*x);
        YAzx2 = (1/4).*w0.^(-1).*x.*(1i.*a1.*ej3.*ej6.*w0.*x+2.*ej1.*(1i.*ej4.*r.*w0+(-2).*j0nr0.*(-a1y.*e0x.*w0+a1x.*e0y.*w0+a0.*e0y.*r.*wc).*x+a0.*e0y.*((-1i).*a1y.*(j0nr1-j0nr2)+a1x.*(j0nr1+j0nr2)).*wc.*x.^2+1i.*ej3.*w0.*(-1+2.*x.^2)));
        YAzx3 = (1/4).*w0.^(-1).*((-4i).*e0z.*ej1.*j0nr0.*nr.*wc+(2i).*(c0.*ej1.*(ej3-ej4.*r)+e0z.*(a1.*ej6.*j0nr0.*nr+a0.*ej1.*(j0nr1-j0nr2).*r).*wc).*x+((-1i).*a1.*c0.*ej3.*ej6+4.*ej1.*j0nr0.*(-a1y.*c0.*e0x+a1x.*c0.*e0y+(2i).*e0z.*w0+a0.*a1y.*e0z.*wc)).*x.^2);
        YBzx1 = (1i*(-1/2)).*c0.*c1.*ej1.*ej3.*w1.^(-1).*x.^3;
        YBzx2 = (-1i).*c0.*c1.*e0z.*ej1.*j0nr0.*w1.^(-1).*x.^2;
        YBzx3 = (1i*(1/4)).*w1.^(-1).*x.*(2.*ej3.*nr.*(ej2.*w1+(-2i).*a1y.*j1r0.*wc.*x)+ej1.*((-2).*ej3.*w1+(-4).*c1.*e0z.*j0nr0.*x+a0.*ej8.*w1.*x+2.*ej3.*(2.*w1+1i.*a0.*a1y.*wc).*x.^2));
        YBzx4 = (1i*(1/4)).*w1.^(-1).*(-c1.*x.*((-2).*ej1.*ej3+2.*ej2.*ej3.*nr+a0.*ej1.*ej8.*x)+2.*e0z.*(2.*ej2.*j0nr0.*nr.*w1+4.*j0nr0.*x.*((-1i).*a1y.*j1r0.*nr.*wc+ej1.*w1.*x)+a0.*ej1.*x.*(-j0nr1.*w1+j0nr2.*w1+(2i).*a1y.*j0nr0.*wc.*x)));
        YBzx5 = (1i*(-1/2)).*c1.*e0z.*w1.^(-1).*(2.*ej2.*j0nr0.*nr+a0.*ej1.*(-j0nr1+j0nr2).*x);
        
        YAzy1 = (-1/2).*c1.*ej2.*w0.^(-1).*x.^2.*(2.*e0z.*j0nr0.*w0nr+c0.*ej3.*x);
        YAzy2 = (1/4).*w0.^(-1).*x.*(a1.*ej3.*ej5.*w0.*x+2.*ej2.*(ej4.*r.*w0+(2i).*j0nr0.*(-a1y.*e0x.*w0+a1x.*e0y.*w0+a0.*e0y.*r.*wc).*x+a0.*e0y.*(a1y.*(-j0nr1+j0nr2)+(-1i).*a1x.*(j0nr1+j0nr2)).*wc.*x.^2+ej3.*w0.*(-1+2.*x.^2)));
        YAzy3 = (-1/4).*w0.^(-1).*(4.*e0z.*ej2.*j0nr0.*nr.*wc+(-2).*(c0.*ej2.*(ej3-ej4.*r)+e0z.*(a1.*ej5.*j0nr0.*nr+a0.*ej2.*(j0nr1-j0nr2).*r).*wc).*x+(a1.*c0.*ej3.*ej5+(4i).*ej2.*j0nr0.*(-a1y.*c0.*e0x+a1x.*c0.*e0y+(2i).*e0z.*w0+a0.*a1y.*e0z.*wc)).*x.^2);
        YBzy1 = (-1/2).*c0.*c1.*ej2.*ej3.*w1.^(-1).*x.^3;
        YBzy2 = (-1).*c0.*c1.*e0z.*ej2.*j0nr0.*w1.^(-1).*x.^2;
        YBzy3 = (1/4).*w1.^(-1).*x.*(2.*ej3.*(ej1.*nr.*w1+(-2).*j1r0.*(a0.*w1+a1x.*nr.*wc).*x+a0.*a1x.*ej1.*wc.*x.^2)+ej2.*((-4).*c1.*e0z.*j0nr0.*x+a0.*ej8.*w1.*x+2.*ej3.*w1.*(-1+2.*x.^2)));
        YBzy4 = (1/4).*w1.^(-1).*(c1.*x.*(2.*ej2.*ej3+(-2).*ej1.*ej3.*nr-a0.*ej2.*ej8.*x+4.*a0.*ej3.*j1r0.*x)+2.*e0z.*x.*(a0.*(ej2.*(-j0nr1+j0nr2)+(-4).*j0nr0.*j1r0).*w1+(-4).*a1x.*j0nr0.*j1r0.*nr.*wc+4.*ej2.*j0nr0.*w1.*x)+4.*e0z.*ej1.*j0nr0.*(nr.*w1+a0.*a1x.*wc.*x.^2));
        YBzy5 = (1/2).*c1.*e0z.*w1.^(-1).*((-2).*ej1.*j0nr0.*nr+a0.*(ej2.*(j0nr1-j0nr2)+4.*j0nr0.*j1r0).*x);
        
        YAzz1 = (-1i).*c1.*j1r0.*w0.^(-1).*x.*(2.*e0z.*j0nr0.*w0nr+c0.*ej3.*x);
        YAzz2 = (-1i).*j1r0.*w0.^(-1).*x.*(2.*e0z.*j0nr0.*w0nr+c0.*ej3.*x);
        YAzz3 = 1i.*ej4.*j1r0.*r+(1/2).*w0.^(-1).*x.*((-4).*a0.*e0y.*j0nr0.*j1r0.*r.*wc+2.*a0.*e0y.*((-1i).*a1y.*(j0nr1-j0nr2)+a1x.*(j0nr1+j0nr2)).*j1r0.*wc.*x+w0.*(4.*a1y.*e0x.*j0nr0.*j1r0+(-4).*a1x.*e0y.*j0nr0.*j1r0+(-1i).*a1.*ej3.*(j1r1-j1r2)+(4i).*ej3.*j1r0.*x));
        YAzz4 = (1/2).*w0.^(-1).*((2i).*e0z.*(a1.*j0nr0.*(-j1r1+j1r2).*nr+a0.*(j0nr1-j0nr2).*j1r0.*r).*wc+4.*e0z.*j0nr0.*j1r0.*((2i).*w0+a0.*a1y.*wc).*x+c0.*((-2i).*ej4.*j1r0.*r+4.*(-a1y.*e0x+a1x.*e0y).*j0nr0.*j1r0.*x+1i.*a1.*ej3.*(j1r1-j1r2).*x));
        YBzz1 = (-1i).*c0.*ej3.*j1r0.*w1.^(-1).*w1r.*x.^2;
        YBzz2 = (-2i).*c0.*e0z.*j0nr0.*j1r0.*w1.^(-1).*w1r.*x;
        YBzz3 = (-2i).*e0z.*j0nr0.*j1r0.*w1.^(-1).*w1r.*x;
        YBzz4 = (1i*(1/2)).*w1.^(-1).*(a0.*ej8.*j1r0.*r.*wc.*x+ej3.*((-2).*j1r0.*r.*wc+a1.*(j1r1-j1r2).*nr.*wc.*x+2.*j1r0.*(2.*w1+1i.*a0.*a1y.*wc).*x.^2));
        YBzz5 = 1i.*e0z.*w1.^(-1).*(a1.*j0nr0.*(j1r1-j1r2).*nr.*wc+4.*j0nr0.*j1r0.*w1.*x+a0.*j1r0.*wc.*(-j0nr1.*r+j0nr2.*r+(2i).*a1y.*j0nr0.*x));
        
        cmx = exp(-x.^2).*(YAmx1*ZAmx1+YAmx2*ZAmx2+YAmx3*ZAmx3+YBmx1*ZBmx1+YBmx2*ZBmx2+YBmx3*ZBmx3+YBmx4*ZBmx4+YBmx5*ZBmx5);
        cmy = exp(-x.^2).*(YAmy1*ZAmy1+YAmy2*ZAmy2+YAmy3*ZAmy3+YBmy1*ZBmy1+YBmy2*ZBmy2+YBmy3*ZBmy3+YBmy4*ZBmy4+YBmy5*ZBmy5);
        cmz = exp(-x.^2).*(YAmz1*ZAmz1+YAmz2*ZAmz2+YAmz3*ZAmz3+YAmz4*ZAmz4+YBmz1*ZBmz1+YBmz2*ZBmz2+YBmz3*ZBmz3+YBmz4*ZBmz4+YBmz5*ZBmz5);
        czx = exp(-x.^2).*(YAzx1*ZAzx1+YAzx2*ZAzx2+YAzx3*ZAzx3+YBzx1*ZBzx1+YBzx2*ZBzx2+YBzx3*ZBzx3+YBzx4*ZBzx4+YBzx5*ZBzx5);
        czy = exp(-x.^2).*(YAzy1*ZAzy1+YAzy2*ZAzy2+YAzy3*ZAzy3+YBzy1*ZBzy1+YBzy2*ZBzy2+YBzy3*ZBzy3+YBzy4*ZBzy4+YBzy5*ZBzy5);
        czz = exp(-x.^2).*(YAzz1*ZAzz1+YAzz2*ZAzz2+YAzz3*ZAzz3+YAzz4*ZAzz4+YBzz1*ZBzz1+YBzz2*ZBzz2+YBzz3*ZBzz3+YBzz4*ZBzz4+YBzz5*ZBzz5);
        
        cxx = cmx.*(Eisc*jsn1+Eis*jsn2)/2;
        cyx = cmx.*1i.*(Eisc*jsn1-Eis*jsn2)/2;
        cxy = cmy.*(Eisc*jsn1+Eis*jsn2)/2;
        cyy = cmy.*1i.*(Eisc*jsn1-Eis*jsn2)/2;
        cxz = cmz.*(Eisc*jsn1+Eis*jsn2)/2;
        cyz = cmz.*1i.*(Eisc*jsn1-Eis*jsn2)/2;
        czx = czx*jsn0;
        czy = czy*jsn0;
        czz = czz*jsn0;
        
        y = [cxx,cxy,cxz;cyx,cyy,cyz;czx,czy,czz];
    end

mat = integral(@(x) fxyz(x),0, inf,'ArrayValued',true);
res = mat*beta0*wp^2*(w0/ws)*exp(1i*de);

end

