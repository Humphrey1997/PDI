function res = fQLs2(ws,n,r)
global wpe wce ve;
global nsp n2p n0p;
global nsz n2z n0z;
global e0x e0y e0z;
global des de2 de0;
global beta0 w0;

wc = wce;
vt = ve;
wp = wpe;

w2 = ws+w0;
de = n*des-r*de2;
rn = r-n;

e0xc = conj(e0x);
e0yc = conj(e0y);
e0zc = conj(e0z);
e0p = e0x + 1i*e0y;
e0s = e0x - 1i*e0y;
e0pc = conj(e0p);
e0sc = conj(e0s);

as = nsp*vt/wc;
a2 = n2p*vt/wc;
a0 = n0p*vt/wc;
asx = as*cos(des);
asy = as*sin(des);
a2x = a2*cos(de2);
a2y = a2*sin(de2);

Eis = exp(1i*des);
Eisc = conj(Eis);
Ei2 = exp(1i*de2);
Ei2c = conj(Ei2);

cs = nsz*vt;
c2 = n2z*vt;
c0 = n0z*vt;
wsn = ws-n*wc;
w2r = w2-r*wc;
w0rn = w0-rn*wc;
zetasn = wsn/cs;
zeta2r = w2r/c2;
zeta0rn = w0rn/c0;

ZA012=g2(0,zetasn,1,zeta2r,2)/cs^1/c2^2;
ZA011=g2(0,zetasn,1,zeta2r,1)/cs^1/c2^1;
ZA112=g2(1,zetasn,1,zeta2r,2)/cs^1/c2^2;
ZA212=g2(2,zetasn,1,zeta2r,2)/cs^1/c2^2;
ZA111=g2(1,zetasn,1,zeta2r,1)/cs^1/c2^1;
ZA211=g2(2,zetasn,1,zeta2r,1)/cs^1/c2^1;
ZA311=g2(3,zetasn,1,zeta2r,1)/cs^1/c2^1;

ZAmx1=ZA012;ZAmx2=ZA011;ZAmx3=ZA111;
ZAmy1=ZA012;ZAmy2=ZA011;ZAmy3=ZA111;
ZAmz1=ZA112;ZAmz2=ZA011;ZAmz3=ZA111;ZAmz4=ZA211;
ZAzx1=ZA112;ZAzx2=ZA111;ZAzx3=ZA211;
ZAzy1=ZA112;ZAzy2=ZA111;ZAzy3=ZA211;
ZAzz1=ZA212;ZAzz2=ZA111;ZAzz3=ZA211;ZAzz4=ZA311;

ZB012=g2(0,zetasn,1,zeta0rn,2)/cs^1/c0^2;
ZB011=g2(0,zetasn,1,zeta0rn,1)/cs^1/c0^1;
ZB112=g2(1,zetasn,1,zeta0rn,2)/cs^1/c0^2;
ZB212=g2(2,zetasn,1,zeta0rn,2)/cs^1/c0^2;
ZB111=g2(1,zetasn,1,zeta0rn,1)/cs^1/c0^1;
ZB211=g2(2,zetasn,1,zeta0rn,1)/cs^1/c0^1;
ZB311=g2(3,zetasn,1,zeta0rn,1)/cs^1/c0^1;

ZBmx1=ZB012;ZBmx2=ZB112;ZBmx3=ZB011;ZBmx4=ZB111;ZBmx5=ZB211;
ZBmy1=ZB012;ZBmy2=ZB112;ZBmy3=ZB011;ZBmy4=ZB111;ZBmy5=ZB211;
ZBmz1=ZB012;ZBmz2=ZB112;ZBmz3=ZB011;ZBmz4=ZB111;ZBmz5=ZB211;
ZBzx1=ZB112;ZBzx2=ZB212;ZBzx3=ZB111;ZBzx4=ZB211;ZBzx5=ZB311;
ZBzy1=ZB112;ZBzy2=ZB212;ZBzy3=ZB111;ZBzy4=ZB211;ZBzy5=ZB311;
ZBzz1=ZB112;ZBzz2=ZB212;ZBzz3=ZB111;ZBzz4=ZB211;ZBzz5=ZB311;

    function y = fxyz(x)
        j2r0 = besselj(r,a2*x);
        j2r1 = besselj(r-1,a2*x);
        j2r2 = besselj(r+1,a2*x);
        j2r3 = besselj(r-2,a2*x);
        j2r4 = besselj(r+2,a2*x);

        j0rn0 = besselj(rn,a0*x);
        j0rn1 = besselj(rn-1,a0*x);
        j0rn2 = besselj(rn+1,a0*x);
        j0rn3 = besselj(rn-2,a0*x);
        j0rn4 = besselj(rn+2,a0*x);

        jsn0 = besselj(n,as*x);
        jsn1 = besselj(n-1,as*x);
        jsn2 = besselj(n+1,as*x);

        ej1 = Ei2*j2r1+Ei2c*j2r2;
        ej2 = Ei2*j2r1-Ei2c*j2r2;
        ej3 = e0sc*j0rn1+e0pc*j0rn2;
        ej4 = e0sc*j0rn1-e0pc*j0rn2;
        ej5 = Ei2*(j2r0-j2r3)+Ei2c*(j2r0-j2r4);
        ej6 = Ei2*(j2r0-j2r3)-Ei2c*(j2r0-j2r4);
        ej7 = e0sc*(j0rn0-j0rn3)+e0pc*(j0rn0-j0rn4);
        ej8 = e0sc*(j0rn0-j0rn3)-e0pc*(j0rn0-j0rn4);

        YAmx1 = (1i*(-1/2)).*c2.*ej1.*w0.^(-1).*x.^3.*(2.*e0zc.*j0rn0.*(w0-rn.*wc)+c0.*ej3.*x);
        YAmx2 = (1i*(1/4)).*x.^2.*(a2.*ej3.*ej6.*x+(-2).*ej1.*w0.^(-1).*(ej4.*r.*w0+(-2i).*a0.*e0yc.*j0rn0.*r.*wc.*x+a0.*e0yc.*(a2y.*(-j0rn1+j0rn2)+1i.*a2x.*(j0rn1+j0rn2)).*wc.*x.^2+w0.*(ej3+(2i).*(a2y.*e0xc-a2x.*e0yc).*j0rn0.*x+(-2).*ej3.*x.^2)));
        YAmx3 = (1/4).*w0.^(-1).*x.*((-4i).*e0zc.*ej1.*j0rn0.*rn.*wc+(2i).*(c0.*ej1.*(ej3+ej4.*r)+e0zc.*(a0.*ej1.*(-j0rn1+j0rn2).*r+a2.*ej6.*j0rn0.*rn).*wc).*x+((-1i).*a2.*c0.*ej3.*ej6+4.*ej1.*j0rn0.*(-a2y.*c0.*e0xc+a2x.*c0.*e0yc+(2i).*e0zc.*w0+a0.*a2y.*e0zc.*wc)).*x.^2);
        YBmx1 = (1i*(1/2)).*c0.*c2.*ej1.*ej3.*w2.^(-1).*x.^4;
        YBmx2 = 1i.*c0.*c2.*e0zc.*ej1.*j0rn0.*w2.^(-1).*x.^3;
        YBmx3 = (1/4).*w2.^(-1).*x.^2.*((2i).*ej3.*(ej1+ej2.*rn).*w2+((4i).*c2.*e0zc.*ej1.*j0rn0+(-1i).*a0.*ej1.*ej8.*w2+4.*a2y.*ej3.*j2r0.*rn.*wc).*x+2.*ej1.*ej3.*((-2i).*w2-a0.*a2y.*wc).*x.^2);
        YBmx4 = (1i*(1/4)).*w2.^(-1).*x.*(4.*e0zc.*ej2.*j0rn0.*rn.*w2+(-2).*(c2.*ej3.*(ej1+ej2.*rn)+a0.*e0zc.*ej1.*(-j0rn1+j0rn2).*w2+(4i).*a2y.*e0zc.*j0rn0.*j2r0.*rn.*wc).*x+ej1.*(a0.*c2.*ej8+(-8).*e0zc.*j0rn0.*w2+(4i).*a0.*a2y.*e0zc.*j0rn0.*wc).*x.^2);
        YBmx5 = (1i*(-1/2)).*c2.*e0zc.*w2.^(-1).*x.*(2.*ej2.*j0rn0.*rn+a0.*ej1.*(j0rn1-j0rn2).*x);

        YAmy1 = (-1/2).*c2.*ej2.*w0.^(-1).*x.^3.*(2.*e0zc.*j0rn0.*(w0-rn.*wc)+c0.*ej3.*x);
        YAmy2 = (1/4).*w0.^(-1).*x.^2.*(a2.*ej3.*ej5.*w0.*x+(-2).*ej2.*(ej4.*r.*w0+(-2i).*a0.*e0yc.*j0rn0.*r.*wc.*x+a0.*e0yc.*(a2y.*(-j0rn1+j0rn2)+1i.*a2x.*(j0rn1+j0rn2)).*wc.*x.^2+w0.*(ej3+(2i).*(a2y.*e0xc-a2x.*e0yc).*j0rn0.*x+(-2).*ej3.*x.^2)));
        YAmy3 = (1/4).*w0.^(-1).*x.*((-4).*e0zc.*ej2.*j0rn0.*rn.*wc+2.*(c0.*ej2.*(ej3+ej4.*r)+e0zc.*(a0.*ej2.*(-j0rn1+j0rn2).*r+a2.*ej5.*j0rn0.*rn).*wc).*x-(a2.*c0.*ej3.*ej5+(4i).*ej2.*j0rn0.*(-a2y.*c0.*e0xc+a2x.*c0.*e0yc+(2i).*e0zc.*w0+a0.*a2y.*e0zc.*wc)).*x.^2);
        YBmy1 = (1/2).*c0.*c2.*ej2.*ej3.*w2.^(-1).*x.^4;
        YBmy2 = c0.*c2.*e0zc.*ej2.*j0rn0.*w2.^(-1).*x.^3;
        YBmy3 = (1/4).*w2.^(-1).*x.^2.*(ej2.*(4.*c2.*e0zc.*j0rn0-a0.*ej8.*w2).*x+(-4).*ej3.*j2r0.*(a0.*w2+a2x.*rn.*wc).*x+ej2.*ej3.*w2.*(2+(-4).*x.^2)+2.*ej1.*ej3.*(rn.*w2+a0.*a2x.*wc.*x.^2));
        YBmy4 = (1/4).*w2.^(-1).*x.*(4.*e0zc.*ej1.*j0rn0.*rn.*w2+(-2).*(c2.*ej3.*(ej2+ej1.*rn)+a0.*e0zc.*(ej2.*(-j0rn1+j0rn2)+4.*j0rn0.*j2r0).*w2+4.*a2x.*e0zc.*j0rn0.*j2r0.*rn.*wc).*x+((-8).*e0zc.*ej2.*j0rn0.*w2+a0.*(c2.*ej2.*ej8+4.*c2.*ej3.*j2r0+4.*a2x.*e0zc.*ej1.*j0rn0.*wc)).*x.^2);
        YBmy5 = (1/2).*c2.*e0zc.*w2.^(-1).*x.*((-2).*ej1.*j0rn0.*rn+a0.*(ej2.*(-j0rn1+j0rn2)+4.*j0rn0.*j2r0).*x);

        YAmz1 = (-1i).*c2.*j2r0.*w0.^(-1).*x.^2.*(2.*e0zc.*j0rn0.*(w0-rn.*wc)+c0.*ej3.*x);
        YAmz2 = (-1i).*j2r0.*w0.^(-1).*x.^2.*(2.*e0zc.*j0rn0.*(w0-rn.*wc)+c0.*ej3.*x);
        YAmz3 = (1/2).*w0.^(-1).*x.*((-2i).*ej4.*j2r0.*r.*w0+x.*((-4).*a0.*e0yc.*j0rn0.*j2r0.*r.*wc+2.*a0.*e0yc.*(1i.*a2y.*(j0rn1-j0rn2)+a2x.*(j0rn1+j0rn2)).*j2r0.*wc.*x+w0.*(4.*a2y.*e0xc.*j0rn0.*j2r0+(-4).*a2x.*e0yc.*j0rn0.*j2r0+(-1i).*a2.*ej3.*(j2r1-j2r2)+(4i).*ej3.*j2r0.*x)));
        YAmz4 = (1/2).*w0.^(-1).*x.*((2i).*(c0.*ej4.*j2r0.*r+e0zc.*(a0.*(-j0rn1+j0rn2).*j2r0.*r+a2.*j0rn0.*(-j2r1+j2r2).*rn).*wc)+(4.*a2x.*c0.*e0yc.*j0rn0.*j2r0+1i.*(a2.*c0.*ej3.*(j2r1-j2r2)+8.*e0zc.*j0rn0.*j2r0.*w0)+4.*a2y.*j0rn0.*j2r0.*(-c0.*e0xc+a0.*e0zc.*wc)).*x);
        YBmz1 = 1i.*c0.*ej3.*j2r0.*w2.^(-1).*(w2-r.*wc).*x.^3;
        YBmz2 = (2i).*c0.*e0zc.*j0rn0.*j2r0.*w2.^(-1).*(w2-r.*wc).*x.^2;
        YBmz3 = (2i).*e0zc.*j0rn0.*j2r0.*w2.^(-1).*(w2-r.*wc).*x.^2;
        YBmz4 = (1i*(-1/2)).*w2.^(-1).*x.*(a0.*ej8.*j2r0.*r.*wc.*x+a2.*ej3.*(-j2r1+j2r2).*rn.*wc.*x+(-2).*ej3.*j2r0.*(r.*wc+((-2).*w2+1i.*a0.*a2y.*wc).*x.^2));
        YBmz5 = 1i.*e0zc.*w2.^(-1).*x.*(a2.*j0rn0.*(j2r1-j2r2).*rn.*wc+(-4).*j0rn0.*j2r0.*w2.*x+a0.*j2r0.*wc.*(j0rn1.*r-j0rn2.*r+(2i).*a2y.*j0rn0.*x));

        YAzx1 = (1i*(-1/2)).*c2.*ej1.*w0.^(-1).*x.^2.*(2.*e0zc.*j0rn0.*(w0-rn.*wc)+c0.*ej3.*x);
        YAzx2 = (1i*(1/4)).*x.*(a2.*ej3.*ej6.*x+(-2).*ej1.*w0.^(-1).*(ej4.*r.*w0+(-2i).*a0.*e0yc.*j0rn0.*r.*wc.*x+a0.*e0yc.*(a2y.*(-j0rn1+j0rn2)+1i.*a2x.*(j0rn1+j0rn2)).*wc.*x.^2+w0.*(ej3+(2i).*(a2y.*e0xc-a2x.*e0yc).*j0rn0.*x+(-2).*ej3.*x.^2)));
        YAzx3 = (1/4).*w0.^(-1).*((-4i).*e0zc.*ej1.*j0rn0.*rn.*wc+(2i).*(c0.*ej1.*(ej3+ej4.*r)+e0zc.*(a0.*ej1.*(-j0rn1+j0rn2).*r+a2.*ej6.*j0rn0.*rn).*wc).*x+((-1i).*a2.*c0.*ej3.*ej6+4.*ej1.*j0rn0.*(-a2y.*c0.*e0xc+a2x.*c0.*e0yc+(2i).*e0zc.*w0+a0.*a2y.*e0zc.*wc)).*x.^2);
        YBzx1 = (1i*(1/2)).*c0.*c2.*ej1.*ej3.*w2.^(-1).*x.^3;
        YBzx2 = 1i.*c0.*c2.*e0zc.*ej1.*j0rn0.*w2.^(-1).*x.^2;
        YBzx3 = (1/4).*w2.^(-1).*x.*((2i).*ej3.*(ej1+ej2.*rn).*w2+((4i).*c2.*e0zc.*ej1.*j0rn0+(-1i).*a0.*ej1.*ej8.*w2+4.*a2y.*ej3.*j2r0.*rn.*wc).*x+2.*ej1.*ej3.*((-2i).*w2-a0.*a2y.*wc).*x.^2);
        YBzx4 = (1i*(1/4)).*w2.^(-1).*(4.*e0zc.*ej2.*j0rn0.*rn.*w2+(-2).*(c2.*ej3.*(ej1+ej2.*rn)+a0.*e0zc.*ej1.*(-j0rn1+j0rn2).*w2+(4i).*a2y.*e0zc.*j0rn0.*j2r0.*rn.*wc).*x+ej1.*(a0.*c2.*ej8+(-8).*e0zc.*j0rn0.*w2+(4i).*a0.*a2y.*e0zc.*j0rn0.*wc).*x.^2);
        YBzx5 = (1i*(-1/2)).*c2.*e0zc.*w2.^(-1).*(2.*ej2.*j0rn0.*rn+a0.*ej1.*(j0rn1-j0rn2).*x);

        YAzy1 = (-1/2).*c2.*ej2.*w0.^(-1).*x.^2.*(2.*e0zc.*j0rn0.*(w0-rn.*wc)+c0.*ej3.*x);
        YAzy2 = (1/4).*w0.^(-1).*x.*(a2.*ej3.*ej5.*w0.*x+(-2).*ej2.*(ej4.*r.*w0+(-2i).*a0.*e0yc.*j0rn0.*r.*wc.*x+a0.*e0yc.*(a2y.*(-j0rn1+j0rn2)+1i.*a2x.*(j0rn1+j0rn2)).*wc.*x.^2+w0.*(ej3+(2i).*(a2y.*e0xc-a2x.*e0yc).*j0rn0.*x+(-2).*ej3.*x.^2)));
        YAzy3 = (-1/4).*w0.^(-1).*(4.*e0zc.*ej2.*j0rn0.*rn.*wc+(-2).*(c0.*ej2.*(ej3+ej4.*r)+e0zc.*(a0.*ej2.*(-j0rn1+j0rn2).*r+a2.*ej5.*j0rn0.*rn).*wc).*x+(a2.*c0.*ej3.*ej5+(4i).*ej2.*j0rn0.*(-a2y.*c0.*e0xc+a2x.*c0.*e0yc+(2i).*e0zc.*w0+a0.*a2y.*e0zc.*wc)).*x.^2);
        YBzy1 = (1/2).*c0.*c2.*ej2.*ej3.*w2.^(-1).*x.^3;
        YBzy2 = c0.*c2.*e0zc.*ej2.*j0rn0.*w2.^(-1).*x.^2;
        YBzy3 = (1/4).*w2.^(-1).*x.*(ej2.*(4.*c2.*e0zc.*j0rn0-a0.*ej8.*w2).*x+(-4).*ej3.*j2r0.*(a0.*w2+a2x.*rn.*wc).*x+ej2.*ej3.*w2.*(2+(-4).*x.^2)+2.*ej1.*ej3.*(rn.*w2+a0.*a2x.*wc.*x.^2));
        YBzy4 = (1/4).*w2.^(-1).*(4.*e0zc.*ej1.*j0rn0.*rn.*w2+(-2).*(c2.*ej3.*(ej2+ej1.*rn)+a0.*e0zc.*(ej2.*(-j0rn1+j0rn2)+4.*j0rn0.*j2r0).*w2+4.*a2x.*e0zc.*j0rn0.*j2r0.*rn.*wc).*x+((-8).*e0zc.*ej2.*j0rn0.*w2+a0.*(c2.*ej2.*ej8+4.*c2.*ej3.*j2r0+4.*a2x.*e0zc.*ej1.*j0rn0.*wc)).*x.^2);
        YBzy5 = (1/2).*c2.*e0zc.*w2.^(-1).*((-2).*ej1.*j0rn0.*rn+a0.*(ej2.*(-j0rn1+j0rn2)+4.*j0rn0.*j2r0).*x);

        YAzz1 = (-1i).*c2.*j2r0.*w0.^(-1).*x.*(2.*e0zc.*j0rn0.*(w0-rn.*wc)+c0.*ej3.*x);
        YAzz2 = (-1i).*j2r0.*w0.^(-1).*x.*(2.*e0zc.*j0rn0.*(w0-rn.*wc)+c0.*ej3.*x);
        YAzz3 = (-1i).*ej4.*j2r0.*r+(1/2).*w0.^(-1).*x.*((-4).*a0.*e0yc.*j0rn0.*j2r0.*r.*wc+2.*a0.*e0yc.*(1i.*a2y.*(j0rn1-j0rn2)+a2x.*(j0rn1+j0rn2)).*j2r0.*wc.*x+w0.*(4.*a2y.*e0xc.*j0rn0.*j2r0+(-4).*a2x.*e0yc.*j0rn0.*j2r0+(-1i).*a2.*ej3.*(j2r1-j2r2)+(4i).*ej3.*j2r0.*x));
        YAzz4 = (1/2).*w0.^(-1).*((2i).*(c0.*ej4.*j2r0.*r+e0zc.*(a0.*(-j0rn1+j0rn2).*j2r0.*r+a2.*j0rn0.*(-j2r1+j2r2).*rn).*wc)+(4.*a2x.*c0.*e0yc.*j0rn0.*j2r0+1i.*(a2.*c0.*ej3.*(j2r1-j2r2)+8.*e0zc.*j0rn0.*j2r0.*w0)+4.*a2y.*j0rn0.*j2r0.*(-c0.*e0xc+a0.*e0zc.*wc)).*x);
        YBzz1 = 1i.*c0.*ej3.*j2r0.*w2.^(-1).*(w2-r.*wc).*x.^2;
        YBzz2 = (2i).*c0.*e0zc.*j0rn0.*j2r0.*w2.^(-1).*(w2-r.*wc).*x;
        YBzz3 = (2i).*e0zc.*j0rn0.*j2r0.*w2.^(-1).*(w2-r.*wc).*x;
        YBzz4 = (1i*(1/2)).*w2.^(-1).*(-a0.*ej8.*j2r0.*r.*wc.*x+a2.*ej3.*(j2r1-j2r2).*rn.*wc.*x+2.*ej3.*j2r0.*(r.*wc+((-2).*w2+1i.*a0.*a2y.*wc).*x.^2));
        YBzz5 = 1i.*e0zc.*w2.^(-1).*(a2.*j0rn0.*(j2r1-j2r2).*rn.*wc+(-4).*j0rn0.*j2r0.*w2.*x+a0.*j2r0.*wc.*(j0rn1.*r-j0rn2.*r+(2i).*a2y.*j0rn0.*x));

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

