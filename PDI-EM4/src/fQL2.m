function res = fQL2(ws,n,r)
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
de = n*de2-r*des;
nr = n - r;

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
w2n = w2-n*wc;
wsr = ws-r*wc;
w0nr = w0-nr*wc;
zeta2n = w2n/c2;
zetasr = wsr/cs;
zeta0nr = w0nr/c0;

ZA012=g2(0,zeta2n,1,zetasr,2)/c2^1/cs^2;
ZA011=g2(0,zeta2n,1,zetasr,1)/c2^1/cs^1;
ZA111=g2(1,zeta2n,1,zetasr,1)/c2^1/cs^1;
ZA112=g2(1,zeta2n,1,zetasr,2)/c2^1/cs^2;
ZA212=g2(2,zeta2n,1,zetasr,2)/c2^1/cs^2;
ZA211=g2(2,zeta2n,1,zetasr,1)/c2^1/cs^1;
ZA311=g2(3,zeta2n,1,zetasr,1)/c2^1/cs^1;

ZAmx1=ZA012;ZAmx2=ZA011;ZAmx3=ZA111;
ZAmy1=ZA012;ZAmy2=ZA011;ZAmy3=ZA111;
ZAmz1=ZA112;ZAmz2=ZA011;ZAmz3=ZA111;ZAmz4=ZA211;
ZAzx1=ZA112;ZAzx2=ZA111;ZAzx3=ZA211;
ZAzy1=ZA112;ZAzy2=ZA111;ZAzy3=ZA211;
ZAzz1=ZA212;ZAzz2=ZA111;ZAzz3=ZA211;ZAzz4=ZA311;

ZB012=g2(0,zeta2n,1,zeta0nr,2)/c2^1/c0^2;
ZB011=g2(0,zeta2n,1,zeta0nr,1)/c2^1/c0^1;
ZB112=g2(1,zeta2n,1,zeta0nr,2)/c2^1/c0^2;
ZB212=g2(2,zeta2n,1,zeta0nr,2)/c2^1/c0^2;
ZB111=g2(1,zeta2n,1,zeta0nr,1)/c2^1/c0^1;
ZB211=g2(2,zeta2n,1,zeta0nr,1)/c2^1/c0^1;
ZB311=g2(3,zeta2n,1,zeta0nr,1)/c2^1/c0^1;

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

        j0nr0 = besselj(nr,a0*x);
        j0nr1 = besselj(nr-1,a0*x);
        j0nr2 = besselj(nr+1,a0*x);
        j0nr3 = besselj(nr-2,a0*x);
        j0nr4 = besselj(nr+2,a0*x);

        j2n0 = besselj(n,a2*x);
        j2n1 = besselj(n-1,a2*x);
        j2n2 = besselj(n+1,a2*x);

        ej1 = Eis*jsr1+Eisc*jsr2;
        ej2 = Eis*jsr1-Eisc*jsr2;
        ej3 = e0s*j0nr1+e0p*j0nr2;
        ej4 = e0s*j0nr1-e0p*j0nr2;
        ej5 = Eis*(jsr0-jsr3)+Eisc*(jsr0-jsr4);
        ej6 = Eis*(jsr0-jsr3)-Eisc*(jsr0-jsr4);
        ej7 = e0s*(j0nr0-j0nr3)+e0p*(j0nr0-j0nr4);
        ej8 = e0s*(j0nr0-j0nr3)-e0p*(j0nr0-j0nr4);

        YAmx1 = (1i*(-1/2)).*cs.*ej1.*w0.^(-1).*x.^3.*(2.*e0z.*j0nr0.*(w0-nr.*wc)+c0.*ej3.*x);
        YAmx2 = (1/4).*w0.^(-1).*x.^2.*(1i.*as.*ej3.*ej6.*w0.*x+2.*ej1.*((-1i).*(ej3-ej4.*r).*w0+2.*j0nr0.*(asy.*e0x.*w0-e0y.*(asx.*w0+a0.*r.*wc)).*x+((2i).*ej3.*w0+a0.*e0y.*((-1i).*asy.*(j0nr1-j0nr2)+asx.*(j0nr1+j0nr2)).*wc).*x.^2));
        YAmx3 = (1/4).*w0.^(-1).*x.*((-4i).*e0z.*ej1.*j0nr0.*nr.*wc+(2i).*(c0.*ej1.*(ej3-ej4.*r)+e0z.*(as.*ej6.*j0nr0.*nr+a0.*ej1.*(j0nr1-j0nr2).*r).*wc).*x+((-1i).*as.*c0.*ej3.*ej6+4.*ej1.*j0nr0.*(-asy.*c0.*e0x+asx.*c0.*e0y+(2i).*e0z.*w0+a0.*asy.*e0z.*wc)).*x.^2);
        YBmx1 = (1i*(-1/2)).*c0.*cs.*ej1.*ej3.*ws.^(-1).*x.^4;
        YBmx2 = (-1i).*c0.*cs.*e0z.*ej1.*j0nr0.*ws.^(-1).*x.^3;
        YBmx3 = (1i*(1/4)).*ws.^(-1).*x.^2.*(2.*ej3.*nr.*(ej2.*ws+(-2i).*asy.*jsr0.*wc.*x)+ej1.*((-2).*ej3.*ws+(-4).*cs.*e0z.*j0nr0.*x+a0.*ej8.*ws.*x+2.*ej3.*(1i.*a0.*asy.*wc+2.*ws).*x.^2));
        YBmx4 = (1i*(1/4)).*ws.^(-1).*x.*(-cs.*x.*((-2).*ej1.*ej3+2.*ej2.*ej3.*nr+a0.*ej1.*ej8.*x)+2.*e0z.*(2.*ej2.*j0nr0.*nr.*ws+(2i).*asy.*j0nr0.*wc.*x.*((-2).*jsr0.*nr+a0.*ej1.*x)+ej1.*ws.*x.*(a0.*(-j0nr1+j0nr2)+4.*j0nr0.*x)));
        YBmx5 = (1i*(-1/2)).*cs.*e0z.*ws.^(-1).*x.*(2.*ej2.*j0nr0.*nr+a0.*ej1.*(-j0nr1+j0nr2).*x);

        YAmy1 = (-1/2).*cs.*ej2.*w0.^(-1).*x.^3.*(2.*e0z.*j0nr0.*(w0-nr.*wc)+c0.*ej3.*x);
        YAmy2 = (1/4).*w0.^(-1).*x.^2.*(as.*ej3.*ej5.*w0.*x+2.*ej2.*(-ej3.*w0+ej4.*r.*w0+(2i).*j0nr0.*(-asy.*e0x.*w0+asx.*e0y.*w0+a0.*e0y.*r.*wc).*x+(2.*ej3.*w0+a0.*e0y.*(asy.*(-j0nr1+j0nr2)+(-1i).*asx.*(j0nr1+j0nr2)).*wc).*x.^2));
        YAmy3 = (1/4).*w0.^(-1).*x.*((-4).*e0z.*ej2.*j0nr0.*nr.*wc+2.*(c0.*ej2.*(ej3-ej4.*r)+e0z.*(as.*ej5.*j0nr0.*nr+a0.*ej2.*(j0nr1-j0nr2).*r).*wc).*x-(as.*c0.*ej3.*ej5+(4i).*ej2.*j0nr0.*(-asy.*c0.*e0x+asx.*c0.*e0y+(2i).*e0z.*w0+a0.*asy.*e0z.*wc)).*x.^2);
        YBmy1 = (-1/2).*c0.*cs.*ej2.*ej3.*ws.^(-1).*x.^4;
        YBmy2 = (-1).*c0.*cs.*e0z.*ej2.*j0nr0.*ws.^(-1).*x.^3;
        YBmy3 = (1/4).*ws.^(-1).*x.^2.*(2.*ej3.*(ej1.*nr.*ws+(-2).*jsr0.*(asx.*nr.*wc+a0.*ws).*x+a0.*asx.*ej1.*wc.*x.^2)+ej2.*((-4).*cs.*e0z.*j0nr0.*x+a0.*ej8.*ws.*x+2.*ej3.*ws.*(-1+2.*x.^2)));
        YBmy4 = (1/4).*ws.^(-1).*x.*(cs.*x.*(2.*ej2.*ej3+(-2).*ej1.*ej3.*nr-a0.*ej2.*ej8.*x+4.*a0.*ej3.*jsr0.*x)+2.*e0z.*(a0.*ej2.*(-j0nr1+j0nr2).*ws.*x+(-4).*j0nr0.*x.*(asx.*jsr0.*nr.*wc+a0.*jsr0.*ws-ej2.*ws.*x)+2.*ej1.*j0nr0.*(nr.*ws+a0.*asx.*wc.*x.^2)));
        YBmy5 = (1/2).*cs.*e0z.*ws.^(-1).*x.*((-2).*ej1.*j0nr0.*nr+a0.*(ej2.*(j0nr1-j0nr2)+4.*j0nr0.*jsr0).*x);

        YAmz1 = (-1i).*cs.*jsr0.*w0.^(-1).*x.^2.*(2.*e0z.*j0nr0.*(w0-nr.*wc)+c0.*ej3.*x);
        YAmz2 = (-1i).*jsr0.*w0.^(-1).*x.^2.*(2.*e0z.*j0nr0.*(w0-nr.*wc)+c0.*ej3.*x);
        YAmz3 = (1/2).*w0.^(-1).*x.*((2i).*ej4.*jsr0.*r.*w0+(4.*asy.*e0x.*j0nr0.*jsr0.*w0+(-4).*asx.*e0y.*j0nr0.*jsr0.*w0+(-1i).*as.*ej3.*(jsr1-jsr2).*w0+(-4).*a0.*e0y.*j0nr0.*jsr0.*r.*wc).*x+2.*jsr0.*((2i).*ej3.*w0+a0.*e0y.*((-1i).*asy.*(j0nr1-j0nr2)+asx.*(j0nr1+j0nr2)).*wc).*x.^2);
        YAmz4 = (1/2).*w0.^(-1).*x.*((2i).*e0z.*(as.*j0nr0.*(-jsr1+jsr2).*nr+a0.*(j0nr1-j0nr2).*jsr0.*r).*wc+4.*e0z.*j0nr0.*jsr0.*((2i).*w0+a0.*asy.*wc).*x+c0.*((-2i).*ej4.*jsr0.*r+4.*(-asy.*e0x+asx.*e0y).*j0nr0.*jsr0.*x+1i.*as.*ej3.*(jsr1-jsr2).*x));
        YBmz1 = 1i.*c0.*ej3.*jsr0.*(r.*wc-ws).*ws.^(-1).*x.^3;
        YBmz2 = (2i).*c0.*e0z.*j0nr0.*jsr0.*(r.*wc-ws).*ws.^(-1).*x.^2;
        YBmz3 = (2i).*e0z.*j0nr0.*jsr0.*(r.*wc-ws).*ws.^(-1).*x.^2;
        YBmz4 = (1i*(1/2)).*ws.^(-1).*x.*(a0.*ej8.*jsr0.*r.*wc.*x+ej3.*((-2).*jsr0.*r.*wc+as.*(jsr1-jsr2).*nr.*wc.*x+2.*jsr0.*(1i.*a0.*asy.*wc+2.*ws).*x.^2));
        YBmz5 = 1i.*e0z.*ws.^(-1).*x.*(as.*j0nr0.*(jsr1-jsr2).*nr.*wc+4.*j0nr0.*jsr0.*ws.*x+a0.*jsr0.*wc.*(-j0nr1.*r+j0nr2.*r+(2i).*asy.*j0nr0.*x));

        YAzx1 = (1i*(-1/2)).*cs.*ej1.*w0.^(-1).*x.^2.*(2.*e0z.*j0nr0.*(w0-nr.*wc)+c0.*ej3.*x);
        YAzx2 = (1/4).*w0.^(-1).*x.*(1i.*as.*ej3.*ej6.*w0.*x+2.*ej1.*((-1i).*(ej3-ej4.*r).*w0+2.*j0nr0.*(asy.*e0x.*w0-e0y.*(asx.*w0+a0.*r.*wc)).*x+((2i).*ej3.*w0+a0.*e0y.*((-1i).*asy.*(j0nr1-j0nr2)+asx.*(j0nr1+j0nr2)).*wc).*x.^2));
        YAzx3 = (1/4).*w0.^(-1).*((-4i).*e0z.*ej1.*j0nr0.*nr.*wc+(2i).*(c0.*ej1.*(ej3-ej4.*r)+e0z.*(as.*ej6.*j0nr0.*nr+a0.*ej1.*(j0nr1-j0nr2).*r).*wc).*x+((-1i).*as.*c0.*ej3.*ej6+4.*ej1.*j0nr0.*(-asy.*c0.*e0x+asx.*c0.*e0y+(2i).*e0z.*w0+a0.*asy.*e0z.*wc)).*x.^2);
        YBzx1 = (1i*(-1/2)).*c0.*cs.*ej1.*ej3.*ws.^(-1).*x.^3;
        YBzx2 = (-1i).*c0.*cs.*e0z.*ej1.*j0nr0.*ws.^(-1).*x.^2;
        YBzx3 = (1i*(1/4)).*ws.^(-1).*x.*(2.*ej3.*nr.*(ej2.*ws+(-2i).*asy.*jsr0.*wc.*x)+ej1.*((-2).*ej3.*ws+(-4).*cs.*e0z.*j0nr0.*x+a0.*ej8.*ws.*x+2.*ej3.*(1i.*a0.*asy.*wc+2.*ws).*x.^2));
        YBzx4 = (1i*(1/4)).*ws.^(-1).*(-cs.*x.*((-2).*ej1.*ej3+2.*ej2.*ej3.*nr+a0.*ej1.*ej8.*x)+2.*e0z.*(2.*ej2.*j0nr0.*nr.*ws+(2i).*asy.*j0nr0.*wc.*x.*((-2).*jsr0.*nr+a0.*ej1.*x)+ej1.*ws.*x.*(a0.*(-j0nr1+j0nr2)+4.*j0nr0.*x)));
        YBzx5 = (1i*(-1/2)).*cs.*e0z.*ws.^(-1).*(2.*ej2.*j0nr0.*nr+a0.*ej1.*(-j0nr1+j0nr2).*x);

        YAzy1 = (-1/2).*cs.*ej2.*w0.^(-1).*x.^2.*(2.*e0z.*j0nr0.*(w0-nr.*wc)+c0.*ej3.*x);
        YAzy2 = (1/4).*w0.^(-1).*x.*(as.*ej3.*ej5.*w0.*x+2.*ej2.*(-ej3.*w0+ej4.*r.*w0+(2i).*j0nr0.*(-asy.*e0x.*w0+asx.*e0y.*w0+a0.*e0y.*r.*wc).*x+(2.*ej3.*w0+a0.*e0y.*(asy.*(-j0nr1+j0nr2)+(-1i).*asx.*(j0nr1+j0nr2)).*wc).*x.^2));
        YAzy3 = (-1/4).*w0.^(-1).*(4.*e0z.*ej2.*j0nr0.*nr.*wc+(-2).*(c0.*ej2.*(ej3-ej4.*r)+e0z.*(as.*ej5.*j0nr0.*nr+a0.*ej2.*(j0nr1-j0nr2).*r).*wc).*x+(as.*c0.*ej3.*ej5+(4i).*ej2.*j0nr0.*(-asy.*c0.*e0x+asx.*c0.*e0y+(2i).*e0z.*w0+a0.*asy.*e0z.*wc)).*x.^2);
        YBzy1 = (-1/2).*c0.*cs.*ej2.*ej3.*ws.^(-1).*x.^3;
        YBzy2 = (-1).*c0.*cs.*e0z.*ej2.*j0nr0.*ws.^(-1).*x.^2;
        YBzy3 = (1/4).*ws.^(-1).*x.*(2.*ej3.*(ej1.*nr.*ws+(-2).*jsr0.*(asx.*nr.*wc+a0.*ws).*x+a0.*asx.*ej1.*wc.*x.^2)+ej2.*((-4).*cs.*e0z.*j0nr0.*x+a0.*ej8.*ws.*x+2.*ej3.*ws.*(-1+2.*x.^2)));
        YBzy4 = (1/4).*ws.^(-1).*(cs.*x.*(2.*ej2.*ej3+(-2).*ej1.*ej3.*nr-a0.*ej2.*ej8.*x+4.*a0.*ej3.*jsr0.*x)+2.*e0z.*(a0.*ej2.*(-j0nr1+j0nr2).*ws.*x+(-4).*j0nr0.*x.*(asx.*jsr0.*nr.*wc+a0.*jsr0.*ws-ej2.*ws.*x)+2.*ej1.*j0nr0.*(nr.*ws+a0.*asx.*wc.*x.^2)));
        YBzy5 = (1/2).*cs.*e0z.*ws.^(-1).*((-2).*ej1.*j0nr0.*nr+a0.*(ej2.*(j0nr1-j0nr2)+4.*j0nr0.*jsr0).*x);

        YAzz1 = (-1i).*cs.*jsr0.*w0.^(-1).*x.*(2.*e0z.*j0nr0.*(w0-nr.*wc)+c0.*ej3.*x);
        YAzz2 = (-1i).*jsr0.*w0.^(-1).*x.*(2.*e0z.*j0nr0.*(w0-nr.*wc)+c0.*ej3.*x);
        YAzz3 = (1/2).*w0.^(-1).*((2i).*ej4.*jsr0.*r.*w0+(4.*asy.*e0x.*j0nr0.*jsr0.*w0+(-4).*asx.*e0y.*j0nr0.*jsr0.*w0+(-1i).*as.*ej3.*(jsr1-jsr2).*w0+(-4).*a0.*e0y.*j0nr0.*jsr0.*r.*wc).*x+2.*jsr0.*((2i).*ej3.*w0+a0.*e0y.*((-1i).*asy.*(j0nr1-j0nr2)+asx.*(j0nr1+j0nr2)).*wc).*x.^2);
        YAzz4 = (1/2).*w0.^(-1).*((2i).*e0z.*(as.*j0nr0.*(-jsr1+jsr2).*nr+a0.*(j0nr1-j0nr2).*jsr0.*r).*wc+4.*e0z.*j0nr0.*jsr0.*((2i).*w0+a0.*asy.*wc).*x+c0.*((-2i).*ej4.*jsr0.*r+4.*(-asy.*e0x+asx.*e0y).*j0nr0.*jsr0.*x+1i.*as.*ej3.*(jsr1-jsr2).*x));
        YBzz1 = 1i.*c0.*ej3.*jsr0.*(r.*wc-ws).*ws.^(-1).*x.^2;
        YBzz2 = (2i).*c0.*e0z.*j0nr0.*jsr0.*(r.*wc-ws).*ws.^(-1).*x;
        YBzz3 = (2i).*e0z.*j0nr0.*jsr0.*(r.*wc-ws).*ws.^(-1).*x;
        YBzz4 = (1i*(1/2)).*ws.^(-1).*(a0.*ej8.*jsr0.*r.*wc.*x+ej3.*((-2).*jsr0.*r.*wc+as.*(jsr1-jsr2).*nr.*wc.*x+2.*jsr0.*(1i.*a0.*asy.*wc+2.*ws).*x.^2));
        YBzz5 = 1i.*e0z.*ws.^(-1).*(as.*j0nr0.*(jsr1-jsr2).*nr.*wc+4.*j0nr0.*jsr0.*ws.*x+a0.*jsr0.*wc.*(-j0nr1.*r+j0nr2.*r+(2i).*asy.*j0nr0.*x));

        cmx = exp(-x.^2).*(YAmx1*ZAmx1+YAmx2*ZAmx2+YAmx3*ZAmx3+YBmx1*ZBmx1+YBmx2*ZBmx2+YBmx3*ZBmx3+YBmx4*ZBmx4+YBmx5*ZBmx5);
        cmy = exp(-x.^2).*(YAmy1*ZAmy1+YAmy2*ZAmy2+YAmy3*ZAmy3+YBmy1*ZBmy1+YBmy2*ZBmy2+YBmy3*ZBmy3+YBmy4*ZBmy4+YBmy5*ZBmy5);
        cmz = exp(-x.^2).*(YAmz1*ZAmz1+YAmz2*ZAmz2+YAmz3*ZAmz3+YAmz4*ZAmz4+YBmz1*ZBmz1+YBmz2*ZBmz2+YBmz3*ZBmz3+YBmz4*ZBmz4+YBmz5*ZBmz5);
        czx = exp(-x.^2).*(YAzx1*ZAzx1+YAzx2*ZAzx2+YAzx3*ZAzx3+YBzx1*ZBzx1+YBzx2*ZBzx2+YBzx3*ZBzx3+YBzx4*ZBzx4+YBzx5*ZBzx5);
        czy = exp(-x.^2).*(YAzy1*ZAzy1+YAzy2*ZAzy2+YAzy3*ZAzy3+YBzy1*ZBzy1+YBzy2*ZBzy2+YBzy3*ZBzy3+YBzy4*ZBzy4+YBzy5*ZBzy5);
        czz = exp(-x.^2).*(YAzz1*ZAzz1+YAzz2*ZAzz2+YAzz3*ZAzz3+YAzz4*ZAzz4+YBzz1*ZBzz1+YBzz2*ZBzz2+YBzz3*ZBzz3+YBzz4*ZBzz4+YBzz5*ZBzz5);

        cxx = cmx.*(Ei2c*j2n1+Ei2*j2n2)/2;
        cyx = cmx.*1i.*(Ei2c*j2n1-Ei2*j2n2)/2;
        cxy = cmy.*(Ei2c*j2n1+Ei2*j2n2)/2;
        cyy = cmy.*1i.*(Ei2c*j2n1-Ei2*j2n2)/2;
        cxz = cmz.*(Ei2c*j2n1+Ei2*j2n2)/2;
        cyz = cmz.*1i.*(Ei2c*j2n1-Ei2*j2n2)/2;
        czx = czx*j2n0;
        czy = czy*j2n0;
        czz = czz*j2n0;

        y = [cxx,cxy,cxz;cyx,cyy,cyz;czx,czy,czz];
    end

mat = integral(@(x) fxyz(x),0, inf,'ArrayValued',true);
res = mat*beta0*wp^2*(w0/w2)*exp(1i*de);

end

