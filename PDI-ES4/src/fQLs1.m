function res = fQLs1(ws,n,r)
global wpe wce ve
global nsp n1p n0p
global nsz n1z n0z
global des de1 de0
global yita0 w0

w1 = ws-w0;
wc = wce;
vt = ve;
wp = wpe;

de = n*des-r*de1;
nr = n-r;

as = nsp*vt/wc;
a1 = n1p*vt/wc;
a0 = n0p*vt/wc;
a1x = a1*cos(de1);
a1y = a1*sin(de1);

cs = nsz*vt;
c1 = n1z*vt;
c0 = n0z*vt;
zetasn = (ws-n*wc)/cs;
zeta1r = (w1-r*wc)/c1;
zeta0nr = (w0-nr*wc)/c0;

ZA1=g2(0,zetasn,1,zeta1r,2)/cs^1/c1^2;
ZA2=g2(1,zetasn,1,zeta1r,2)/cs^1/c1^2;
ZA3=g2(0,zetasn,1,zeta1r,1)/cs^1/c1^1;
ZA4=g2(1,zetasn,1,zeta1r,1)/cs^1/c1^1;
ZA5=g2(2,zetasn,1,zeta1r,1)/cs^1/c1^1;
ZB1=g2(0,zetasn,1,zeta0nr,2)/cs^1/c0^2;
ZB2=g2(1,zetasn,1,zeta0nr,2)/cs^1/c0^2;
ZB3=g2(0,zetasn,1,zeta0nr,1)/cs^1/c0^1;
ZB4=g2(1,zetasn,1,zeta0nr,1)/cs^1/c0^1;
ZB5=g2(2,zetasn,1,zeta0nr,1)/cs^1/c0^1;

    function y=fx(x)
        j1r0 = besselj(r,a1*x);
        j1r1 = besselj(r-1,a1*x);
        j1r2 = besselj(r+1,a1*x);
        
        j0nr0 = besselj(nr,a0*x);
        j0nr1 = besselj(nr-1,a0*x);
        j0nr2 = besselj(nr+1,a0*x);
        
        jsn0 = besselj(n,as*x);
        
        YA1 = (-1).*c0.*c1.*j0nr0.*j1r0.*jsn0.*r.*wc.*x;
        YA2 = (-1).*c0.*c1.^2.*j0nr0.*j1r0.*jsn0.*x;
        YA3 = (1/2).*jsn0.*r.*(a1.*j0nr0.*(-j1r1+j1r2).*nr+a0.*(j0nr1-j0nr2).*j1r0.*r).*wc.^2-j0nr0.*j1r0.*jsn0.*(c0.*c1+(1i.*a0.*a1y+(-2).*nr).*r.*wc.^2).*x;
        YA4 = (1/2).*jsn0.*wc.*(a1.*c1.*j0nr0.*(-j1r1+j1r2).*nr+4.*j0nr0.*j1r0.*(c1.*nr+c0.*r).*x+a0.*c1.*j1r0.*(j0nr1.*r-j0nr2.*r+(-2i).*a1y.*j0nr0.*x));
        YA5 = 2.*c0.*c1.*j0nr0.*j1r0.*jsn0.*x;
        YB1 = (-1).*c0.*c1.*j0nr0.*j1r0.*jsn0.*nr.*wc.*x;
        YB2 = (-1).*c0.^2.*c1.*j0nr0.*j1r0.*jsn0.*x;
        YB3 = (1/2).*jsn0.*nr.*(a1.*j0nr0.*(j1r1-j1r2).*nr+a0.*(-j0nr1+j0nr2).*j1r0.*r).*wc.^2+j0nr0.*j1r0.*jsn0.*(-c0.*c1+nr.*(1i.*a0.*a1y+2.*r).*wc.^2).*x;
        YB4 = (1/2).*jsn0.*wc.*(a1.*c0.*j0nr0.*(j1r1-j1r2).*nr+4.*j0nr0.*j1r0.*(c1.*nr+c0.*r).*x+a0.*c0.*j1r0.*(-j0nr1.*r+j0nr2.*r+2i.*a1y.*j0nr0.*x));
        YB5 = 2.*c0.*c1.*j0nr0.*j1r0.*jsn0.*x;
        
        y = exp(-x.^2).*(YA1*ZA1+YA2*ZA2+YA3*ZA3+YA4*ZA4+YA5*ZA5+YB1*ZB1+YB2*ZB2+YB3*ZB3+YB4*ZB4+YB5*ZB5);
    end

in = integral(@(x)fx(x),0,inf);
ns = sqrt(nsp^2+nsz^2);
res = in*yita0*(wp/ns/vt)^2*exp(1i*de);
end