function res = fQL2(ws,n,r)
global wpe wce ve
global nsp n2p n0p
global nsz n2z n0z
global des de2 de0
global yita0 w0

w2 = ws+w0;
wc = wce;
vt = ve;
wp = wpe;

de = n*de2-r*des;
nr = n - r;

as = nsp*vt/wc;
a2 = n2p*vt/wc;
a0 = n0p*vt/wc;
asx = as*cos(des);
asy = as*sin(des);

cs = nsz*vt;
c2 = n2z*vt;
c0 = n0z*vt;
zeta2n = (w2-n*wc)/c2;
zetasr = (ws-r*wc)/cs;
zeta0nr = (w0-nr*wc)/c0;

ZA1=g2(0,zeta2n,1,zetasr,2)/c2^1/cs^2;
ZA2=g2(1,zeta2n,1,zetasr,2)/c2^1/cs^2;
ZA3=g2(0,zeta2n,1,zetasr,1)/c2^1/cs^1;
ZA4=g2(1,zeta2n,1,zetasr,1)/c2^1/cs^1;
ZA5=g2(2,zeta2n,1,zetasr,1)/c2^1/cs^1;
ZB1=g2(0,zeta2n,1,zeta0nr,2)/c2^1/c0^2;
ZB2=g2(1,zeta2n,1,zeta0nr,2)/c2^1/c0^2;
ZB3=g2(0,zeta2n,1,zeta0nr,1)/c2^1/c0^1;
ZB4=g2(1,zeta2n,1,zeta0nr,1)/c2^1/c0^1;
ZB5=g2(2,zeta2n,1,zeta0nr,1)/c2^1/c0^1;

    function y=fx(x)
        jsr0 = besselj(r,as*x);
        jsr1 = besselj(r-1,as*x);
        jsr2 = besselj(r+1,as*x);

        j0nr0 = besselj(nr,a0*x);
        j0nr1 = besselj(nr-1,a0*x);
        j0nr2 = besselj(nr+1,a0*x);

        j2n0 = besselj(n,a2*x);

        YA1 = (-1).*c0.*cs.*j0nr0.*j2n0.*jsr0.*r.*wc.*x;
        YA2 = (-1).*c0.*cs.^2.*j0nr0.*j2n0.*jsr0.*x;
        YA3 = (-1/2).*j2n0.*(as.*j0nr0.*(jsr1-jsr2).*nr.*r.*wc.^2+2.*j0nr0.*jsr0.*(c0.*cs+(-2).*nr.*r.*wc.^2).*x+a0.*jsr0.*r.*wc.^2.*(-j0nr1.*r+j0nr2.*r+2i.*asy.*j0nr0.*x));
        YA4 = (-1/2).*j2n0.*wc.*(as.*cs.*j0nr0.*(jsr1-jsr2).*nr+(-4).*j0nr0.*jsr0.*(cs.*nr+c0.*r).*x+a0.*cs.*jsr0.*(-j0nr1.*r+j0nr2.*r+2i.*asy.*j0nr0.*x));
        YA5 = 2.*c0.*cs.*j0nr0.*j2n0.*jsr0.*x;
        YB1 = (-1).*c0.*cs.*j0nr0.*j2n0.*jsr0.*nr.*wc.*x;
        YB2 = (-1).*c0.^2.*cs.*j0nr0.*j2n0.*jsr0.*x;
        YB3 = (1/2).*j2n0.*(as.*j0nr0.*(jsr1-jsr2).*nr.^2.*wc.^2+(-2).*j0nr0.*jsr0.*(c0.*cs+(-2).*nr.*r.*wc.^2).*x+a0.*jsr0.*nr.*wc.^2.*(-j0nr1.*r+j0nr2.*r+2i.*asy.*j0nr0.*x));
        YB4 = (1/2).*j2n0.*wc.*(as.*c0.*j0nr0.*(jsr1-jsr2).*nr+4.*j0nr0.*jsr0.*(cs.*nr+c0.*r).*x+a0.*c0.*jsr0.*(-j0nr1.*r+j0nr2.*r+2i.*asy.*j0nr0.*x));
        YB5 = 2.*c0.*cs.*j0nr0.*j2n0.*jsr0.*x;

        y = exp(-x.^2).*(YA1*ZA1+YA2*ZA2+YA3*ZA3+YA4*ZA4+YA5*ZA5+YB1*ZB1+YB2*ZB2+YB3*ZB3+YB4*ZB4+YB5*ZB5);
    end

in = integral(@(x)fx(x),0,inf);
n2 = sqrt(n2p^2+n2z^2);
res = in*yita0*(wp/n2/vt)^2*exp(1i*de);
end

