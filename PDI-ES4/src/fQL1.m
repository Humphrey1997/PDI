function res = fQL1(ws,n,r)
global wpe wce ve
global nsp n1p n0p
global nsz n1z n0z
global des de1 de0
global yita0 w0

w1 = ws-w0;
wc = wce;
vt = ve;
wp = wpe;

de = n*de1-r*des;
rn = r-n;

as = nsp*vt/wc;
a1 = n1p*vt/wc;
a0 = n0p*vt/wc;
asx = as*cos(des);
asy = as*sin(des);

cs = nsz*vt;
c1 = n1z*vt;
c0 = n0z*vt;
zeta1n = (w1-n*wc)/c1;
zetasr = (ws-r*wc)/cs;
zeta0rn = (w0-rn*wc)/c0;

ZA1=g2(0,zeta1n,1,zetasr,2)/c1^1/cs^2;
ZA2=g2(1,zeta1n,1,zetasr,2)/c1^1/cs^2;
ZA3=g2(0,zeta1n,1,zetasr,1)/c1^1/cs^1;
ZA4=g2(1,zeta1n,1,zetasr,1)/c1^1/cs^1;
ZA5=g2(2,zeta1n,1,zetasr,1)/c1^1/cs^1;
ZB1=g2(0,zeta1n,1,zeta0rn,2)/c1^1/c0^2;
ZB2=g2(1,zeta1n,1,zeta0rn,2)/c1^1/c0^2;
ZB3=g2(0,zeta1n,1,zeta0rn,1)/c1^1/c0^1;
ZB4=g2(1,zeta1n,1,zeta0rn,1)/c1^1/c0^1;
ZB5=g2(2,zeta1n,1,zeta0rn,1)/c1^1/c0^1;

    function y=fx(x)
        jsr0 = besselj(r,as*x);
        jsr1 = besselj(r-1,as*x);
        jsr2 = besselj(r+1,as*x);
        
        j0rn0 = besselj(rn,a0*x);
        j0rn1 = besselj(rn-1,a0*x);
        j0rn2 = besselj(rn+1,a0*x);
        
        j1n0 = besselj(n,a1*x);
        
        YA1 = c0.*cs.*j0rn0.*j1n0.*jsr0.*r.*wc.*x;
        YA2 = c0.*cs.^2.*j0rn0.*j1n0.*jsr0.*x;
        YA3 = (1/2).*j1n0.*r.*(a0.*(j0rn1-j0rn2).*jsr0.*r+as.*j0rn0.*(jsr1-jsr2).*rn).*wc.^2+j0rn0.*j1n0.*jsr0.*(c0.*cs+r.*(1i.*a0.*asy+(-2).*rn).*wc.^2).*x;
        YA4 = (1/2).*j1n0.*wc.*(as.*cs.*j0rn0.*(jsr1-jsr2).*rn+(-4).*j0rn0.*jsr0.*(c0.*r+cs.*rn).*x+a0.*cs.*jsr0.*(j0rn1.*r-j0rn2.*r+2i.*asy.*j0rn0.*x));
        YA5 = (-2).*c0.*cs.*j0rn0.*j1n0.*jsr0.*x;
        YB1 = (-1).*c0.*cs.*j0rn0.*j1n0.*jsr0.*rn.*wc.*x;
        YB2 = (-1).*c0.^2.*cs.*j0rn0.*j1n0.*jsr0.*x;
        YB3 = (1/2).*j1n0.*rn.*(a0.*(-j0rn1+j0rn2).*jsr0.*r+as.*j0rn0.*(-jsr1+jsr2).*rn).*wc.^2-j0rn0.*j1n0.*jsr0.*(c0.*cs+(1i.*a0.*asy+(-2).*r).*rn.*wc.^2).*x;
        YB4 = (1/2).*j1n0.*wc.*(as.*c0.*j0rn0.*(-jsr1+jsr2).*rn+4.*j0rn0.*jsr0.*(c0.*r+cs.*rn).*x+a0.*c0.*jsr0.*(-j0rn1.*r+j0rn2.*r+(-2i).*asy.*j0rn0.*x));
        YB5 = 2.*c0.*cs.*j0rn0.*j1n0.*jsr0.*x;
        
        y = exp(-x.^2).*(YA1*ZA1+YA2*ZA2+YA3*ZA3+YA4*ZA4+YA5*ZA5+YB1*ZB1+YB2*ZB2+YB3*ZB3+YB4*ZB4+YB5*ZB5);
    end

in = integral(@(x)fx(x),0,inf);
n1 = sqrt(n1p^2+n1z^2);
res = in*yita0*(wp/n1/vt)^2*exp(1i*de);
end

