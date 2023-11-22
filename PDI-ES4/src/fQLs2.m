function res = fQLs2(ws,n,r)
global wpe wce ve
global nsp n2p n0p
global nsz n2z n0z
global des de2 de0
global yita0 w0

w2 = ws+w0;
wc = wce;
vt = ve;
wp = wpe;

de = n*des-r*de2;
rn = r-n;

as = nsp*vt/wc;
a2 = n2p*vt/wc;
a0 = n0p*vt/wc;
a2x = a2*cos(de2);
a2y = a2*sin(de2);

cs = nsz*vt;
c2 = n2z*vt;
c0 = n0z*vt;
zetasn = (ws-n*wc)/cs;
zeta2r = (w2-r*wc)/c2;
zeta0rn = (w0-rn*wc)/c0;

ZA1=g2(0,zetasn,1,zeta2r,2)/cs^1/c2^2;
ZA2=g2(1,zetasn,1,zeta2r,2)/cs^1/c2^2;
ZA3=g2(0,zetasn,1,zeta2r,1)/cs^1/c2^1;
ZA4=g2(1,zetasn,1,zeta2r,1)/cs^1/c2^1;
ZA5=g2(2,zetasn,1,zeta2r,1)/cs^1/c2^1;
ZB1=g2(0,zetasn,1,zeta0rn,2)/cs^1/c0^2;
ZB2=g2(1,zetasn,1,zeta0rn,2)/cs^1/c0^2;
ZB3=g2(0,zetasn,1,zeta0rn,1)/cs^1/c0^1;
ZB4=g2(1,zetasn,1,zeta0rn,1)/cs^1/c0^1;
ZB5=g2(2,zetasn,1,zeta0rn,1)/cs^1/c0^1;

    function y=fx(x)
        j2r0 = besselj(r,a2*x);
        j2r1 = besselj(r-1,a2*x);
        j2r2 = besselj(r+1,a2*x);

        j0rn0 = besselj(rn,a0*x);
        j0rn1 = besselj(rn-1,a0*x);
        j0rn2 = besselj(rn+1,a0*x);

        jsn0 = besselj(n,as*x);

        YA1 = c0.*c2.*j0rn0.*j2r0.*jsn0.*r.*wc.*x;
        YA2 = c0.*c2.^2.*j0rn0.*j2r0.*jsn0.*x;
        YA3 = (1/2).*jsn0.*(a2.*j0rn0.*(j2r1-j2r2).*r.*rn.*wc.^2+2.*j0rn0.*j2r0.*(c0.*c2+(-2).*r.*rn.*wc.^2).*x+a0.*j2r0.*r.*wc.^2.*(j0rn1.*r-j0rn2.*r+2i.*a2y.*j0rn0.*x));
        YA4 = (1/2).*jsn0.*wc.*(a2.*c2.*j0rn0.*(j2r1-j2r2).*rn+(-4).*j0rn0.*j2r0.*(c0.*r+c2.*rn).*x+a0.*c2.*j2r0.*(j0rn1.*r-j0rn2.*r+2i.*a2y.*j0rn0.*x));
        YA5 = (-2).*c0.*c2.*j0rn0.*j2r0.*jsn0.*x;
        YB1 = (-1).*c0.*c2.*j0rn0.*j2r0.*jsn0.*rn.*wc.*x;
        YB2 = (-1).*c0.^2.*c2.*j0rn0.*j2r0.*jsn0.*x;
        YB3 = (-1/2).*jsn0.*(a2.*j0rn0.*(j2r1-j2r2).*rn.^2.*wc.^2+2.*j0rn0.*j2r0.*(c0.*c2+(-2).*r.*rn.*wc.^2).*x+a0.*j2r0.*rn.*wc.^2.*(j0rn1.*r-j0rn2.*r+2i.*a2y.*j0rn0.*x));
        YB4 = (-1/2).*jsn0.*wc.*(a2.*c0.*j0rn0.*(j2r1-j2r2).*rn+(-4).*j0rn0.*j2r0.*(c0.*r+c2.*rn).*x+a0.*c0.*j2r0.*(j0rn1.*r-j0rn2.*r+2i.*a2y.*j0rn0.*x));
        YB5 = 2.*c0.*c2.*j0rn0.*j2r0.*jsn0.*x;

        y = exp(-x.^2).*(YA1*ZA1+YA2*ZA2+YA3*ZA3+YA4*ZA4+YA5*ZA5+YB1*ZB1+YB2*ZB2+YB3*ZB3+YB4*ZB4+YB5*ZB5);
    end

in = integral(@(x)fx(x),0,inf);
ns = sqrt(nsp^2+nsz^2);
res = in*yita0*(wp/ns/vt)^2*exp(1i*de);
end