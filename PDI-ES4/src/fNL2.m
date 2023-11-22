function res = fNL2(ws,n,p,r)
global wpe wce ve
global nsp n2p n0p
global nsz n2z n0z
global des de2 de0
global yita0 w0

w2 = ws+w0;
wc = wce;
vt = ve;
wp = wpe;

de = (n-r)*de2;
rp = r-p;
np = n-p;

as = nsp*vt/wc;
a2 = n2p*vt/wc;
a0 = n0p*vt/wc;
asx = as*cos(des);
asy = as*sin(des);
a2x = a2*cos(de2);
a2y = a2*sin(de2);

cs = nsz*vt;
c2 = n2z*vt;
c0 = n0z*vt;
zeta2n = (w2-n*wc)/c2;
zetasp = (ws-p*wc)/cs;
zeta2r = (w2-r*wc)/c2;
zeta0rp = (w0-rp*wc)/c0;

ZA1=g3(0,zeta2n,1,zetasp,2,zeta2r,2)/c2^1/cs^2/c2^2;
ZA2=g3(1,zeta2n,1,zetasp,2,zeta2r,2)/c2^1/cs^2/c2^2;
ZA3=g3(0,zeta2n,1,zetasp,2,zeta2r,1)/c2^1/cs^2/c2^1;
ZA4=g3(1,zeta2n,1,zetasp,2,zeta2r,1)/c2^1/cs^2/c2^1;
ZA5=g3(2,zeta2n,1,zetasp,2,zeta2r,1)/c2^1/cs^2/c2^1;
ZA6=g3(0,zeta2n,1,zetasp,1,zeta2r,3)/c2^1/cs^1/c2^3;
ZA7=g3(1,zeta2n,1,zetasp,1,zeta2r,3)/c2^1/cs^1/c2^3;
ZA8=g3(0,zeta2n,1,zetasp,1,zeta2r,2)/c2^1/cs^1/c2^2;
ZA9=g3(1,zeta2n,1,zetasp,1,zeta2r,2)/c2^1/cs^1/c2^2;
ZA10=g3(2,zeta2n,1,zetasp,1,zeta2r,2)/c2^1/cs^1/c2^2;
ZA11=g3(0,zeta2n,1,zetasp,1,zeta2r,1)/c2^1/cs^1/c2^1;
ZA12=g3(1,zeta2n,1,zetasp,1,zeta2r,1)/c2^1/cs^1/c2^1;
ZA13=g3(2,zeta2n,1,zetasp,1,zeta2r,1)/c2^1/cs^1/c2^1;
ZA14=g3(3,zeta2n,1,zetasp,1,zeta2r,1)/c2^1/cs^1/c2^1;

ZB1=g3(0,zeta2n,1,zetasp,2,zeta0rp,2)/c2^1/cs^2/c0^2;
ZB2=g3(1,zeta2n,1,zetasp,2,zeta0rp,2)/c2^1/cs^2/c0^2;
ZB3=g3(0,zeta2n,1,zetasp,2,zeta0rp,1)/c2^1/cs^2/c0^1;
ZB4=g3(1,zeta2n,1,zetasp,2,zeta0rp,1)/c2^1/cs^2/c0^1;
ZB5=g3(2,zeta2n,1,zetasp,2,zeta0rp,1)/c2^1/cs^2/c0^1;
ZB6=g3(0,zeta2n,1,zetasp,1,zeta0rp,3)/c2^1/cs^1/c0^3;
ZB7=g3(1,zeta2n,1,zetasp,1,zeta0rp,3)/c2^1/cs^1/c0^3;
ZB8=g3(0,zeta2n,1,zetasp,1,zeta0rp,2)/c2^1/cs^1/c0^2;
ZB9=g3(1,zeta2n,1,zetasp,1,zeta0rp,2)/c2^1/cs^1/c0^2;
ZB10=g3(2,zeta2n,1,zetasp,1,zeta0rp,2)/c2^1/cs^1/c0^2;
ZB11=g3(0,zeta2n,1,zetasp,1,zeta0rp,1)/c2^1/cs^1/c0^1;
ZB12=g3(1,zeta2n,1,zetasp,1,zeta0rp,1)/c2^1/cs^1/c0^1;
ZB13=g3(2,zeta2n,1,zetasp,1,zeta0rp,1)/c2^1/cs^1/c0^1;
ZB14=g3(3,zeta2n,1,zetasp,1,zeta0rp,1)/c2^1/cs^1/c0^1;


    function y=fx(x)
        j2r0 = besselj(r,a2*x);
        j2r1 = besselj(r-1,a2*x);
        j2r2 = besselj(r+1,a2*x);
        j2r3 = besselj(r-2,a2*x);
        j2r4 = besselj(r+2,a2*x);

        j0rp0 = besselj(rp,a0*x);
        j0rp1 = besselj(rp-1,a0*x);
        j0rp2 = besselj(rp+1,a0*x);
        j0rp3 = besselj(rp-2,a0*x);
        j0rp4 = besselj(rp+2,a0*x);

        j0np0 = besselj(np,a0*x);
        j0np1 = besselj(np-1,a0*x);
        j0np2 = besselj(np+1,a0*x);

        j2n0 = besselj(n,a2*x);

        YA1 = (1/4).*c0.^2.*c2.*cs.*j0np0.*j0rp0.*j2n0.*j2r0.*r.*wc.*x;
        YA2 = (1/4).*c0.^2.*c2.^2.*cs.*j0np0.*j0rp0.*j2n0.*j2r0.*x;
        YA3 = (1/8).*c0.*cs.*j0np0.*j2n0.*(r.*(a0.*(j0rp1-j0rp2).*j2r0.*r+a2.*j0rp0.*(j2r1-j2r2).*rp).*wc.^2+2.*j0rp0.*j2r0.*(c0.*c2+r.*(1i.*a0.*a2y+(-2).*rp).*wc.^2).*x);
        YA4 = (1/8).*c0.*cs.*j0np0.*j2n0.*wc.*(a2.*c2.*j0rp0.*(j2r1-j2r2).*rp+(-4).*j0rp0.*j2r0.*(c0.*r+c2.*rp).*x+a0.*c2.*j2r0.*(j0rp1.*r-j0rp2.*r+2i.*a2y.*j0rp0.*x));
        YA5 = (-1/2).*c0.^2.*c2.*cs.*j0np0.*j0rp0.*j2n0.*j2r0.*x;
        YA6 = (1/2).*c0.^2.*c2.^2.*j0np0.*j0rp0.*j2n0.*j2r0.*r.*wc.*x;
        YA7 = (1/2).*c0.^2.*c2.^3.*j0np0.*j0rp0.*j2n0.*j2r0.*x;
        YA8 = (1/8).*c0.*c2.*j2n0.*(a0.*j2r0.*r.*wc.^2.*((-j0np1+j0np2).*j0rp0.*p+j0np0.*(j0rp1-j0rp2).*(np+r)+2i.*(a2y+asy).*j0np0.*j0rp0.*x)+j0np0.*j0rp0.*(4.*c0.*c2.*j2r0.*x+r.*(np+rp).*wc.^2.*(a2.*(j2r1-j2r2)+(-4).*j2r0.*x)));
        YA9 = (1/8).*c0.*c2.*j2n0.*wc.*(a0.*c2.*j2r0.*((-j0np1+j0np2).*j0rp0.*p+j0np0.*(j0rp1-j0rp2).*(np+r)+2i.*(a2y+asy).*j0np0.*j0rp0.*x)+j0np0.*j0rp0.*(a2.*c2.*(j2r1-j2r2).*(np+rp)+(-4).*j2r0.*(2.*c0.*r+c2.*(np+rp)).*x));
        YA10 = (-1).*c0.^2.*c2.^2.*j0np0.*j0rp0.*j2n0.*j2r0.*x;
        YA11 = (1/16).*j2n0.*wc.*x.^(-2).*(a0.^2.*j2r0.*r.*wc.^2.*x.*(-(j0np1-j0np2).*p.*(j0rp1.*r-j0rp2.*r+2i.*a2y.*j0rp0.*x)+j0np0.*(((-2).*j0rp0+j0rp3+j0rp4).*np.*r+2i.*(j0rp1-j0rp2).*(a2y.*np+asy.*r).*x+(-4).*a2y.*asy.*j0rp0.*x.^2))+j0np0.*j0rp0.*(a2.^2.*((-2).*j2r0+j2r3+j2r4).*np.*r.*rp.*wc.^2.*x+(-8).*j2r0.*(c0.^2.*r+c0.*c2.*(np+rp)+(-2).*np.*r.*rp.*wc.^2).*x.^3+(-2).*a2.*(j2r1-j2r2).*(np.*r.*rp.*wc.^2-(c0.*c2.*(np+rp)+(-4).*np.*r.*rp.*wc.^2).*x.^2))+a0.*(j0np0.*(2.*(-j0rp1+j0rp2).*j2r0.*np.*r.^2.*wc.^2+a2.*(j0rp1-j0rp2).*(j2r1-j2r2).*np.*r.*(r+rp).*wc.^2.*x+2.*(c0.*c2.*(j0rp1-j0rp2).*j2r0.*(np+r)+r.*((-2).*(j0rp1-j0rp2).*j2r0.*np.*(r+rp)+1i.*a2.*j0rp0.*(j2r1-j2r2).*(a2y.*np+asy.*rp)).*wc.^2).*x.^2+4i.*j0rp0.*j2r0.*((a2y+asy).*c0.*c2+(-2).*r.*(a2y.*np+asy.*rp).*wc.^2).*x.^3)+(j0np1-j0np2).*j0rp0.*p.*x.*((-2).*c0.*c2.*j2r0.*x+r.*rp.*wc.^2.*(a2.*(-j2r1+j2r2)+4.*j2r0.*x))));
        YA12 = (1/16).*j2n0.*x.^(-2).*(a0.^2.*c2.*j2r0.*wc.^2.*x.*(-(j0np1-j0np2).*p.*(j0rp1.*r-j0rp2.*r+2i.*a2y.*j0rp0.*x)+j0np0.*(((-2).*j0rp0+j0rp3+j0rp4).*np.*r+2i.*(j0rp1-j0rp2).*(a2y.*np+asy.*r).*x+(-4).*a2y.*asy.*j0rp0.*x.^2))+j0np0.*j0rp0.*(a2.^2.*c2.*((-2).*j2r0+j2r3+j2r4).*np.*rp.*wc.^2.*x+8.*j2r0.*((-3).*c0.^2.*c2+2.*(c2.*np.*rp+c0.*r.*(np+rp)).*wc.^2).*x.^3+(-2).*a2.*(j2r1-j2r2).*wc.^2.*(2.*c0.*r.*(np+rp).*x.^2+c2.*np.*rp.*(1+4.*x.^2)))-a0.*wc.^2.*(4.*c0.*j2r0.*r.*x.^2.*((-j0np1+j0np2).*j0rp0.*p+j0np0.*(j0rp1-j0rp2).*(np+r)+2i.*(a2y+asy).*j0np0.*j0rp0.*x)+c2.*((j0np1-j0np2).*j0rp0.*p.*rp.*x.*(a2.*(j2r1-j2r2)+(-4).*j2r0.*x)+j0np0.*(2.*(j0rp1-j0rp2).*j2r0.*np.*r-a2.*(j0rp1-j0rp2).*(j2r1-j2r2).*np.*(r+rp).*x+2.*(2.*(j0rp1-j0rp2).*j2r0.*np.*(r+rp)+(-1i).*a2.*j0rp0.*(j2r1-j2r2).*(a2y.*np+asy.*rp)).*x.^2+(1i*8).*j0rp0.*j2r0.*(a2y.*np+asy.*rp).*x.^3))));
        YA13 = (1/4).*c0.*j2n0.*wc.*(a0.*c2.*j2r0.*((j0np1-j0np2).*j0rp0.*p-j0np0.*(j0rp1-j0rp2).*(np+r)+(-2i).*(a2y+asy).*j0np0.*j0rp0.*x)+j0np0.*j0rp0.*(-a2.*c2.*(j2r1-j2r2).*(np+rp)+4.*j2r0.*(c0.*r+c2.*(np+rp)).*x));
        YA14 = c0.^2.*c2.*j0np0.*j0rp0.*j2n0.*j2r0.*x;
        YB1 = (-1/4).*c0.^2.*c2.*cs.*j0np0.*j0rp0.*j2n0.*j2r0.*rp.*wc.*x;
        YB2 = (-1/4).*c0.^3.*c2.*cs.*j0np0.*j0rp0.*j2n0.*j2r0.*x;
        YB3 = (1/8).*c0.*cs.*j0np0.*j2n0.*(rp.*(a0.*(-j0rp1+j0rp2).*j2r0.*r+a2.*j0rp0.*(-j2r1+j2r2).*rp).*wc.^2+(-2).*j0rp0.*j2r0.*(c0.*c2+(1i.*a0.*a2y+(-2).*r).*rp.*wc.^2).*x);
        YB4 = (1/8).*c0.*cs.*j0np0.*j2n0.*wc.*(a2.*c0.*j0rp0.*(-j2r1+j2r2).*rp+4.*j0rp0.*j2r0.*(c0.*r+c2.*rp).*x+a0.*c0.*j2r0.*(-j0rp1.*r+j0rp2.*r+(-2i).*a2y.*j0rp0.*x));
        YB5 = (1/2).*c0.^2.*c2.*cs.*j0np0.*j0rp0.*j2n0.*j2r0.*x;
        YB6 = (-1/2).*c0.^3.*c2.*j0np0.*j0rp0.*j2n0.*j2r0.*rp.*wc.*x;
        YB7 = (-1/2).*c0.^4.*c2.*j0np0.*j0rp0.*j2n0.*j2r0.*x;
        YB8 = (1/8).*c0.*j2n0.*(rp.*(a0.*j2r0.*(c2.*j0np0.*(-j0rp1+j0rp2).*np+c2.*(j0np1-j0np2).*j0rp0.*p+c0.*j0np0.*(-j0rp1+j0rp2).*r)-a2.*j0np0.*j0rp0.*(j2r1-j2r2).*(c2.*np+c0.*rp)).*wc.^2+(-2).*j0np0.*j0rp0.*j2r0.*(2.*c0.^2.*c2+1i.*(a0.*(a2y.*c0+asy.*c2)+2i.*(c2.*np+c0.*r)).*rp.*wc.^2).*x);
        YB9 = (-1/8).*c0.^2.*j2n0.*wc.*(j0np0.*j0rp0.*(a2.*(j2r1-j2r2).*(c2.*np+c0.*rp)+(-4).*j2r0.*(c2.*np+c0.*r+2.*c2.*rp).*x)+a0.*j2r0.*(c2.*(-j0np1+j0np2).*j0rp0.*p+c0.*j0np0.*(j0rp1.*r-j0rp2.*r+2i.*a2y.*j0rp0.*x)+c2.*j0np0.*(j0rp1.*np-j0rp2.*np+2i.*asy.*j0rp0.*x)));
        YB10 = c0.^3.*c2.*j0np0.*j0rp0.*j2n0.*j2r0.*x;
        YB11 = (1/16).*j2n0.*wc.*x.^(-2).*(a0.^2.*j2r0.*rp.*wc.^2.*x.*((j0np1-j0np2).*p.*(j0rp1.*r-j0rp2.*r+2i.*a2y.*j0rp0.*x)+j0np0.*(-((-2).*j0rp0+j0rp3+j0rp4).*np.*r+(-2i).*(j0rp1-j0rp2).*(a2y.*np+asy.*r).*x+4.*a2y.*asy.*j0rp0.*x.^2))+j0np0.*j0rp0.*(a2.^2.*(2.*j2r0-j2r3-j2r4).*np.*rp.^2.*wc.^2.*x+8.*j2r0.*(c0.^2.*r+c0.*c2.*(np+rp)+(-2).*np.*r.*rp.*wc.^2).*x.^3+2.*a2.*(j2r1-j2r2).*(np.*rp.^2.*wc.^2-(c0.*c2.*np+c0.^2.*rp+(-2).*np.*rp.*(r+rp).*wc.^2).*x.^2))+a0.*(j0np0.*(2.*(j0rp1-j0rp2).*j2r0.*np.*r.*rp.*wc.^2-a2.*(j0rp1-j0rp2).*(j2r1-j2r2).*np.*rp.*(r+rp).*wc.^2.*x+2.*(-c0.*(j0rp1-j0rp2).*j2r0.*(c2.*np+c0.*r)+(-1i).*rp.*(4i.*(j0rp1-j0rp2).*j2r0.*np.*r+a2.*j0rp0.*(j2r1-j2r2).*(a2y.*np+asy.*rp)).*wc.^2).*x.^2+(-4i).*j0rp0.*j2r0.*(c0.*(a2y.*c0+asy.*c2)+(-2).*(a2y.*np+asy.*r).*rp.*wc.^2).*x.^3)+(j0np1-j0np2).*j0rp0.*p.*x.*(2.*c0.*c2.*j2r0.*x+rp.*wc.^2.*(a2.*(j2r1-j2r2).*rp+(-4).*j2r0.*r.*x))));
        YB12 = (1/16).*j2n0.*x.^(-2).*(a0.*wc.^2.*(2.*c0.*j0np0.*(j0rp1-j0rp2).*j2r0.*np.*r-a2.*c0.*(j2r1-j2r2).*((-j0np1+j0np2).*j0rp0.*p.*rp+j0np0.*(j0rp1-j0rp2).*np.*(r+rp)).*x+((-2i).*a2.*c0.*j0np0.*j0rp0.*(j2r1-j2r2).*(a2y.*np+asy.*rp)+4.*j2r0.*(c2.*(j0np0.*(j0rp1-j0rp2).*np+(-j0np1+j0np2).*j0rp0.*p).*rp+c0.*r.*((-j0np1+j0np2).*j0rp0.*p+j0np0.*(j0rp1-j0rp2).*(2.*np+rp)))).*x.^2+(1i*8).*j0np0.*j0rp0.*j2r0.*(a2y.*c0.*(np+rp)+asy.*(c0.*r+c2.*rp)).*x.^3)+a0.^2.*c0.*j2r0.*wc.^2.*x.*((j0np1-j0np2).*p.*(j0rp1.*r-j0rp2.*r+2i.*a2y.*j0rp0.*x)+j0np0.*(-((-2).*j0rp0+j0rp3+j0rp4).*np.*r+(-2i).*(j0rp1-j0rp2).*(a2y.*np+asy.*r).*x+4.*a2y.*asy.*j0rp0.*x.^2))+j0np0.*j0rp0.*(a2.^2.*c0.*(2.*j2r0-j2r3-j2r4).*np.*rp.*wc.^2.*x+8.*j2r0.*(3.*c0.^2.*c2+(-2).*(c2.*np.*rp+c0.*r.*(np+rp)).*wc.^2).*x.^3+2.*a2.*(j2r1-j2r2).*wc.^2.*(c0.*np.*rp+2.*(c0.*np.*r+(c0+c2).*np.*rp+c0.*rp.^2).*x.^2)));
        YB13 = (1/4).*c0.*j2n0.*wc.*(j0np0.*j0rp0.*(a2.*(j2r1-j2r2).*(c2.*np+c0.*rp)+(-4).*j2r0.*(c0.*r+c2.*(np+rp)).*x)+a0.*j2r0.*(c2.*(-j0np1+j0np2).*j0rp0.*p+c0.*j0np0.*(j0rp1.*r-j0rp2.*r+2i.*a2y.*j0rp0.*x)+c2.*j0np0.*(j0rp1.*np-j0rp2.*np+2i.*asy.*j0rp0.*x)));
        YB14 = (-1).*c0.^2.*c2.*j0np0.*j0rp0.*j2n0.*j2r0.*x;

        y = exp(-x.^2).*(YA1*ZA1+YA2*ZA2+YA3*ZA3+YA4*ZA4+YA5*ZA5+YA6*ZA6+YA7*ZA7+YA8*ZA8+YA9*ZA9+YA10*ZA10+YA11*ZA11+YA12*ZA12+YA13*ZA13+YA14*ZA14+YB1*ZB1+YB2*ZB2+YB3*ZB3+YB4*ZB4+YB5*ZB5+YB6*ZB6+YB7*ZB7+YB8*ZB8+YB9*ZB9+YB10*ZB10+YB11*ZB11+YB12*ZB12+YB13*ZB13+YB14*ZB14);
    end

in = integral(@(x)fx(x),0,inf);
n2 = sqrt(n2p^2+n2z^2);
res = in*yita0^2*(wp/n2/vt)^2*exp(1i*de);
end

