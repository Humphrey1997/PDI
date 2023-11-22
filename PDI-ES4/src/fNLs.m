function res = fNLs(ws,n,p,r)
global wpe wce ve
global nsp n1p n0p
global nsz n1z n0z
global des de1 de0
global yita0 w0

w1 = ws-w0;
wc = wce;
vt = ve;
wp = wpe;

de = (n-r)*des;
rp = r-p;
np = n-p;
as = nsp*vt/wc;
a1 = n1p*vt/wc;
a0 = n0p*vt/wc;
asx = as*cos(des);
asy = as*sin(des);
a1x = a1*cos(de1);
a1y = a1*sin(de1);

cs = nsz*vt;
c1 = n1z*vt;
c0 = n0z*vt;
zetasn = (ws-n*wc)/cs;
zeta1p = (w1-p*wc)/c1;
zetasr = (ws-r*wc)/cs;
zeta0rp = (w0-rp*wc)/c0;

ZA1=g3(0,zetasn,1,zeta1p,2,zetasr,2)/cs^1/c1^2/cs^2;
ZA2=g3(1,zetasn,1,zeta1p,2,zetasr,2)/cs^1/c1^2/cs^2;
ZA3=g3(0,zetasn,1,zeta1p,2,zetasr,1)/cs^1/c1^2/cs^1;
ZA4=g3(1,zetasn,1,zeta1p,2,zetasr,1)/cs^1/c1^2/cs^1;
ZA5=g3(2,zetasn,1,zeta1p,2,zetasr,1)/cs^1/c1^2/cs^1;
ZA6=g3(0,zetasn,1,zeta1p,1,zetasr,3)/cs^1/c1^1/cs^3;
ZA7=g3(1,zetasn,1,zeta1p,1,zetasr,3)/cs^1/c1^1/cs^3;
ZA8=g3(0,zetasn,1,zeta1p,1,zetasr,2)/cs^1/c1^1/cs^2;
ZA9=g3(1,zetasn,1,zeta1p,1,zetasr,2)/cs^1/c1^1/cs^2;
ZA10=g3(2,zetasn,1,zeta1p,1,zetasr,2)/cs^1/c1^1/cs^2;
ZA11=g3(0,zetasn,1,zeta1p,1,zetasr,1)/cs^1/c1^1/cs^1;
ZA12=g3(1,zetasn,1,zeta1p,1,zetasr,1)/cs^1/c1^1/cs^1;
ZA13=g3(2,zetasn,1,zeta1p,1,zetasr,1)/cs^1/c1^1/cs^1;
ZA14=g3(3,zetasn,1,zeta1p,1,zetasr,1)/cs^1/c1^1/cs^1;

ZB1=g3(0,zetasn,1,zeta1p,2,zeta0rp,2)/cs^1/c1^2/c0^2;
ZB2=g3(1,zetasn,1,zeta1p,2,zeta0rp,2)/cs^1/c1^2/c0^2;
ZB3=g3(0,zetasn,1,zeta1p,2,zeta0rp,1)/cs^1/c1^2/c0^1;
ZB4=g3(1,zetasn,1,zeta1p,2,zeta0rp,1)/cs^1/c1^2/c0^1;
ZB5=g3(2,zetasn,1,zeta1p,2,zeta0rp,1)/cs^1/c1^2/c0^1;
ZB6=g3(0,zetasn,1,zeta1p,1,zeta0rp,3)/cs^1/c1^1/c0^3;
ZB7=g3(1,zetasn,1,zeta1p,1,zeta0rp,3)/cs^1/c1^1/c0^3;
ZB8=g3(0,zetasn,1,zeta1p,1,zeta0rp,2)/cs^1/c1^1/c0^2;
ZB9=g3(1,zetasn,1,zeta1p,1,zeta0rp,2)/cs^1/c1^1/c0^2;
ZB10=g3(2,zetasn,1,zeta1p,1,zeta0rp,2)/cs^1/c1^1/c0^2;
ZB11=g3(0,zetasn,1,zeta1p,1,zeta0rp,1)/cs^1/c1^1/c0^1;
ZB12=g3(1,zetasn,1,zeta1p,1,zeta0rp,1)/cs^1/c1^1/c0^1;
ZB13=g3(2,zetasn,1,zeta1p,1,zeta0rp,1)/cs^1/c1^1/c0^1;
ZB14=g3(3,zetasn,1,zeta1p,1,zeta0rp,1)/cs^1/c1^1/c0^1;

    function y=fx(x)
        jsr0 = besselj(r,as*x);
        jsr1 = besselj(r-1,as*x);
        jsr2 = besselj(r+1,as*x);
        jsr3 = besselj(r-2,as*x);
        jsr4 = besselj(r+2,as*x);
        
        j0rp0 = besselj(rp,a0*x);
        j0rp1 = besselj(rp-1,a0*x);
        j0rp2 = besselj(rp+1,a0*x);
        j0rp3 = besselj(rp-2,a0*x);
        j0rp4 = besselj(rp+2,a0*x);
        
        j0np0 = besselj(np,a0*x);
        j0np1 = besselj(np-1,a0*x);
        j0np2 = besselj(np+1,a0*x);
        
        jsn0 = besselj(n,as*x);
        
        YA1 = (1/4).*c0.^2.*c1.*cs.*j0np0.*j0rp0.*jsn0.*jsr0.*r.*wc.*x;
        YA2 = (1/4).*c0.^2.*c1.*cs.^2.*j0np0.*j0rp0.*jsn0.*jsr0.*x;
        YA3 = (1/8).*c0.*c1.*j0np0.*jsn0.*(r.*(a0.*(j0rp1-j0rp2).*jsr0.*r+as.*j0rp0.*(jsr1-jsr2).*rp).*wc.^2+2.*j0rp0.*jsr0.*(c0.*cs+r.*(1i.*a0.*asy+(-2).*rp).*wc.^2).*x);
        YA4 = (1/8).*c0.*c1.*j0np0.*jsn0.*wc.*(as.*cs.*j0rp0.*(jsr1-jsr2).*rp+(-4).*j0rp0.*jsr0.*(c0.*r+cs.*rp).*x+a0.*cs.*jsr0.*(j0rp1.*r-j0rp2.*r+2i.*asy.*j0rp0.*x));
        YA5 = (-1/2).*c0.^2.*c1.*cs.*j0np0.*j0rp0.*jsn0.*jsr0.*x;
        YA6 = (1/2).*c0.^2.*cs.^2.*j0np0.*j0rp0.*jsn0.*jsr0.*r.*wc.*x;
        YA7 = (1/2).*c0.^2.*cs.^3.*j0np0.*j0rp0.*jsn0.*jsr0.*x;
        YA8 = (1/8).*c0.*cs.*jsn0.*(a0.*jsr0.*r.*wc.^2.*((-j0np1+j0np2).*j0rp0.*p+j0np0.*(j0rp1-j0rp2).*(np+r)+2i.*(a1y+asy).*j0np0.*j0rp0.*x)+j0np0.*j0rp0.*(4.*c0.*cs.*jsr0.*x+r.*(np+rp).*wc.^2.*(as.*(jsr1-jsr2)+(-4).*jsr0.*x)));
        YA9 = (1/8).*c0.*cs.*jsn0.*wc.*(a0.*cs.*jsr0.*((-j0np1+j0np2).*j0rp0.*p+j0np0.*(j0rp1-j0rp2).*(np+r)+2i.*(a1y+asy).*j0np0.*j0rp0.*x)+j0np0.*j0rp0.*(as.*cs.*(jsr1-jsr2).*(np+rp)+(-4).*jsr0.*(2.*c0.*r+cs.*(np+rp)).*x));
        YA10 = (-1).*c0.^2.*cs.^2.*j0np0.*j0rp0.*jsn0.*jsr0.*x;
        YA11 = (1/16).*jsn0.*wc.*x.^(-2).*(a0.^2.*jsr0.*r.*wc.^2.*x.*(-(j0np1-j0np2).*p.*(j0rp1.*r-j0rp2.*r+2i.*asy.*j0rp0.*x)+j0np0.*(((-2).*j0rp0+j0rp3+j0rp4).*np.*r+2i.*(j0rp1-j0rp2).*(asy.*np+a1y.*r).*x+(-4).*a1y.*asy.*j0rp0.*x.^2))+j0np0.*j0rp0.*(as.^2.*((-2).*jsr0+jsr3+jsr4).*np.*r.*rp.*wc.^2.*x+(-8).*jsr0.*(c0.^2.*r+c0.*cs.*(np+rp)+(-2).*np.*r.*rp.*wc.^2).*x.^3+(-2).*as.*(jsr1-jsr2).*(np.*r.*rp.*wc.^2-(c0.*cs.*(np+rp)+(-4).*np.*r.*rp.*wc.^2).*x.^2))+a0.*(j0np0.*(2.*(-j0rp1+j0rp2).*jsr0.*np.*r.^2.*wc.^2+as.*(j0rp1-j0rp2).*(jsr1-jsr2).*np.*r.*(r+rp).*wc.^2.*x+2.*(c0.*cs.*(j0rp1-j0rp2).*jsr0.*(np+r)+r.*((-2).*(j0rp1-j0rp2).*jsr0.*np.*(r+rp)+1i.*as.*j0rp0.*(jsr1-jsr2).*(asy.*np+a1y.*rp)).*wc.^2).*x.^2+4i.*j0rp0.*jsr0.*((a1y+asy).*c0.*cs+(-2).*r.*(asy.*np+a1y.*rp).*wc.^2).*x.^3)+(j0np1-j0np2).*j0rp0.*p.*x.*((-2).*c0.*cs.*jsr0.*x+r.*rp.*wc.^2.*(as.*(-jsr1+jsr2)+4.*jsr0.*x))));
        YA12 = (1/16).*jsn0.*x.^(-2).*(a0.^2.*cs.*jsr0.*wc.^2.*x.*(-(j0np1-j0np2).*p.*(j0rp1.*r-j0rp2.*r+2i.*asy.*j0rp0.*x)+j0np0.*(((-2).*j0rp0+j0rp3+j0rp4).*np.*r+2i.*(j0rp1-j0rp2).*(asy.*np+a1y.*r).*x+(-4).*a1y.*asy.*j0rp0.*x.^2))+j0np0.*j0rp0.*(as.^2.*cs.*((-2).*jsr0+jsr3+jsr4).*np.*rp.*wc.^2.*x+8.*jsr0.*((-3).*c0.^2.*cs+2.*(cs.*np.*rp+c0.*r.*(np+rp)).*wc.^2).*x.^3+(-2).*as.*(jsr1-jsr2).*wc.^2.*(2.*c0.*r.*(np+rp).*x.^2+cs.*np.*rp.*(1+4.*x.^2)))-a0.*wc.^2.*(4.*c0.*jsr0.*r.*x.^2.*((-j0np1+j0np2).*j0rp0.*p+j0np0.*(j0rp1-j0rp2).*(np+r)+2i.*(a1y+asy).*j0np0.*j0rp0.*x)+cs.*((j0np1-j0np2).*j0rp0.*p.*rp.*x.*(as.*(jsr1-jsr2)+(-4).*jsr0.*x)+j0np0.*(2.*(j0rp1-j0rp2).*jsr0.*np.*r-as.*(j0rp1-j0rp2).*(jsr1-jsr2).*np.*(r+rp).*x+2.*(2.*(j0rp1-j0rp2).*jsr0.*np.*(r+rp)+(-1i).*as.*j0rp0.*(jsr1-jsr2).*(asy.*np+a1y.*rp)).*x.^2+(1i*8).*j0rp0.*jsr0.*(asy.*np+a1y.*rp).*x.^3))));
        YA13 = (1/4).*c0.*jsn0.*wc.*(a0.*cs.*jsr0.*((j0np1-j0np2).*j0rp0.*p-j0np0.*(j0rp1-j0rp2).*(np+r)+(-2i).*(a1y+asy).*j0np0.*j0rp0.*x)+j0np0.*j0rp0.*(-as.*cs.*(jsr1-jsr2).*(np+rp)+4.*jsr0.*(c0.*r+cs.*(np+rp)).*x));
        YA14 = c0.^2.*cs.*j0np0.*j0rp0.*jsn0.*jsr0.*x;
        YB1 = (-1/4).*c0.^2.*c1.*cs.*j0np0.*j0rp0.*jsn0.*jsr0.*rp.*wc.*x;
        YB2 = (-1/4).*c0.^3.*c1.*cs.*j0np0.*j0rp0.*jsn0.*jsr0.*x;
        YB3 = (1/8).*c0.*c1.*j0np0.*jsn0.*(rp.*(a0.*(-j0rp1+j0rp2).*jsr0.*r+as.*j0rp0.*(-jsr1+jsr2).*rp).*wc.^2+(-2).*j0rp0.*jsr0.*(c0.*cs+(1i.*a0.*asy+(-2).*r).*rp.*wc.^2).*x);
        YB4 = (1/8).*c0.*c1.*j0np0.*jsn0.*wc.*(as.*c0.*j0rp0.*(-jsr1+jsr2).*rp+4.*j0rp0.*jsr0.*(c0.*r+cs.*rp).*x+a0.*c0.*jsr0.*(-j0rp1.*r+j0rp2.*r+(-2i).*asy.*j0rp0.*x));
        YB5 = (1/2).*c0.^2.*c1.*cs.*j0np0.*j0rp0.*jsn0.*jsr0.*x;
        YB6 = (-1/2).*c0.^3.*cs.*j0np0.*j0rp0.*jsn0.*jsr0.*rp.*wc.*x;
        YB7 = (-1/2).*c0.^4.*cs.*j0np0.*j0rp0.*jsn0.*jsr0.*x;
        YB8 = (1/8).*c0.*jsn0.*(rp.*(a0.*jsr0.*(cs.*j0np0.*(-j0rp1+j0rp2).*np+cs.*(j0np1-j0np2).*j0rp0.*p+c0.*j0np0.*(-j0rp1+j0rp2).*r)-as.*j0np0.*j0rp0.*(jsr1-jsr2).*(cs.*np+c0.*rp)).*wc.^2+(-2).*j0np0.*j0rp0.*jsr0.*(2.*c0.^2.*cs+1i.*(a0.*(asy.*c0+a1y.*cs)+2i.*(cs.*np+c0.*r)).*rp.*wc.^2).*x);
        YB9 = (-1/8).*c0.^2.*jsn0.*wc.*(j0np0.*j0rp0.*(as.*(jsr1-jsr2).*(cs.*np+c0.*rp)+(-4).*jsr0.*(cs.*np+c0.*r+2.*cs.*rp).*x)+a0.*jsr0.*(cs.*(-j0np1+j0np2).*j0rp0.*p+cs.*j0np0.*(j0rp1.*np-j0rp2.*np+2i.*a1y.*j0rp0.*x)+c0.*j0np0.*(j0rp1.*r-j0rp2.*r+2i.*asy.*j0rp0.*x)));
        YB10 = c0.^3.*cs.*j0np0.*j0rp0.*jsn0.*jsr0.*x;
        YB11 = (1/16).*jsn0.*wc.*x.^(-2).*(a0.^2.*jsr0.*rp.*wc.^2.*x.*((j0np1-j0np2).*p.*(j0rp1.*r-j0rp2.*r+2i.*asy.*j0rp0.*x)+j0np0.*(-((-2).*j0rp0+j0rp3+j0rp4).*np.*r+(-2i).*(j0rp1-j0rp2).*(asy.*np+a1y.*r).*x+4.*a1y.*asy.*j0rp0.*x.^2))+j0np0.*j0rp0.*(as.^2.*(2.*jsr0-jsr3-jsr4).*np.*rp.^2.*wc.^2.*x+8.*jsr0.*(c0.^2.*r+c0.*cs.*(np+rp)+(-2).*np.*r.*rp.*wc.^2).*x.^3+2.*as.*(jsr1-jsr2).*(np.*rp.^2.*wc.^2-(c0.*cs.*np+c0.^2.*rp+(-2).*np.*rp.*(r+rp).*wc.^2).*x.^2))+a0.*(j0np0.*(2.*(j0rp1-j0rp2).*jsr0.*np.*r.*rp.*wc.^2-as.*(j0rp1-j0rp2).*(jsr1-jsr2).*np.*rp.*(r+rp).*wc.^2.*x+2.*(-c0.*(j0rp1-j0rp2).*jsr0.*(cs.*np+c0.*r)+(-1i).*rp.*(4i.*(j0rp1-j0rp2).*jsr0.*np.*r+as.*j0rp0.*(jsr1-jsr2).*(asy.*np+a1y.*rp)).*wc.^2).*x.^2+(-4i).*j0rp0.*jsr0.*(c0.*(asy.*c0+a1y.*cs)+(-2).*(asy.*np+a1y.*r).*rp.*wc.^2).*x.^3)+(j0np1-j0np2).*j0rp0.*p.*x.*(2.*c0.*cs.*jsr0.*x+rp.*wc.^2.*(as.*(jsr1-jsr2).*rp+(-4).*jsr0.*r.*x))));
        YB12 = (1/16).*jsn0.*x.^(-2).*(a0.*wc.^2.*(2.*c0.*j0np0.*(j0rp1-j0rp2).*jsr0.*np.*r-as.*c0.*(jsr1-jsr2).*((-j0np1+j0np2).*j0rp0.*p.*rp+j0np0.*(j0rp1-j0rp2).*np.*(r+rp)).*x+((-2i).*as.*c0.*j0np0.*j0rp0.*(jsr1-jsr2).*(asy.*np+a1y.*rp)+4.*jsr0.*(cs.*(j0np0.*(j0rp1-j0rp2).*np+(-j0np1+j0np2).*j0rp0.*p).*rp+c0.*r.*((-j0np1+j0np2).*j0rp0.*p+j0np0.*(j0rp1-j0rp2).*(2.*np+rp)))).*x.^2+(1i*8).*j0np0.*j0rp0.*jsr0.*(asy.*c0.*(np+rp)+a1y.*(c0.*r+cs.*rp)).*x.^3)+a0.^2.*c0.*jsr0.*wc.^2.*x.*((j0np1-j0np2).*p.*(j0rp1.*r-j0rp2.*r+2i.*asy.*j0rp0.*x)+j0np0.*(-((-2).*j0rp0+j0rp3+j0rp4).*np.*r+(-2i).*(j0rp1-j0rp2).*(asy.*np+a1y.*r).*x+4.*a1y.*asy.*j0rp0.*x.^2))+j0np0.*j0rp0.*(as.^2.*c0.*(2.*jsr0-jsr3-jsr4).*np.*rp.*wc.^2.*x+8.*jsr0.*(3.*c0.^2.*cs+(-2).*(cs.*np.*rp+c0.*r.*(np+rp)).*wc.^2).*x.^3+2.*as.*(jsr1-jsr2).*wc.^2.*(c0.*np.*rp+2.*(c0.*np.*r+(c0+cs).*np.*rp+c0.*rp.^2).*x.^2)));
        YB13 = (1/4).*c0.*jsn0.*wc.*(j0np0.*j0rp0.*(as.*(jsr1-jsr2).*(cs.*np+c0.*rp)+(-4).*jsr0.*(c0.*r+cs.*(np+rp)).*x)+a0.*jsr0.*(cs.*(-j0np1+j0np2).*j0rp0.*p+cs.*j0np0.*(j0rp1.*np-j0rp2.*np+2i.*a1y.*j0rp0.*x)+c0.*j0np0.*(j0rp1.*r-j0rp2.*r+2i.*asy.*j0rp0.*x)));
        YB14 = (-1).*c0.^2.*cs.*j0np0.*j0rp0.*jsn0.*jsr0.*x;
        
        y = exp(-x.^2).*(YA1*ZA1+YA2*ZA2+YA3*ZA3+YA4*ZA4+YA5*ZA5+YA6*ZA6+YA7*ZA7+YA8*ZA8+YA9*ZA9+YA10*ZA10+YA11*ZA11+YA12*ZA12+YA13*ZA13+YA14*ZA14+YB1*ZB1+YB2*ZB2+YB3*ZB3+YB4*ZB4+YB5*ZB5+YB6*ZB6+YB7*ZB7+YB8*ZB8+YB9*ZB9+YB10*ZB10+YB11*ZB11+YB12*ZB12+YB13*ZB13+YB14*ZB14);
    end

in = integral(@(x)fx(x),0,inf);
ns = sqrt(nsp^2+nsz^2);
res = in*yita0^2*(wp/ns/vt)^2*exp(1i*de);
end

