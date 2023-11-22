function res = fNL1(ws,n,p,r)
global wpe wce ve
global nsp n1p n0p
global nsz n1z n0z
global des de1 de0
global yita0 w0

w1 = ws-w0;
wc = wce;
vt = ve;
wp = wpe;

de = (n-r)*de1;
pr = p-r;
pn = p-n;

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
zeta1n = (w1-n*wc)/c1;
zetasp = (ws-p*wc)/cs;
zeta1r = (w1-r*wc)/c1;
zeta0pr = (w0-pr*wc)/c0;

ZA1=g3(0,zeta1n,1,zetasp,2,zeta1r,2)/c1^1/cs^2/c1^2;
ZA2=g3(1,zeta1n,1,zetasp,2,zeta1r,2)/c1^1/cs^2/c1^2;
ZA3=g3(0,zeta1n,1,zetasp,2,zeta1r,1)/c1^1/cs^2/c1^1;
ZA4=g3(1,zeta1n,1,zetasp,2,zeta1r,1)/c1^1/cs^2/c1^1;
ZA5=g3(2,zeta1n,1,zetasp,2,zeta1r,1)/c1^1/cs^2/c1^1;
ZA6=g3(0,zeta1n,1,zetasp,1,zeta1r,3)/c1^1/cs^1/c1^3;
ZA7=g3(1,zeta1n,1,zetasp,1,zeta1r,3)/c1^1/cs^1/c1^3;
ZA8=g3(0,zeta1n,1,zetasp,1,zeta1r,2)/c1^1/cs^1/c1^2;
ZA9=g3(1,zeta1n,1,zetasp,1,zeta1r,2)/c1^1/cs^1/c1^2;
ZA10=g3(2,zeta1n,1,zetasp,1,zeta1r,2)/c1^1/cs^1/c1^2;
ZA11=g3(0,zeta1n,1,zetasp,1,zeta1r,1)/c1^1/cs^1/c1^1;
ZA12=g3(1,zeta1n,1,zetasp,1,zeta1r,1)/c1^1/cs^1/c1^1;
ZA13=g3(2,zeta1n,1,zetasp,1,zeta1r,1)/c1^1/cs^1/c1^1;
ZA14=g3(3,zeta1n,1,zetasp,1,zeta1r,1)/c1^1/cs^1/c1^1;

ZB1=g3(0,zeta1n,1,zetasp,2,zeta0pr,2)/c1^1/cs^2/c0^2;
ZB2=g3(1,zeta1n,1,zetasp,2,zeta0pr,2)/c1^1/cs^2/c0^2;
ZB3=g3(0,zeta1n,1,zetasp,2,zeta0pr,1)/c1^1/cs^2/c0^1;
ZB4=g3(1,zeta1n,1,zetasp,2,zeta0pr,1)/c1^1/cs^2/c0^1;
ZB5=g3(2,zeta1n,1,zetasp,2,zeta0pr,1)/c1^1/cs^2/c0^1;
ZB6=g3(0,zeta1n,1,zetasp,1,zeta0pr,3)/c1^1/cs^1/c0^3;
ZB7=g3(1,zeta1n,1,zetasp,1,zeta0pr,3)/c1^1/cs^1/c0^3;
ZB8=g3(0,zeta1n,1,zetasp,1,zeta0pr,2)/c1^1/cs^1/c0^2;
ZB9=g3(1,zeta1n,1,zetasp,1,zeta0pr,2)/c1^1/cs^1/c0^2;
ZB10=g3(2,zeta1n,1,zetasp,1,zeta0pr,2)/c1^1/cs^1/c0^2;
ZB11=g3(0,zeta1n,1,zetasp,1,zeta0pr,1)/c1^1/cs^1/c0^1;
ZB12=g3(1,zeta1n,1,zetasp,1,zeta0pr,1)/c1^1/cs^1/c0^1;
ZB13=g3(2,zeta1n,1,zetasp,1,zeta0pr,1)/c1^1/cs^1/c0^1;
ZB14=g3(3,zeta1n,1,zetasp,1,zeta0pr,1)/c1^1/cs^1/c0^1;

    function y=fx(x)
        j1r0 = besselj(r,a1*x);
        j1r1 = besselj(r-1,a1*x);
        j1r2 = besselj(r+1,a1*x);
        j1r3 = besselj(r-2,a1*x);
        j1r4 = besselj(r+2,a1*x);
        
        j0pr0 = besselj(pr,a0*x);
        j0pr1 = besselj(pr-1,a0*x);
        j0pr2 = besselj(pr+1,a0*x);
        j0pr3 = besselj(pr-2,a0*x);
        j0pr4 = besselj(pr+2,a0*x);
        
        j0pn0 = besselj(pn,a0*x);
        j0pn1 = besselj(pn-1,a0*x);
        j0pn2 = besselj(pn+1,a0*x);
        
        j1n0 = besselj(n,a1*x);
        
        YA1 = (1/4).*c0.^2.*c1.*cs.*j0pn0.*j0pr0.*j1n0.*j1r0.*r.*wc.*x;
        YA2 = (1/4).*c0.^2.*c1.^2.*cs.*j0pn0.*j0pr0.*j1n0.*j1r0.*x;
        YA3 = (1/8).*c0.*cs.*j0pn0.*j1n0.*(r.*(a1.*j0pr0.*(j1r1-j1r2).*pr+a0.*(-j0pr1+j0pr2).*j1r0.*r).*wc.^2+2.*j0pr0.*j1r0.*(c0.*c1+(1i.*a0.*a1y+(-2).*pr).*r.*wc.^2).*x);
        YA4 = (1/8).*c0.*cs.*j0pn0.*j1n0.*wc.*(a1.*c1.*j0pr0.*(j1r1-j1r2).*pr+(-4).*j0pr0.*j1r0.*(c1.*pr+c0.*r).*x+a0.*c1.*j1r0.*(-j0pr1.*r+j0pr2.*r+2i.*a1y.*j0pr0.*x));
        YA5 = (-1/2).*c0.^2.*c1.*cs.*j0pn0.*j0pr0.*j1n0.*j1r0.*x;
        YA6 = (1/2).*c0.^2.*c1.^2.*j0pn0.*j0pr0.*j1n0.*j1r0.*r.*wc.*x;
        YA7 = (1/2).*c0.^2.*c1.^3.*j0pn0.*j0pr0.*j1n0.*j1r0.*x;
        YA8 = (1/8).*c0.*c1.*j1n0.*(a0.*j1r0.*r.*wc.^2.*((j0pn1-j0pn2).*j0pr0.*p+j0pn0.*(j0pr1-j0pr2).*(pn-r)+2i.*(a1y+asy).*j0pn0.*j0pr0.*x)+j0pn0.*j0pr0.*(4.*c0.*c1.*j1r0.*x+(pn+pr).*r.*wc.^2.*(a1.*(j1r1-j1r2)+(-4).*j1r0.*x)));
        YA9 = (1/8).*c0.*c1.*j1n0.*wc.*(a0.*c1.*j1r0.*((j0pn1-j0pn2).*j0pr0.*p+j0pn0.*(j0pr1-j0pr2).*(pn-r)+2i.*(a1y+asy).*j0pn0.*j0pr0.*x)+j0pn0.*j0pr0.*(a1.*c1.*(j1r1-j1r2).*(pn+pr)+(-4).*j1r0.*(c1.*(pn+pr)+2.*c0.*r).*x));
        YA10 = (-1).*c0.^2.*c1.^2.*j0pn0.*j0pr0.*j1n0.*j1r0.*x;
        YA11 = (-1/16).*j1n0.*wc.*x.^(-2).*(a1.^2.*j0pn0.*j0pr0.*(2.*j1r0-j1r3-j1r4).*pn.*pr.*r.*wc.^2.*x+a1.*(j1r1-j1r2).*(2.*j0pn0.*j0pr0.*pn.*pr.*r.*wc.^2+a0.*(-j0pn1.*j0pr0.*p.*pr+j0pn2.*j0pr0.*p.*pr-j0pn0.*(j0pr1-j0pr2).*pn.*(pr-r)).*r.*wc.^2.*x+(-2).*j0pn0.*j0pr0.*(c0.*c1.*(pn+pr)+1i.*(a0.*a1y.*pn+a0.*asy.*pr+4i.*pn.*pr).*r.*wc.^2).*x.^2)+j1r0.*(2.*a0.*j0pn0.*(-j0pr1+j0pr2).*pn.*r.^2.*wc.^2+a0.^2.*((j0pn1-j0pn2).*(j0pr1-j0pr2).*p+j0pn0.*((-2).*j0pr0+j0pr3+j0pr4).*pn).*r.^2.*wc.^2.*x+2.*a0.*(c0.*c1.*((-j0pn1+j0pn2).*j0pr0.*p-j0pn0.*(j0pr1-j0pr2).*(pn-r))+r.*(((j0pn1-j0pn2).*j0pr0.*p+j0pn0.*(j0pr1-j0pr2).*pn).*((-1i).*a0.*a1y+2.*pr)+1i.*j0pn0.*(j0pr1-j0pr2).*(a0.*asy+2i.*pn).*r).*wc.^2).*x.^2+4.*j0pn0.*j0pr0.*(c0.*((-1i).*a0.*(a1y+asy).*c1+2.*c1.*(pn+pr)+2.*c0.*r)+(a0.*asy+2i.*pn).*(a0.*a1y+2i.*pr).*r.*wc.^2).*x.^3));
        YA12 = (-1/16).*j1n0.*x.^(-2).*(a1.^2.*c1.*j0pn0.*j0pr0.*(2.*j1r0-j1r3-j1r4).*pn.*pr.*wc.^2.*x+a1.*(j1r1-j1r2).*wc.^2.*(2.*c1.*j0pn0.*j0pr0.*pn.*pr+a0.*c1.*(-j0pn1.*j0pr0.*p.*pr+j0pn2.*j0pr0.*p.*pr-j0pn0.*(j0pr1-j0pr2).*pn.*(pr-r)).*x+2.*j0pn0.*j0pr0.*(4.*c1.*pn.*pr+(-1i).*a0.*c1.*(a1y.*pn+asy.*pr)+2.*c0.*(pn+pr).*r).*x.^2)+j1r0.*(2.*a0.*c1.*j0pn0.*(-j0pr1+j0pr2).*pn.*r.*wc.^2+a0.^2.*c1.*((j0pn1-j0pn2).*(j0pr1-j0pr2).*p+j0pn0.*((-2).*j0pr0+j0pr3+j0pr4).*pn).*r.*wc.^2.*x+2.*a0.*(c1.*((j0pn1-j0pn2).*j0pr0.*p+j0pn0.*(j0pr1-j0pr2).*pn).*((-1i).*a0.*a1y+2.*pr)+(2.*c0.*(j0pn1-j0pn2).*j0pr0.*p+j0pn0.*(j0pr1-j0pr2).*(1i.*a0.*asy.*c1+2.*(c0-c1).*pn)).*r+2.*c0.*j0pn0.*(-j0pr1+j0pr2).*r.^2).*wc.^2.*x.^2+4.*j0pn0.*j0pr0.*(6.*c0.^2.*c1+(c1.*(a0.*asy+2i.*pn).*(a0.*a1y+2i.*pr)+2i.*a0.*(a1y+asy).*c0.*r+(-4).*c0.*(pn+pr).*r).*wc.^2).*x.^3));
        YA13 = (1/4).*c0.*j1n0.*wc.*(a0.*c1.*j1r0.*((-j0pn1+j0pn2).*j0pr0.*p-j0pn0.*(j0pr1-j0pr2).*(pn-r)+(-2i).*(a1y+asy).*j0pn0.*j0pr0.*x)+j0pn0.*j0pr0.*(-a1.*c1.*(j1r1-j1r2).*(pn+pr)+4.*j1r0.*(c1.*(pn+pr)+c0.*r).*x));
        YA14 = c0.^2.*c1.*j0pn0.*j0pr0.*j1n0.*j1r0.*x;
        YB1 = (1/4).*c0.^2.*c1.*cs.*j0pn0.*j0pr0.*j1n0.*j1r0.*pr.*wc.*x;
        YB2 = (1/4).*c0.^3.*c1.*cs.*j0pn0.*j0pr0.*j1n0.*j1r0.*x;
        YB3 = (1/8).*c0.*cs.*j0pn0.*j1n0.*(pr.*(a1.*j0pr0.*(-j1r1+j1r2).*pr+a0.*(j0pr1-j0pr2).*j1r0.*r).*wc.^2+2.*j0pr0.*j1r0.*(c0.*c1+pr.*((-1i).*a0.*a1y+(-2).*r).*wc.^2).*x);
        YB4 = (-1/8).*c0.*cs.*j0pn0.*j1n0.*wc.*(a1.*c0.*j0pr0.*(j1r1-j1r2).*pr+4.*j0pr0.*j1r0.*(c1.*pr+c0.*r).*x+a0.*c0.*j1r0.*(-j0pr1.*r+j0pr2.*r+2i.*a1y.*j0pr0.*x));
        YB5 = (-1/2).*c0.^2.*c1.*cs.*j0pn0.*j0pr0.*j1n0.*j1r0.*x;
        YB6 = (1/2).*c0.^3.*c1.*j0pn0.*j0pr0.*j1n0.*j1r0.*pr.*wc.*x;
        YB7 = (1/2).*c0.^4.*c1.*j0pn0.*j0pr0.*j1n0.*j1r0.*x;
        YB8 = (1/8).*c0.*j1n0.*(pr.*(a1.*j0pn0.*j0pr0.*(j1r1-j1r2).*(c1.*pn-c0.*pr)+a0.*j1r0.*(c1.*(j0pn1-j0pn2).*j0pr0.*p+c1.*j0pn0.*(j0pr1-j0pr2).*pn+c0.*j0pn0.*(j0pr1-j0pr2).*r)).*wc.^2+2.*j0pn0.*j0pr0.*j1r0.*(2.*c0.^2.*c1+pr.*((-1i).*a0.*(a1y.*c0-asy.*c1)+(-2).*(c1.*pn+c0.*r)).*wc.^2).*x);
        YB9 = (1/8).*c0.^2.*j1n0.*wc.*(j0pn0.*j0pr0.*(a1.*(j1r1-j1r2).*(c1.*pn-c0.*pr)+(-4).*j1r0.*(c1.*pn+2.*c1.*pr+c0.*r).*x)+a0.*j1r0.*(c1.*j0pn0.*(j0pr1-j0pr2).*pn+c1.*j0pr0.*(j0pn1.*p-j0pn2.*p+2i.*asy.*j0pn0.*x)+c0.*j0pn0.*(j0pr1.*r-j0pr2.*r+(-2i).*a1y.*j0pr0.*x)));
        YB10 = (-1).*c0.^3.*c1.*j0pn0.*j0pr0.*j1n0.*j1r0.*x;
        YB11 = (1/16).*j1n0.*wc.*x.^(-2).*(a1.^2.*j0pn0.*j0pr0.*(2.*j1r0-j1r3-j1r4).*pn.*pr.^2.*wc.^2.*x+a1.*(j1r1-j1r2).*(2.*j0pn0.*j0pr0.*pn.*pr.^2.*wc.^2+a0.*pr.*(-j0pn1.*j0pr0.*p.*pr+j0pn2.*j0pr0.*p.*pr-j0pn0.*(j0pr1-j0pr2).*pn.*(pr-r)).*wc.^2.*x+(-2).*j0pn0.*j0pr0.*(-c0.*c1.*pn+c0.^2.*pr+pr.*(1i.*a0.*(a1y.*pn+asy.*pr)+2.*pn.*(-pr+r)).*wc.^2).*x.^2)+j1r0.*(2.*a0.*j0pn0.*(-j0pr1+j0pr2).*pn.*pr.*r.*wc.^2+a0.^2.*((j0pn1-j0pn2).*(j0pr1-j0pr2).*p+j0pn0.*((-2).*j0pr0+j0pr3+j0pr4).*pn).*pr.*r.*wc.^2.*x+2.*a0.*(c0.*c1.*((j0pn1-j0pn2).*j0pr0.*p+j0pn0.*(j0pr1-j0pr2).*pn)+c0.^2.*j0pn0.*(j0pr1-j0pr2).*r+pr.*((-1i).*a0.*a1y.*((j0pn1-j0pn2).*j0pr0.*p+j0pn0.*(j0pr1-j0pr2).*pn)+(2.*(-j0pn1+j0pn2).*j0pr0.*p+j0pn0.*(j0pr1-j0pr2).*(1i.*a0.*asy+(-4).*pn)).*r).*wc.^2).*x.^2+4.*j0pn0.*j0pr0.*((-1i).*a0.*c0.*(a1y.*c0-asy.*c1)+(-2).*c0.*(c1.*(pn+pr)+c0.*r)+(a0.*asy+2i.*pn).*pr.*(a0.*a1y+(-2i).*r).*wc.^2).*x.^3));
        YB12 = (1/16).*j1n0.*x.^(-2).*(a1.^2.*c0.*j0pn0.*j0pr0.*(2.*j1r0-j1r3-j1r4).*pn.*pr.*wc.^2.*x+a1.*(j1r1-j1r2).*wc.^2.*(2.*c0.*j0pn0.*j0pr0.*pn.*pr+a0.*c0.*(-j0pn1.*j0pr0.*p.*pr+j0pn2.*j0pr0.*p.*pr-j0pn0.*(j0pr1-j0pr2).*pn.*(pr-r)).*x+2.*j0pn0.*j0pr0.*((-2).*c1.*pn.*pr+2.*c0.*pr.*(pn+pr)+(-1i).*a0.*c0.*(a1y.*pn+asy.*pr)+(-2).*c0.*pn.*r).*x.^2)+j1r0.*(2.*a0.*c0.*j0pn0.*(-j0pr1+j0pr2).*pn.*r.*wc.^2+a0.^2.*c0.*((j0pn1-j0pn2).*(j0pr1-j0pr2).*p+j0pn0.*((-2).*j0pr0+j0pr3+j0pr4).*pn).*r.*wc.^2.*x+2.*a0.*((-1i).*((j0pn1-j0pn2).*j0pr0.*p+j0pn0.*(j0pr1-j0pr2).*pn).*(a0.*a1y.*c0+(-2i).*c1.*pr)+c0.*(2.*(-j0pn1+j0pn2).*j0pr0.*p+1i.*j0pn0.*(j0pr1-j0pr2).*(a0.*asy+2i.*(2.*pn+pr))).*r).*wc.^2.*x.^2+4.*j0pn0.*j0pr0.*((-6).*c0.^2.*c1+2.*c1.*((-1i).*a0.*asy+2.*pn).*pr.*wc.^2+c0.*(a0.*asy+2i.*(pn+pr)).*(a0.*a1y+(-2i).*r).*wc.^2).*x.^3));
        YB13 = (1/4).*c0.*j1n0.*wc.*(j0pn0.*j0pr0.*(-a1.*(j1r1-j1r2).*(c1.*pn-c0.*pr)+4.*j1r0.*(c1.*(pn+pr)+c0.*r).*x)+a0.*j1r0.*(c1.*j0pn0.*(-j0pr1+j0pr2).*pn+c1.*j0pr0.*(-j0pn1.*p+j0pn2.*p+(-2i).*asy.*j0pn0.*x)+c0.*j0pn0.*(-j0pr1.*r+j0pr2.*r+2i.*a1y.*j0pr0.*x)));
        YB14 = c0.^2.*c1.*j0pn0.*j0pr0.*j1n0.*j1r0.*x;
        
        y = exp(-x.^2).*(YA1*ZA1+YA2*ZA2+YA3*ZA3+YA4*ZA4+YA5*ZA5+YA6*ZA6+YA7*ZA7+YA8*ZA8+YA9*ZA9+YA10*ZA10+YA11*ZA11+YA12*ZA12+YA13*ZA13+YA14*ZA14+YB1*ZB1+YB2*ZB2+YB3*ZB3+YB4*ZB4+YB5*ZB5+YB6*ZB6+YB7*ZB7+YB8*ZB8+YB9*ZB9+YB10*ZB10+YB11*ZB11+YB12*ZB12+YB13*ZB13+YB14*ZB14);
    end

in = integral(@(x)fx(x),0,inf);
n1 = sqrt(n1p^2+n1z^2);
res = in*yita0^2*(wp/n1/vt)^2*exp(1i*de);
end

