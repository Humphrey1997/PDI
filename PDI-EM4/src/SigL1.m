function y = SigL1(w)
%下变频 线性
global wpe wce ve n1p n1z de1;
w1=w;
num = 3;
SigLxx=0;
SigLxy=0;
SigLxz=0;
SigLyx=0;
SigLyy=0;
SigLyz=0;
SigLzx=0;
SigLzy=0;
SigLzz=0;
A = n1p*ve/wce;
B = (n1p.*ve./wce).^2./2;
C = n1z*ve;
Sin = sin(de1);
Cos = cos(de1);
for n=-num:num
    zeta = (w1-n*wce)/n1z/ve;
    z = Z(zeta);
    
    %I0=In*exp(-b)
    I0 = besseli(n,B,1);
    %Id=d[In]/db*exp(-b)
    Id = (besseli(n-1,B,1)/2+besseli(n+1,B,1)/2);
    %I1=d[In*exp(-b)]/db
    I1 = (besseli(n-1,B,1)/2+besseli(n+1,B,1)/2-besseli(n,B,1));
    %I2=d2[In*exp(-b)]/db2
    I2 = (besseli(n-2,B,1)/4-besseli(n-1,B,1)+besseli(n,B,1)*3/2-besseli(n+1,B,1)+besseli(n+2,B,1)/4);
    
    SigLxx = SigLxx + (n^2/B*Cos^2*I0+Sin^2*(B*I2+Id))*z/C;
    SigLxy = SigLxy - 1i*(n*I1+1i*Cos*Sin*(n^2/B*I0-B*I2-Id))*z/C;
    SigLxz = SigLxz + 2*(n/A*Cos*I0-1i*A/2*Sin*I1)/C*(1+zeta*z);
    SigLyx = SigLyx + 1i*(n*I1-1i*Sin*Cos*(n^2/B*I0-B*I2-Id))*z/C;
    SigLyy = SigLyy + (n^2/B*Sin^2*I0+Cos^2*(B*I2+Id))*z/C;
    SigLyz = SigLyz + 2*(n/A*Sin*I0+1i*A/2*Cos*I1)/C*(1+zeta*z);
    SigLzx = SigLzx + 2*(n/A*Cos*I0+1i*A/2*Sin*I1)/C*(1+zeta*z);
    SigLzy = SigLzy + 2*(n/A*Sin*I0-1i*A/2*Cos*I1)/C*(1+zeta*z);
    SigLzz = SigLzz + 2*I0*zeta/C*(1+zeta*z);
end

y = wpe.^2./w1.*[SigLxx SigLxy SigLxz
    SigLyx SigLyy SigLyz
    SigLzx SigLzy SigLzz];