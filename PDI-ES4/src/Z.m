function [zs] = Z(zeta)
%Z 计算等离子体色散函数
%   此处显示详细说明
tm = 15;
N = 30;
a0 = 2*sqrt(pi)/tm;
z = 0;
zeta0 = zeta;
if imag(zeta)<0
    zeta = conj(zeta);
    p = 1;
else
    p = 0;
end
for n=0:N
    an = a0*exp(-n^2*pi^2/tm^2);
    z = z+an*tm*((1-exp(1i*(n*pi+tm*zeta)))/(n*pi+tm*zeta)-(1-exp(1i*(-n*pi+tm*zeta)))/(n*pi-tm*zeta));
end
z = z-a0*(1-exp(1i*tm*zeta))/zeta;
z = -z/2;
if p==0
    zs = z;
else
    zs = conj(z)+1i*2*sqrt(pi)*exp(-zeta0*zeta0);
end
end

