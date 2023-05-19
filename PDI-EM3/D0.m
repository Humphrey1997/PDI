function y = D0(n0)
%D0 泵浦色散
%   此处显示详细说明
global n0p de0 n0z
global w0

n0p = n0;

n0x = n0*cos(de0);
n0y = n0*sin(de0);

n0Mat = [n0x n0y n0z];
n0n0Mat = n0Mat' * n0Mat;
n0squ = n0x^2+n0y^2+n0z^2;

epsilon = eye(3)+SigL0(w0);

% 电磁
y = epsilon+(n0n0Mat-n0squ.*eye(3))./w0.^2;
end

