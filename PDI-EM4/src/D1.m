%下边频波的包含介电张量的线性项
function y = D1(w)
global n1p n1z de1;
global w0

ws = w;
w1 = ws-w0;
n1x = n1p*cos(de1);
n1y = n1p*sin(de1);

n1Mat = [n1x n1y n1z];
n1n1Mat = n1Mat' * n1Mat;
n1squ = n1x^2+n1y^2+n1z^2;

epsilon = eye(3)+SigL1(w1);
% 电磁
y = epsilon+(n1n1Mat-n1squ.*eye(3))./w1.^2;