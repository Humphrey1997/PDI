%上边频波的包含介电张量的线性项
function y = D2(w)
global n2p n2z de2;
global w0

ws = w;
w2 = ws+w0;
n2x = n2p*cos(de2);
n2y = n2p*sin(de2);

n2Mat = [n2x n2y n2z];
n2n2Mat = n2Mat' * n2Mat;
n2squ = n2x^2+n2y^2+n2z^2;

epsilon = eye(3)+SigL1(w2);
% 电磁
y = epsilon+(n2n2Mat-n2squ.*eye(3))./w2.^2;