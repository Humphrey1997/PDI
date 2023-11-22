%低频模的包含介电张量的线性项
function y = Ds(w)
global nsp nsz des;

ws = w;
nsx = nsp*cos(des);
nsy = nsp*sin(des);

nMat = [nsx nsy nsz];
nnMat = nMat' * nMat;
ns2 = nsx^2+nsy^2+nsz^2;

epsilon = eye(3)+SigLs(ws)+SigLsi(ws);
% 电磁
y = epsilon+(nnMat-ns2.*eye(3))./ws.^2;