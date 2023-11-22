%��Ƶģ�İ������������������
function y = Ds(w)
global nsp nsz des;

ws = w;
nsx = nsp*cos(des);
nsy = nsp*sin(des);

nMat = [nsx nsy nsz];
nnMat = nMat' * nMat;
ns2 = nsx^2+nsy^2+nsz^2;

epsilon = eye(3)+SigLs(ws)+SigLsi(ws);
% ���
y = epsilon+(nnMat-ns2.*eye(3))./ws.^2;