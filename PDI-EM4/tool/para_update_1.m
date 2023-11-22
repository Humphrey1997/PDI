%给定n1p时求nsp,用于扫描时更新参数
%应用时一般给定nz
de1 = Del/180*pi;
nsp = sqrt(n1p^2+n0p^2+2*n1p*n0p*cos(de1));
des = de1-asin(n0p/nsp*sin(de1));
n2p = sqrt(nsp^2+n0p^2+2*nsp*n0p*cos(des));
de2 = asin(nsp/n2p*sin(des));