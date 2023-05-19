%给定nsp时求n1p,用于扫描时更新参数
de1 = Del/180*pi;
des = de1-asin(n0p/nsp*sin(de1));
n1p = sqrt(nsp^2+n0p^2-2*nsp*n0p*cos(des));
n2p = sqrt(nsp^2+n0p^2+2*nsp*n0p*cos(des));
de2 = asin(nsp/n2p*sin(des));
n1z = nsz-n0z;
n2z = nsz+n0z;