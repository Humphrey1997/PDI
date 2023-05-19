function Sig = SigNL1_par(ws)
%H+H
global wpe wce ve
global nsp n1p n0p
global nsz n1z n0z
global beta0 e0x e0y e0z
global des de1
global w0

num = 1;
Sig = complex(zeros(3,3));
[nL,pL,rL] = meshgrid(-num: num, -num: num, -num: num);
parfor i = 1:numel(nL)
    n = nL(i);
    p = pL(i);
    r = rL(i);
    Sig = Sig + fNL1_par(ws,n,p,r,w0,wpe,wce,ve,nsp,n1p,n0p,nsz,n1z,n0z,e0x,e0y,e0z,des,de1,beta0,num);
end

