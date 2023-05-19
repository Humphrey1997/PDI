% delete(gcp('nocreate'));
% MyPar = parpool;
profile off
profile on
sfind
% sig  = fNL1_par(wb,-1,0,-1,w0,wpe,wce,ve,nsp,n1p,n0p,nsz,n1z,n0z,e0x,e0y,e0z,des,de1,beta0,1)
% sig = SigNL1_par(wb)
% sigc = sig/sig(1,1)
profile viewer

% delete(MyPar);