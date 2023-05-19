%单点计算 并行计算 QN代码
delete(gcp('nocreate'));
MyPar = parpool;

UHW
iva = 9.2;
ivb = 3;

lb = [iva-0.2,ivb-0.5];
ub = [iva+0.2,ivb+0.5];
wiv = [iva,ivb];
options = optimoptions('patternsearch','TolMesh',1e-3,'Display','iter');
tic
[wx,fval,exitflag] = patternsearch(@fw_par,wiv,[],[],[],[],lb,ub,[],options);
t = toc
wx

delete(MyPar);