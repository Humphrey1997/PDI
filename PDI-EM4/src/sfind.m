%%单点计算 QL代码
UHW
iva = 8.8;
ivb = 0.5;

lb = [iva-0.1,ivb-0.5];
ub = [iva+0.1,ivb+0.5];
wiv = [iva,ivb];
options = optimoptions('patternsearch','TolMesh',1e-4,'Display','iter');
tic
[wx,fval,exitflag] = patternsearch(@fw,wiv,[],[],[],[],lb,ub,[],options);
t = toc
wx