tic
EBW

iva = 4.8;
ivb = 2;

lb = [iva-0.1,ivb-1];
ub = [iva+0.1,ivb+1];
wiv = [iva,ivb];
options = optimoptions('patternsearch','TolMesh',1e-4,'Display','iter');
[wx,fval,exitflag] = patternsearch(@fw,wiv,[],[],[],[],lb,ub,[],options);
ww = (wx(1)+wx(2)*1i)
t=toc