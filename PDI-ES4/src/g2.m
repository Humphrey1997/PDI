function s = g2(p,zeta1,k1,zeta2,k2)
%g2
% 1/sqrt(pi)*integrate[u^p*exp(-u^2)/(zeta1-u)^k1/(zeta2-u)^k2,-inf,inf]
s = 0;

if(abs(zeta1)>=5&&abs(zeta2)>=5)
    fx2 = @(x) x.^p.*exp(-x.^2)./(zeta1-x).^k1./(zeta2-x).^k2;
    s = integral(fx2,-inf,inf,'RelTol',0,'AbsTol',1e-8)/sqrt(pi);
    return
else
    
if k1==1&&k2==1
    s = (g1(p,zeta1,1)-g1(p,zeta2,1))/(zeta2-zeta1);
end
if k1==1&&k2>=2
    s = (g2(p,zeta1,1,zeta2,k2-1)-g1(p,zeta2,k2))/(zeta2-zeta1);
end
if k1>=2&&k2==1
    s = (g1(p,zeta1,k1)-g2(p,zeta1,k1-1,zeta2,1))/(zeta2-zeta1);
end
if k1>=2&&k2>=2
    s = (g2(p,zeta1,k1,zeta2,k2-1)-g2(p,zeta1,k1-1,zeta2,k2))/(zeta2-zeta1);
end
end