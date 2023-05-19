function s = g3(p,zeta1,k1,zeta2,k2,zeta3,k3)
%g3
% 1/sqrt(pi)*integrate[u^p*exp(-u^2)/(zeta1-u)^k1/(zeta2-u)^k2/(zeta3-u)^k3,-inf,inf]
s = 0;
if(abs(zeta1)>=5&&abs(zeta2)>=5&&abs(zeta3)>=5)
    fx3 = @(x) x.^p.*exp(-x.^2)./(zeta1-x).^k1./(zeta2-x).^k2./(zeta3-x).^k3;
    s = integral(fx3,-inf,inf,'RelTol',0,'AbsTol',1e-8)/sqrt(pi);
    return
else
    
if(zeta1==zeta3)
    s = g2(p,zeta1,k1+k3,zeta2,k2);
    return
end

if k1==1&&k2==1
    s = (g2(p,zeta1,1,zeta3,k3)-g2(p,zeta2,1,zeta3,k3))/(zeta2-zeta1);
end
if k1==1&&k2>=2
    s = (g3(p,zeta1,1,zeta2,k2-1,zeta3,k3)-g2(p,zeta2,k2,zeta3,k3))/(zeta2-zeta1);
end
if k1>=2&&k2==1
    s = (g2(p,zeta1,k1,zeta3,k3)-g3(p,zeta1,k1-1,zeta2,1,zeta3,k3))/(zeta2-zeta1);
end
if k1>=2&&k2>=2
    s = (g3(p,zeta1,k1,zeta2,k2-1,zeta3,k3)-g3(p,zeta1,k1-1,zeta2,k2,zeta3,k3))/(zeta2-zeta1);
end
end