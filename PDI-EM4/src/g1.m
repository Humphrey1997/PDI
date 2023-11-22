function s = g1(p,zeta,k)
% z = Z(zeta)=1/sqrt(pi)*integrate[exp(-u^2)/(u-zeta),-inf,inf]
% zp = 1/sqrt(pi)*integrate[exp(-u^2)/(u-zeta)^2,-inf,inf]
% zpp = 1/sqrt(pi)*integrate[exp(-u^2)/(u-zeta)^3,-inf,inf]
% zppp = 1/sqrt(pi)*integrate[exp(-u^2)/(u-zeta)^4,-inf,inf]
% 1/sqrt(pi)*integrate[u^p*exp(-u^2)/(zeta-u)^k,-inf,inf]

s = 0;

if(abs(zeta)>=5)
    fx1 = @(x) x.^p.*exp(-x.^2)./(zeta-x).^k;
    s = integral(fx1,-inf,inf,'RelTol',0,'AbsTol',1e-8)/sqrt(pi);
else
    z = Z(zeta);
    if k==1
        if p==0
            s = z;
        end
        if p==1
            s = zeta * z + 1;
        end
        if p==2
            s = zeta^2 * z + zeta;
        end
        if p==3
            s = zeta^3 * z + zeta^2 + 1/2;
        end
        if p==4
            s = zeta^4 * z + zeta^3 + 1/2 * zeta;
        end
    end
    
    if k==2
        zp = -2*(1+zeta*z);
        if p==0
            s = zp;
        end
        if p==1
            s = zeta * zp + z;
        end
        if p==2
            s = zeta^2 * zp + 2 * zeta * z + 1;
        end
        if p==3
            s = zeta^3 * zp + 3 * zeta^2 * z+ 2 * zeta;
        end
        if p==4
            s = zeta^4 * zp + 4 * zeta^3 * z+ 3 * zeta^2 + 1/2;
        end
    end
    
    if k==3
        zp = -2*(1+zeta*z);
        zpp = -(z+zeta*zp);
        if p==0
            s = zpp;
        end
        if p==1
            s = zeta * zpp + zp;
        end
        if p==2
            s = zeta^2 * zpp + 2 * zeta * zp + z;
        end
        if p==3
            s = zeta^3 * zpp + 3 * zeta^2 * zp + 3 * zeta * z + 1;
        end
        if p==4
            s = zeta^4 * zpp + 4 * zeta^3 * zp + 6 * zeta^2 * z + 3 * zeta;
        end
    end
    
    
    if k==4
        zp = -2*(1+zeta*z);
        zpp = -(z+zeta*zp);
        zppp = -2/3*(zp+zeta*zpp);
        if p==0
            s = zppp;
        end
        if p==1
            s = zeta * zppp + zpp;
        end
        if p==2
            s = zeta^2 * zppp + 2 * zeta * zpp + zp;
        end
        if p==3
            s = zeta^3 * zppp + 3 * zeta^2 * zpp + 3 * zeta * zp + z;
        end
        if p==4
            s = zeta^4 * zppp + 4 * zeta^3 * zpp + 6 * zeta^2 * zp + 4 * zeta * z + 1;
        end
    end
    
    if k==5
        zp = -2*(1+zeta*z);
        zpp = -(z+zeta*zp);
        zppp = -2/3*(zp+zeta*zpp);
        zpppp = -1/2*(zpp+zeta*zppp);
        if p==0
            s = zpppp;
        end
        if p==1
            s = zeta * zpppp + zppp;
        end
        if p==2
            s = zeta^2 * zpppp + 2 * zeta * zppp + zpp;
        end
        if p==3
            s = zeta^3 * zpppp + 3 * zeta^2 * zppp + 3 * zeta * zpp + zp;
        end
        if p==4
            s = zeta^4 * zpppp + 4 * zeta^3 * zppp + 6 * zeta^2 * zpp + 4 * zeta * zp + z;
        end
    end
    
    s = s*(-1)^k;
end
end