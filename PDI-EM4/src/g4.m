function s = g4(p,zeta1,k1,zeta2,k2,zeta3,k3,zeta4,k4)
%g4
% 1/sqrt(pi)*integrate[u^p*exp(-u^2)/(zeta1-u)^k1/(zeta2-u)^k2/(zeta3-u)^k3/(zeta4-u)^k4,-inf,inf]
if((zeta1==zeta3)&&(zeta2==zeta4))
    s = g2(p,zeta1,k1+k3,zeta2,k2+k4);
    return
end

if(zeta1==zeta3)
    s = g3(p,zeta2,k2,zeta3,k1+k3,zeta4,k4);
    return
end

if(zeta2==zeta4)
    s = g3(p,zeta1,k1,zeta3,k3,zeta4,k2+k4);
    return
end

zetam = max([abs(zeta1),abs(zeta2),abs(zeta3),abs(zeta4)]);
if(zetam>100)
    if(zetam == abs(zeta1))
        s = g3(p,zeta2,k2,zeta3,k3,zeta4,k4)/zeta1^k1;
    elseif(zetam == abs(zeta2))
        s = g3(p,zeta1,k1,zeta3,k3,zeta4,k4)/zeta2^k2;
    elseif(zetam == abs(zeta3))
        s = g3(p,zeta1,k1,zeta2,k2,zeta4,k4)/zeta3^k3;
    else
        s = g3(p,zeta1,k1,zeta2,k2,zeta3,k3)/zeta4^k4;
    end
    return
end

if k1==1&&k2==1
    s = (g3(p,zeta1,1,zeta3,k3,zeta4,k4)-g3(p,zeta2,1,zeta3,k3,zeta4,k4))/(zeta2-zeta1);
end
if k1==1&&k2>=2
    s = (g4(p,zeta1,1,zeta2,k2-1,zeta3,k3,zeta4,k4)-g3(p,zeta2,k2,zeta3,k3,zeta4,k4))/(zeta2-zeta1);
end
if k1>=2&&k2==1
    s = (g3(p,zeta1,k1,zeta3,k3,zeta4,k4)-g4(p,zeta1,k1-1,zeta2,1,zeta3,k3,zeta4,k4))/(zeta2-zeta1);
end
if k1>=2&&k2>=2
    s = (g4(p,zeta1,k1,zeta2,k2-1,zeta3,k3,zeta4,k4)-g4(p,zeta1,k1-1,zeta2,k2,zeta3,k3,zeta4,k4))/(zeta2-zeta1);
end
end