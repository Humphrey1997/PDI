function Sig = SigQLs1(ws)
%低频 准线性
num = 3;
Sig = complex(zeros(3,3));
for n = -num:num
    for r = -num:num
        if(abs(n-r)>num)
            continue
        end
        Sig = Sig + fQLs1(ws,n,r);
    end
end
end