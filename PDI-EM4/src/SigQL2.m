function Sig = SigQL2(ws)
%上边频 准线性
num = 3;
Sig = complex(zeros(3,3));
for n = -num:num
    for r = -num:num
        if(abs(n-r)>num)
            continue
        end
        Sig = Sig + fQL2(ws,n,r);
    end
end
end