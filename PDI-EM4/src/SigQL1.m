function Sig = SigQL1(ws)
%下边频 准线性
num = 3;
Sig = complex(zeros(3,3));
for n = -num:num
    for r = -num:num
        if(abs(n-r)>num)
            continue
        end
        Sig = Sig + fQL1(ws,n,r);
    end
end
end