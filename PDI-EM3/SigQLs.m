function Sig = SigQLs(ws)
%H+H
num = 3;
Sig = complex(zeros(3,3));
for n = -num:num
    for r = -num:num
        if(abs(n-r)>num)
            continue
        end
        Sig = Sig + fQLs(ws,n,r);
    end
end
end