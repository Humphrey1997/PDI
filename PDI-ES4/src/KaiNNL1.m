function kai = KaiNNL1(ws)
num = 2;
kai = 0;
for n=-num:num
    for q=-num:num
        for p=-num:num
            for r=-num:num
                if((abs(n-q)>num)||(abs(q-p)>num)||(abs(p-r)>num))
                    continue
                end
                kai = kai + fNNL1(ws,n,q,p,r);
            end
        end
    end
end
end