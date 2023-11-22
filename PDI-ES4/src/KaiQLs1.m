function kai = KaiQLs1(ws)
num = 3;
kai = 0;
for n=-num:num
    for r=-num:num
        if(abs(n-r)>num)
            continue
        end
        kai = kai + fQLs1(ws,n,r);
    end
end
end