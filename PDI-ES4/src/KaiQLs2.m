function kai = KaiQLs2(ws)
num = 3;
kai = 0;
for n=-num:num
    for r=-num:num
        if(abs(n-r)>num)
            continue
        end
        kai = kai + fQLs2(ws,n,r);
    end
end
end