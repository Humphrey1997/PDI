function kai = KaiNLs(ws)
num = 2;
kai = 0;
for n=-num:num
    for p=-num:num
        for r=-num:num
            if(abs(n-p)>num||abs(p-r)>num)
                continue
            end
            kai = kai + fNLs(ws,n,p,r);
        end
    end
end

end