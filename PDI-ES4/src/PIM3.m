function mat = PIM3(ws,order)
%electron\w1
eps = 1+KaiLs(ws)+KaiLsi(ws);
ep1 = 1+KaiL1(ws);

if(order==0)
    mat=[eps,0;
        0,ep1];
    return
end

kaiQL1 = KaiQL1(ws);
kaiQLs = KaiQLs1(ws);

if(order==1)
    mat=[eps,kaiQLs;
        kaiQL1,ep1];
end

if(order==2)
    kaiNL1 = KaiNL1(ws);
    mat=[eps,kaiQLs;
        kaiQL1,ep1+kaiNL1];
end

if(order==3)
    kaiNLs = KaiNLs(ws);
    mat=[eps+kaiNLs,kaiQLs;
        kaiQL1,ep1];
end

if(order==4)
    kaiNLs = KaiNLs(ws);
    kaiNL1 = KaiNL1(ws);
    kaiNNL1 = KaiNNL1(ws);
    mat=[eps+kaiNLs,kaiQLs;
        kaiQL1+kaiNNL1,ep1+kaiNL1];
end

end