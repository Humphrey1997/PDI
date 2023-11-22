function mat = PIM4(ws,order)
%electron\w1+w2
eps = 1+KaiLs(ws)+KaiLsi(ws);
ep1 = 1+KaiL1(ws);
ep2 = 1+KaiL2(ws);

if(order==0)
    mat=[eps,0,0;
        0,ep1,0;
        0,0,ep2];
    return
end

kaiQL_1 = KaiQL1(ws);
kaiQL_2 = KaiQL2(ws);
kaiQL_s1 = KaiQLs1(ws);
kaiQL_s2 = KaiQLs2(ws);

%QL
if(order==1)
    mat=[eps,kaiQL_s1,kaiQL_s2;
        kaiQL_1,ep1,0;
        kaiQL_2,0,ep2];
end

%QN
if(order==2)
    kaiNL_1 = KaiNL1(ws);
    kaiNL_2 = KaiNL2(ws);
    mat=[eps,kaiQL_s1,kaiQL_s2;
        kaiQL_1,ep1+kaiNL_1,0;
        kaiQL_2,0,ep2+kaiNL_2];
end

end