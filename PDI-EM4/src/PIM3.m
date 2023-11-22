%参量不稳定性色散关系矩阵
%低频-泵浦-下边频(6x6)
function y = PIM3(w,order)
% 电磁
ds = Ds(w);
d1 = D1(w);
if order==0
    y = [ds, zeros(3,3);
        zeros(3,3), d1];
    return
end

sig_QL_s = SigQLs1(w);
sig_QL_1 = SigQL1(w);

if order==1
    y = [ds, sig_QL_s;
        sig_QL_1, d1];
end
if order==2
    sig_NL_1 = SigNL1_par(w);
    y = [ds, sig_QL_s;
        sig_QL_1, d1+sig_NL_1];
end
if order==3
    sig_NL_s = SigNLs(w);
    y = [ds+sig_NL_s, sig_QL_s;
        sig_QL_1, d1];
end
end