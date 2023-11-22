%参量不稳定性色散关系矩阵
%低频-泵浦-下边频-上边频(9x9)
function y = PIM4(w,order)
% 电磁
ds = Ds(w);
d1 = D1(w);
d2 = D2(w);
if order==0
    y = [ds, zeros(3,3),zeros(3,3);
        zeros(3,3), d1,zeros(3,3);
        zeros(3,3),zeros(3,3),d2];
    return
end

sig_QL_s1 = SigQLs1(w);
sig_QL_s2 = SigQLs2(w);
sig_QL_1 = SigQL1(w);
sig_QL_2 = SigQL2(w);

if order==1
    y = [ds, sig_QL_s1,sig_QL_s2;
        sig_QL_1, d1,zeros(3,3);
        sig_QL_2,zeros(3,3),d2];
end

end
