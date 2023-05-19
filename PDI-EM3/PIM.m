%�������ȶ���ɫɢ��ϵ����
function y = PIM(w,order)
% ���
ds = Ds(w);
d1 = D1(w);
if order==0
    y = [ds, zeros(3,3);
        zeros(3,3), d1];
    return
end

sig_QL_s = SigQLs(w);
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