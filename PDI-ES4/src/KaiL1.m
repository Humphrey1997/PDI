function kai = KaiL1(ws)
global wpe wce ve n1p n1z w0;

num = 3;
w1 = ws-w0;
s = 0;
A = n1p*ve/wce;
B = (n1p*ve/wce)^2/2;

for n=-num:num
    ksin = (w1-n*wce)/n1z/ve;
    In = besseli(n,B,1);
    s = s + In*Z(ksin);
end

ksi0 = w1/n1z/ve;
n1 = sqrt(n1p^2+n1z^2);
kai = 2*(wpe/ve/n1)^2*(1+ksi0*s);
end

