function kai = KaiL2(ws)
global wpe wce ve n2p n2z w0;

num = 3;
w2 = ws+w0;
s = 0;
A = n2p*ve/wce;
B = (n2p*ve/wce)^2/2;

for n=-num:num
    ksin = (w2-n*wce)/n2z/ve;
    In = besseli(n,B,1);
    s = s + In*Z(ksin);
end

ksi0 = w2/n2z/ve;
n2 = sqrt(n2p^2+n2z^2);
kai = 2*(wpe/ve/n2)^2*(1+ksi0*s);
end

