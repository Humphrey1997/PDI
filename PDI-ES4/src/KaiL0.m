function kai = KaiL0(w0)
global wpe wce ve n0p n0z;

num = 3;
s = 0;
B = (n0p*ve/wce)^2/2;

for n=-num:num
    ksin = (w0-n*wce)/n0z/ve;
    In = besseli(n,B,1);
    s = s + In*Z(ksin);
end

ksi0 = w0/n0z/ve;
n0 = sqrt(n0p^2+n0z^2);
kai = 2*(wpe/ve/n0)^2*(1+ksi0*s);
end
