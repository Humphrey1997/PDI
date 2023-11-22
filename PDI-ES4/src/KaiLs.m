function kai = KaiLs(ws)
global wpe wce ve nsp nsz;

num = 3;
s = 0;
B = (nsp*ve/wce)^2/2;

for n=-num:num
    ksin = (ws-n*wce)/nsz/ve;
    In = besseli(n,B,1);
    s = s + In*Z(ksin);
end

ksi0 = ws/nsz/ve;
ns = sqrt(nsp^2+nsz^2);
kai = 2*(wpe/ve/ns)^2*(1+ksi0*s);
end

