function kai = KaiLsi(ws)
global wpi wci vi nsp nsz;

num = 150;
s = 0;
B = (nsp*vi/wci)^2/2;

for n=-num:num
    ksin = (ws-n*wci)/nsz/vi;
    In = besseli(n,B,1);
    s = s + In*Z(ksin);
end

ksi0 = ws/nsz/vi;
ns = sqrt(nsp^2+nsz^2);
kai = 2*(wpi/vi/ns)^2*(1+ksi0*s);
end

