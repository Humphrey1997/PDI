UHW

num = 1001;
wws = linspace(2.0,3.0,num);
dss = zeros(1,num);
for i=1:num
    ws = wws(i)*1e-4;
    ds = Ds(ws);
    dss(i) = cond(ds);

end
plot(wws,dss,'linewidth',2)
xlabel('ws');
ylabel('|D|')
set(gca,'fontsize',30)