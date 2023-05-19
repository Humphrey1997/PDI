UHW

num = 1001;
wws = linspace(8.0,12.0,num);
d1s = zeros(1,num);
for i=1:num
    ws = wws(i)*1e-3;
    d1 = D1(ws);
    d1s(i) = cond(d1);

end
plot(wws,d1s,'linewidth',2)
xlabel('ws');
ylabel('|D|')
set(gca,'fontsize',30)