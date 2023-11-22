EBW

num = 1001;
wws = linspace(2.0,3.0,num);
d1s = zeros(1,num);
for i=1:num
    ws = wws(i)*1e-4;
    d1s(i) = abs(1/(1+KaiL1(ws)));
end
plot(wws,d1s,'linewidth',2)
xlabel('ws');
ylabel('|D|')
set(gca,'fontsize',30)