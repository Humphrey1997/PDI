UHW

num = 1001;
n0s = linspace(1,11,num);
d0s = zeros(1,num);
for i=1:num
    d0 = D0(n0s(i));
    d0s(i) = cond(d0);
    
end
% plot(n0s,d0s,'linewidth',2)
% xlabel('n0');
% ylabel('|D|')
% set(gca,'fontsize',30)

xs = find(diff(sign(diff(d0s)))==-2)+1;
n0A = n0s(xs(1))
n0B = n0s(xs(2))