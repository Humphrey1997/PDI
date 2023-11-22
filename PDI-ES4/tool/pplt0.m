EBW

num = 2001;
n0s = linspace(9,11,num);
d0s = zeros(1,num);
for i=1:num
    n0p = n0s(i);
    d0s(i) = abs(1/(1+KaiL0(w0)));
end
plot(n0s,d0s,'linewidth',2)