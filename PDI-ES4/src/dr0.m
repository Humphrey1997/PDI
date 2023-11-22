EBW
f = fopen('../res/DR-ES.txt','w');
numw = 101;
w0s = linspace(0.99,1,numw);

n0c = 10;

for iw = 1:numw
    num = 2001;
    n0s = linspace(n0c-1,n0c+1,num);
    d0s = zeros(1,num);
    w0 = w0s(iw)
    for i=1:num
        n0p = n0s(i);
        d0s(i) = abs(1/(1+KaiL0(w0)));
    end
    xs = find(diff(sign(diff(d0s)))==-2)+1;
    n0t = n0s(xs(1));
    n0c = round(n0t);
    fprintf(f,'%.3f\t%.3f\r\n',w0*105,n0t);
end

fclose(f);