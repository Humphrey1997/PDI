UHW

data = importdata('./res/ref.txt');
[M,N] = size(data);

ks = linspace(1,M,M);
wre = zeros(1,M);
gamma = zeros(1,M);

S_ns = zeros(1,M);
S_eps = zeros(1,M);
S_ep1 = zeros(1,M);


for k=1:M
    k
    wre(k) = data(k,2)*1e-3;
    gamma(k) = data(k,3)*1e-4;
    ws = wre(k)+gamma(k)*1i;
    nsp = data(k,1);
    des = de1-asin(n0p/nsp*sin(de1));
    n1p = sqrt(nsp^2+n0p^2-2*nsp*n0p*cos(des));
    
    S_ns(k) = nsp;
    S_eps(k) = cond(Ds(ws));
    S_ep1(k) = cond(D1(ws));
    
end
