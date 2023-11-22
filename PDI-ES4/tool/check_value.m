EBW

data = importdata('./res/ref.txt');
[M,N] = size(data);

ks = linspace(1,M,M);
wre = zeros(1,M);
gamma = zeros(1,M);

K_se = zeros(1,M);
K_si = zeros(1,M);
S_ns = zeros(1,M);
S_eps = zeros(1,M);
S_ep1 = zeros(1,M);
S_QL1 = zeros(1,M);
S_QLs = zeros(1,M);
S_NL1 = zeros(1,M);
S_NLs = zeros(1,M);
S_NNL1 = zeros(1,M);
S_NNLs = zeros(1,M);

miu = zeros(1,M);
rat = zeros(1,M);

for k=1:M
    k
    wre(k) = data(k,2)*1e-4;
    gamma(k) = data(k,3)*1e-6;
    ws = wre(k)+gamma(k)*1i;
    nsp = data(k,1);
    des = de1-asin(n0p/nsp*sin(de1));
    n1p = sqrt(nsp^2+n0p^2-2*nsp*n0p*cos(des));
    
    S_ns(k) = nsp;
    K_se(k) = KaiLs(ws);
    K_si(k) = KaiLsi(ws);
    S_eps(k) = Eps(ws);
    S_ep1(k) = Ep1(ws);
    S_QL1(k) = KaiQL1(ws);
    S_QLs(k) = KaiQLs(ws);
    S_NL1(k) = KaiNL1(ws);
    S_NLs(k) = KaiNLs(ws);
    S_NNL1(k) = KaiNNL1(ws);
    S_NNLs(k) = KaiNNLs(ws);
    
    miu(k) = S_QL1(k)*S_QLs(k);
    rat(k) = miu(k)/S_eps(k)/(S_ep1(k)+S_NL1(k));
end
