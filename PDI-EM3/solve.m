%实部虚部扫描 寻找初始解 QL
UHW

tic
nre = 11;
nim = 11;
are = 8;
bre = 9;
aim = 0.01;
bim = 5.01;
dotre = linspace(are,bre,nre);
dotim = linspace(aim,bim,nim);
conds = zeros(nre,nim);
for i=1:nre
    for j=1:nim
        [i j]
        ww = [dotre(i) dotim(j)];
        conds(i,j) = -fw(ww);
    end
end
[X,Y]=meshgrid(dotim,dotre);
surf(X,Y,abs(conds));
xlabel('Im(\omega)');
ylabel('Re(\omega)');
zlabel('Cond number');
title([num2str(nsp),'-',num2str(de1/pi*180)]);
t = toc

