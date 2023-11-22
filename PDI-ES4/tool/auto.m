% delete(gcp('nocreate'));
% MyPar = parpool;
% 
EBW

n_start = 0.07;
n_gap = 0.001;
num = 40;
sign = 1;

for k=1:num
    nsp = n_start+sign*(k-1)*n_gap;
    des = de1-asin(n0p/nsp*sin(de1));
    n1p = sqrt(nsp^2+n0p^2-2*nsp*n0p*cos(des));
    [nsp n1p]
    nre = 11;
    nim = 11;
    are = 2.3;
    bre = 2.6;
    aim = 0.01;
    bim = 1.01;
    dotre = linspace(are,bre,nre);
    dotim = linspace(aim,bim,nim);
    conds = zeros(nre,nim);
    for i=1:nre
        for j=1:nim
        ww = [dotre(i) dotim(j)];
        conds(i,j) = -fw(ww);
        end
    end
    [X,Y]=meshgrid(dotim,dotre);
    surf(X,Y,abs(conds));
    xlabel('Im(\omega)');
    ylabel('Re(\omega)');
    zlabel('Cond number');
    
    saveas(gcf,[num2str(nsp),'.png'],'png');
    close(gcf);
    clock
end

% delete(MyPar);