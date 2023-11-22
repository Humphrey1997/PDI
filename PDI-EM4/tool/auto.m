delete(gcp('nocreate'));
MyPar = parpool;

UHW

n_start = 13;
n_gap = 0.1;
num = 31;
sign = 1;

for k=1:num
    nsp = n_start+sign*(k-1)*n_gap
    des = de1-asin(n0p/nsp*sin(de1));
    n1p = sqrt(nsp^2+n0p^2-2*nsp*n0p*cos(des));
    
    nre = 11;
    nim = 11;
    are = 6;
    bre = 10;
    aim = 0.01;
    bim = 4.01;
    dotre = linspace(are,bre,nre);
    dotim = linspace(aim,bim,nim);
    conds = zeros(nre,nim);
    for i=1:nre
        for j=1:nim
        ww = [dotre(i) dotim(j)];
        conds(i,j) = -fw_par(ww);
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

delete(MyPar);