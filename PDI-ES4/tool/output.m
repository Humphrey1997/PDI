f_wk = fopen('./res/res_peak.txt','w');
for i=1:M
    for j=1:N
        fprintf(f_wk,'%.3f\t%.3f\t',data_w(i,j),data_g(i,j));
    end
    fprintf(f_wk,'\r\n');
end
fclose(f_wk);

f_nsp = fopen('./res/res_nsp.txt','w');
for i=1:M
    nsp = nP(i);
    para_update_s
    fprintf(f_nsp,'%.4f\t%.4f\t\r\n',nsp,n1p);
end
fclose(f_nsp);

f_nsz = fopen('./res/res_nsz.txt','w');
for j=1:N
    fprintf(f_nsz,'%.4f\t%.4f\t',nZ(j),nZ(j));
end
fclose(f_nsz);