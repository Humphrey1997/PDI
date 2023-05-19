iv_w = mid_w;
iv_g = mid_g;
default_rad_w = center_rad_w;
default_rad_g = center_rad_g;

i = mid_i;
j = mid_j;
nsp = nP(i);
nsz = nZ(j);
para_update_s
single_solve
data_w(i,j) = res_w;
data_g(i,j) = res_g;
mid_w = res_w;
mid_g = res_g;

iv_w = mid_w;
iv_g = mid_g;
sign = -1;
for i=mid_i-1:-1:1
    if(data_g(i,j)>1e-3)
        if(data_w(i,j)<LW||data_g(i,j)<0.02)
            break
        else
            continue;
        end
    end
    
    nsp = nP(i);
    nsz = nZ(j);
    para_update_s
    single_solve
    data_w(i,j) = res_w;
    data_g(i,j) = res_g;
    
    iv_w = max(2*data_w(i,j)-data_w(i+1,j),center_rad_w);
    iv_g = max(2*data_g(i,j)-data_g(i+1,j),center_rad_g);
    
    if(res_w<LW)
        break
    end
    if(res_g<0.02)
        break
    end
end


iv_w = mid_w;
iv_g = mid_g;
sign = 1;
for i=mid_i+1:1:M
    if(data_g(i,j)>1e-3)
        if(data_w(i,j)>RW||data_g(i,j)<0.02)
            break
        else
            continue;
        end
    end
    
    nsp = nP(i);
    nsz = nZ(j);
    para_update_s
    single_solve
    data_w(i,j) = res_w;
    data_g(i,j) = res_g;
    
    iv_w = max(2*data_w(i,j)-data_w(i-1,j),center_rad_w);
    iv_g = max(2*data_g(i,j)-data_g(i-1,j),center_rad_g);
    
    if(res_w>RW)
        break
    end
    if(res_g<0.02)
        break
    end
end

%%
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