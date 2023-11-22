fprintf(f_log,'Start:%d\t%.4f\t%.4f\t%.4f\t%.4f\r\n',Del,nZ(mid_j),nP(mid_i),mid_w,mid_g);

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
        if(data_w(i,j)<LW||data_g(i,j)<0.1)
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
    if(res_g<0.1)
        break
    end
end


iv_w = mid_w;
iv_g = mid_g;
sign = 1;
for i=mid_i+1:1:M
    if(data_g(i,j)>1e-3)
        if(data_w(i,j)>RW||data_g(i,j)<0.1)
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
    if(res_g<0.1)
        break
    end
end

%%
output