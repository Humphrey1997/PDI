%一维扫描，扫描变量nsp
iv_w = mid_w;
iv_g = mid_g;
default_rad_w = center_rad_w;
default_rad_g = center_rad_g;

i = mid_i;
nsp = nP(i);
nsz = nZ(j);
para_update_s
single_solve
data_w(i,j) = res_w;
data_g(i,j) = res_g;
mid_w = res_w;
mid_g = res_g;

%向左
default_rad_w = center_rad_w;
default_rad_g = center_rad_g;
pre_w = mid_w;
pre_g = mid_g;
iv_w = mid_w;
iv_g = mid_g;
sign = -1;
for i=mid_i-1:-1:1
    nsp = nP(i);
    nsz = nZ(j);
    para_update_s
    single_solve
    data_w(i,j) = res_w;
    data_g(i,j) = res_g;
    
    if(res_g<center_rad_g)
        iv_w = res_w-center_rad_w/2;
        iv_g = center_rad_g;
        default_rad_w = center_rad_w;
        default_rad_g = center_rad_g;
    else
        gap_w = res_w-pre_w;
        gap_g = res_g-pre_g;
        iv_w = res_w+gap_w;
        iv_g = res_g+gap_g;
        default_rad_w = max(0.5*abs(gap_w),0.01);
        default_rad_g = max(2*abs(gap_g),0.1);
    end
    pre_w = res_w;
    pre_g = res_g;
    if(res_w<1)
        break
    end
end

%向右
default_rad_w = center_rad_w;
default_rad_g = center_rad_g;
pre_w = mid_w;
pre_g = mid_g;
iv_w = mid_w;
iv_g = mid_g;
sign = 1;
for i=mid_i+1:1:M
    nsp = nP(i);
    nsz = nZ(j);
    para_update_s
    single_solve
    data_w(i,j) = res_w;
    data_g(i,j) = res_g;
    if(res_g<center_rad_g)
        iv_w = res_w+center_rad_w/2;
        iv_g = center_rad_g;
        default_rad_w = center_rad_w;
        default_rad_g = center_rad_g;
    else
        gap_w = res_w-pre_w;
        gap_g = res_g-pre_g;
        iv_w = res_w+gap_w;
        iv_g = res_g+gap_g;
        default_rad_w = max(0.5*abs(gap_w),0.01);
        default_rad_g = max(2*abs(gap_g),0.1);
    end
    pre_w = res_w;
    pre_g = res_g;
    if(res_w>8)
        break
    end
end

