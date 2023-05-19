%一维扫描，扫描变量nsp
iv_w = mid_w;
iv_g = mid_g;
default_rad_w = center_rad_w;
default_rad_g = center_rad_g;

i = mid_i;
nsp = nP(i);
para_update_s
single_solve
data_w(i) = res_w;
data_g(i) = res_g;
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
    para_update_s
    single_solve
    data_w(i) = res_w;
    data_g(i) = res_g;
    
    gap_w = res_w-pre_w;
    gap_g = res_g-pre_g;
    iv_w = res_w;
    iv_g = res_g;
    default_rad_w = max(0.5*abs(gap_w),0.1);
    default_rad_g = max(2*abs(gap_g),0.5);
    pre_w = res_w;
    pre_g = res_g;
    if(res_g<0.02)
        left_i = i+1;
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
    para_update_s
    single_solve
    data_w(i) = res_w;
    data_g(i) = res_g;
    
    gap_w = res_w-pre_w;
    gap_g = res_g-pre_g;
    iv_w = res_w;
    iv_g = res_g;
    default_rad_w = max(0.5*abs(gap_w),0.1);
    default_rad_g = max(2*abs(gap_g),0.5);
    pre_w = res_w;
    pre_g = res_g;
    if(res_g<0.02)
        right_i = i-1;
        break
    end
end

