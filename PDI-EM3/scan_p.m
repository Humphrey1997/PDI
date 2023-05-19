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
gap_w = 0.03;
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
    iv_w = res_w-gap_w;
    iv_g = res_g;
    gap_w = max(abs(iv_w-pre_w),0.01);
    pre_w = iv_w;  
    if(res_w<8)
        break
    end
end

%向右
% gap_w = 0.03;
% pre_w = mid_w;
% pre_g = mid_g;
% iv_w = mid_w;
% iv_g = mid_g;
% sign = 1;
% for i=mid_i+1:1:M
%     nsp = nP(i);
%     nsz = nZ(j);
%     para_update_s
%     single_solve
%     data_w(i,j) = res_w;
%     data_g(i,j) = res_g;
%     iv_w = res_w+gap_w;
%     iv_g = res_g;
%     gap_w = max(abs(iv_w-pre_w),0.01);
%     pre_w = iv_w;
%     if(res_w>8)
%         break
%     end
% end

