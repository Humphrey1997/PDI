%一维扫描，扫描变量nsp
fclose all;
warning off;
T1 = clock;
UHW
f_log = fopen('./res/log.txt','w');

%%
%参数设置，从点center开始
Del = 180;
center_rad_w = 0.1;
center_rad_g = 0.2;
%网格划分
np_start = 0;
np_end = 50;
np_gap = 0.01;
M = round((np_end-np_start)/(np_gap))+1;
nP = linspace(np_start,np_end,M);

data_w = zeros(M,1);
data_g = zeros(M,1);
left_i = 0;
right_i = 0;
flag = true;

%%
%正式扫描
center_np = 7;
center_w = 15.3;
center_g = 0.2;

mid_i = round((center_np-np_start)/np_gap)+1;
mid_w = center_w;
mid_g = center_g;

scan_pDD

[max_g,max_i] = max(data_g(left_i:right_i));
max_i_mid = max_i+left_i-1;
left_i_mid = left_i;
right_i_mid = right_i;

mid_i = max_i_mid;
gap_i = 7;
while(1)
    mid_i_old = mid_i;
    mid_w = data_w(mid_i_old) + vci;
    mid_g = data_g(mid_i_old);
    mid_i = mid_i_old + gap_i;
    scan_pDD
    tmp = right_i;
    while(1)
        if((data_g(tmp)>=data_g(tmp-1))&&(data_g(tmp)>=data_g(tmp+1)))
            mid_i = tmp;
            break
        end
        tmp = tmp - 1;
    end
    gap_i = min(abs(mid_i - mid_i_old),10);
    if(data_w(right_i)>25)
        break
    end
end

% mid_i = max_i_mid;
% gap_i = 9;
% while(1)
%     mid_i_old = mid_i;
%     mid_w = data_w(mid_i_old) - vci;
%     mid_g = data_g(mid_i_old);
%     mid_i = mid_i_old - gap_i;
%     scan_pDD
%     tmp = left_i;
%     while(1)
%         if((data_g(tmp)>=data_g(tmp-1))&&(data_g(tmp)>=data_g(tmp+1)))
%             mid_i = tmp;
%             break
%         end
%         tmp = tmp + 1;
%     end
%     gap_i = min(abs(mid_i - mid_i_old),10);
%     if(data_w(left_i)<5)
%         break
%     end
% end
%%
%结果输出
f_wg = fopen('./res/res_wg.txt','w');
for i=1:M
    nsp = nP(i);
    if(data_g(i)>1e-3)&&(data_w(i)>1e-3)
        fprintf(f_wg,'%.3f\t%.3f\t%.3f\t\r\n',nsp,data_w(i,1),data_g(i,1));
    end
end
fclose(f_wg);

T2 = clock;
disp(['Start:',num2str(T1(1)),'/',num2str(T1(2)),'/',num2str(T1(3)),' ',num2str(T1(4)),':',num2str(T1(5))])
disp(['End:',num2str(T2(1)),'/',num2str(T2(2)),'/',num2str(T2(3)),' ',num2str(T2(4)),':',num2str(T2(5))])
fclose(f_log);

open('./res/res_wg.txt')