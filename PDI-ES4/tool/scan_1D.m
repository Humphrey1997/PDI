%一维扫描，扫描变量nsp
fclose all;
warning off;
T1 = clock;
EBW
f_log = fopen('./res/log.txt','w');

%%
%参数设置，从点center开始
Del = 180;
center_np = 0.135;
center_nz = 4e-3;
center_w = 4.8;
center_g = 3.0;
center_rad_w = 0.2;
center_rad_g = 0.2;
%网格划分
nz_start = center_nz;
nz_end = center_nz;
nz_gap = center_nz;
N = 1;
nZ = linspace(nz_start,nz_end,N);
np_start = 1e-3;
np_end = 0.5;
np_gap = 1e-3;
M = round((np_end-np_start)/(np_gap))+1;
nP = linspace(np_start,np_end,M);

data_w = zeros(M,N);
data_g = zeros(M,N);
flag = true;

center_i = round((center_np-np_start)/np_gap)+1;
center_j = round((center_nz-nz_start)/nz_gap)+1;

%%
%正式扫描
mid_i = center_i;
mid_w = center_w;
mid_g = center_g;
j = 1;
left_i = 0;
right_i = 0;
scan_p

%%
%结果输出
f_wg = fopen('./res/res_wg.txt','w');
for i=1:M
    nsp = nP(i);
    if(data_g(i)>1e-3)&&(data_w(i)>1e-3)
        para_update_s
        fprintf(f_wg,'%.4f\t%.4f\t%.4f\t\r\n',nsp,data_w(i,1),data_g(i,1));
    end
end
fclose(f_wg);

T2 = clock;
disp(['Start:',num2str(T1(1)),'/',num2str(T1(2)),'/',num2str(T1(3)),' ',num2str(T1(4)),':',num2str(T1(5))])
disp(['End:',num2str(T2(1)),'/',num2str(T2(2)),'/',num2str(T2(3)),' ',num2str(T2(4)),':',num2str(T2(5))])
fclose(f_log);