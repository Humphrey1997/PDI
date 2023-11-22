fclose all;
warning off;
T1 = clock;
UHW
f_log = fopen('./res/log.txt','w');

Del = 180;
center_np = 4.92;
center_nz = 0.015;
center_w = 9;
center_g = 1.2;
center_rad_w = 0.05;
center_rad_g = 0.2;

LW = 8.2;
RW = 10.2;

nz_start = 0;
nz_end = 0.1;
nz_gap = 5e-3;
N = round((nz_end-nz_start)/(nz_gap))+1;
nZ = linspace(nz_start,nz_end,N)+1e-6;
np_start = 4.5;
np_end = 5.5;
np_gap = 5e-3;
M = round((np_end-np_start)/(np_gap))+1;
nP = linspace(np_start,np_end,M);

data_w = zeros(M,N);
data_g = zeros(M,N);
flag = true;

center_i = round((center_np-np_start)/np_gap)+1;
center_j = round((center_nz-nz_start)/nz_gap)+1;

%iM(nP) jN(nZ)
mid_j = center_j;
mid_i = center_i;
mid_w = center_w;
mid_g = center_g;
scan_peak
[max_g,max_index] = max(data_g(:,mid_j));
center_i = max_index;
center_w = data_w(mid_i,mid_j);
center_g = data_g(mid_i,mid_j);

sign_j = 1;
max_j = N;
scan_peak2
sign_j = -1;
min_j = 1;
scan_peak2

T2 = clock;
disp(['Start:',num2str(T1(1)),'/',num2str(T1(2)),'/',num2str(T1(3)),' ',num2str(T1(4)),':',num2str(T1(5))])
disp(['End:',num2str(T2(1)),'/',num2str(T2(2)),'/',num2str(T2(3)),' ',num2str(T2(4)),':',num2str(T2(5))])

fclose(f_log);