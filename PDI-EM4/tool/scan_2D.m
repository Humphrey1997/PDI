%%
fclose all;
warning off;
T1 = clock;
parameter
f_log = fopen('./res/log.txt','w');

%%
Del = 90;
center_np = 10.6;
center_nz = 0.003;
center_w = 4.88;
center_g = 0.1;
center_rad_w = 0.1;
center_rad_g = 0.2;

nz_start = 0.002;
nz_end = 0.03;
nz_gap = 0.002;
N = round((nz_end-nz_start)/(nz_gap))+1;
nZ = linspace(nz_start,nz_end,N);
np_start = 9.5;
np_end = 11.5;
np_gap = 0.001;
M = round((np_end-np_start)/(np_gap))+1;
nP = linspace(np_start,np_end,M);

data_w = zeros(M,N);
data_g = zeros(M,N);
flag = true;

center_i = round((center_np-np_start)/np_gap)+1;
center_j = round((center_nz-nz_start)/nz_gap)+1;

%iM(nP) jN(nZ)
%%
mid_i = center_i;
mid_w = center_w;
mid_g = center_g;
for j=center_j:1:N
    scan_p
    [max_g,max_index] = max(data_g(:,j));
    mid_i = max_index;
    mid_w = data_w(mid_i,j);
    mid_g = data_g(mid_i,j);
    if(j==center_j)
        center_i = mid_i;
    end
end

mid_i = center_i;
mid_w = data_w(center_i,center_j);
mid_g = data_g(center_i,center_j);
for j=center_j-1:-1:1
    scan_p
    [max_g,max_index] = max(data_g(:,j));
    mid_i = max_index;
    mid_w = data_w(mid_i,j);
    mid_g = data_g(mid_i,j);
end

%%
f_wk = fopen('./res/res_scan.txt','w');
for i=1:M
    for j=1:N
        fprintf(f_wk,'%.3f\t%.3f\t',data_w(i,j),data_g(i,j));
    end
    fprintf(f_wk,'\r\n');
end
fclose(f_wk);

f_np = fopen('./res/res_np.txt','w');
for i=1:M
    np = nP(i);
    para_update_s
    fprintf(f_np,'%.4f\t%.4f\t\r\n',nsp,n1p);
end
fclose(f_np);

f_nz = fopen('./res/res_nz.txt','w');
for j=1:N
    fprintf(f_nz,'%.4f\r\n',nZ(j));
end
fclose(f_nz);
%%
T2 = clock;
disp(['Start:',num2str(T1(1)),'/',num2str(T1(2)),'/',num2str(T1(3)),' ',num2str(T1(4)),':',num2str(T1(5))])
disp(['End:',num2str(T2(1)),'/',num2str(T2(2)),'/',num2str(T2(3)),' ',num2str(T2(4)),':',num2str(T2(5))])
fclose(f_log);