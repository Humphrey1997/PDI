fclose all;
warning off;
T1 = clock;
UHW
f_log = fopen('./res/log.txt','w');

Del = 180;
nsz = 1e-6;
nsp_center = 4.91;
yita_center = 3;
w_center = 9.0;
g_center = 3.7;

yita_start = 0;
yita_end = 5;
yita_gap = 0.2;
numy = round((yita_end-yita_start)/(yita_gap))+1;
nY = linspace(yita_start,yita_end,numy);

center_N = round((yita_center-yita_start)/yita_gap)+1;

S_ynF = zeros(1,numy);
S_ywF = zeros(1,numy);
S_ygF = zeros(1,numy);

default_rad_w = 0.1;
default_rad_g = 0.2;
gap_nsp = 1e-3;

for Ny = center_N:numy
    Ny
    yita0 = nY(Ny);
    beta0 = yita0*1e6*(e/me/(ve*c)/W0);
    iv_nsp = nsp_center;
    iv_w = w_center;
    iv_g = g_center;
    
    scan_max
    
    nsp_center = max_nsp;
    w_center = max_w;
    g_center = max_g;
    
    S_ynF(Ny) = max_nsp;
    S_ywF(Ny) = max_w;
    S_ygF(Ny) = max_g;
    if(max_g<0.1)
        break;
    end
    if(Ny==center_N)
        tmp_nsp = max_nsp;
        tmp_w = max_w;
        tmp_g = max_g;
    end
end

nsp_center = tmp_nsp;
w_center = tmp_w;
cemter_g = tmp_g;
for Ny = center_N-1:-1:1
    Ny
    yita0 = nY(Ny);
    beta0 = yita0*1e6*(e/me/(ve*c)/W0);
    iv_nsp = nsp_center;
    iv_w = w_center;
    iv_g = g_center;
    
    scan_max
    
    nsp_center = max_nsp;
    w_center = max_w;
    g_center = max_g;
    
    S_ynF(Ny) = max_nsp;
    S_ywF(Ny) = max_w;
    S_ygF(Ny) = max_g;
    if(max_g<0.1)
        break;
    end
end

fclose(f_log);

f_max = fopen('./res/max_yF.txt','w');
for Ny=1:numy
    fprintf(f_max,'%d\t%.4f\t%.4f\t%.4f\r\n',nY(Ny),S_ynF(Ny),S_ywF(Ny),S_ygF(Ny));
end
fclose(f_max);

open('./res/max_yF.txt');