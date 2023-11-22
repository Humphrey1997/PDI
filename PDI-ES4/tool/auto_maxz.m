fclose all;
warning off;
T1 = clock;
EBW
f_log = fopen('./res/log.txt','w');

Del = 180;
center_nsp = 0.1344;
center_nz = 0.0042;
center_w = 4.8;
center_g = 3.0;

nz_start = 0;
nz_end = 0.01;
nz_gap = 2e-4;
numz = round((nz_end-nz_start)/(nz_gap))+1;
nZ = linspace(nz_start,nz_end,numz)+1e-6;

center_N = round((center_nz-nz_start)/nz_gap)+1;

S_znF = zeros(1,numz);
S_zwF = zeros(1,numz);
S_zgF = zeros(1,numz);

default_rad_w = 0.01;
default_rad_g = 0.02;
gap_nsp = 1e-4;

for Nz = center_N:numz
    Nz
    nsz = nZ(Nz);
    iv_nsp = center_nsp;
    iv_w = center_w;
    iv_g = center_g;
    
    scan_max
    
    center_nsp = max_nsp;
    center_w = max_w;
    center_g = max_g;
    
    S_znF(Nz) = max_nsp;
    S_zwF(Nz) = max_w;
    S_zgF(Nz) = max_g;
    if(max_g<0.02)
        break;
    end
    if(Nz==center_N)
        tmp_nsp = max_nsp;
        tmp_w = max_w;
        tmp_g = max_g;
    end
end

center_nsp = tmp_nsp;
center_w = tmp_w;
cemter_g = tmp_g;
for Nz = center_N-1:-1:1
    Nz
    nsz = nZ(Nz);
    iv_nsp = center_nsp;
    iv_w = center_w;
    iv_g = center_g;
    
    scan_max
    
    center_nsp = max_nsp;
    center_w = max_w;
    center_g = max_g;
    
    S_znF(Nz) = max_nsp;
    S_zwF(Nz) = max_w;
    S_zgF(Nz) = max_g;
    if(max_g<0.02)
        break;
    end
end

fclose(f_log);

f_max = fopen('./res/max_zF.txt','w');
for Nz=1:numz
    fprintf(f_max,'%d\t%.4f\t%.4f\t%.4f\r\n',nZ(Nz),S_znF(Nz),S_zwF(Nz),S_zgF(Nz));
end
fclose(f_max);

open('./res/max_zF.txt');