fclose all;
warning off;
T1 = clock;
EBW
f_log = fopen('./res/log.txt','w');
f_max = fopen('./res/max.txt','w');

maxN = 5;
S_nsp = zeros(1,maxN);
S_w = zeros(1,maxN);
S_g = zeros(1,maxN);

Ni_start = 1;
Ni_min = 1;
Ni_max = maxN;


Del = 0;
start_nsp = 19.898;
start_w = 2.42;
start_g = 5.2;
default_rad_w = 0.1;
default_rad_g = 1;
gap_nsp = 1e-4;

Ni = Ni_start;
while(Ni<=Ni_end)
    iv_w = start_w*Ni;
    iv_g = start_g*sqrt(Ni);
    iv_nsp = start_nsp+0.076*(Ni-1);
    scan_max
    nsp = nsp_max;
    fprintf(f_max,'%d\t%.3f\t%.3f\t%.3f\r\n',Ni,nsp,max_w,max_g);
    Ni = Ni + 1;
end
fclose(f_log);
fclose(f_max);