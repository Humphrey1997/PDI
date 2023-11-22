fclose all;
EBW
data = importdata('./res/res_wk.txt');
LW = 0.1;
RW = 4;

Del = 0;
nz_start = 0;
nz_end = 0.015;
nz_gap = 2e-4;
N = round((nz_end-nz_start)/(nz_gap))+1;
nZ = linspace(nz_start,nz_end,N)+1e-6;
np_start = 19.82;
np_end = 19.898;
np_gap = 2e-3;
M = round((np_end-np_start)/(np_gap))+1;
nP = linspace(np_start,np_end,M);
center_rad_w = 0.05;
center_rad_g = 1;

data_w = zeros(M,N);
data_g = zeros(M,N);

for si=1:M
    for sj=1:N
        data_w(si,sj) = data(si,2*sj-1);
        data_g(si,sj) = data(si,2*sj);
    end
end

f_log = fopen('./res/log.txt','w');
NotFinish = true;
while(NotFinish)
    NotFinish = false;
    for sj=1:N
        for si=1:M
            if(data_g(si,sj)>0.1&&data_w(si,sj)<RW&&data_w(si,sj)>LW)
                if((si>1&&data_g(si-1,sj)<1e-3)||(si<M&&data_g(si+1,sj)<1e-3))
                    NotFinish = true;
                    mid_i = si;
                    mid_j = sj;
                    mid_w = data_w(si,sj);
                    mid_g = data_g(si,sj);
                    scan_peak
                elseif(si==1||si==M||(si>1&&si<M&&(data_g(si,sj)>=data_g(si-1,sj)||data_g(si,sj)>=data_g(si+1,sj))))
                    if(sj>1&&data_g(si,sj-1)<1e-3)
                        NotFinish = true;
                        mid_i = si;
                        mid_j = sj-1;
                        mid_w = data_w(si,sj);
                        mid_g = data_g(si,sj);
                        scan_peak
                    end
                    if(sj<N&&data_g(si,sj+1)<1e-3)
                        NotFinish = true;
                        mid_i = si;
                        mid_j = sj+1;
                        mid_w = data_w(si,sj);
                        mid_g = data_g(si,sj);
                        scan_peak
                    end
                end
            end
        end
    end
end

fclose(f_log);