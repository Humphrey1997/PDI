st_nsp = iv_nsp;

nsp = st_nsp;
para_update_s
single_solve
st_w = res_w;
st_g = res_g;

nsp = st_nsp - gap_nsp;
para_update_s
single_solve
left_w = res_w;
left_g = res_g;

nsp_right = nsp*2;
nsp_left = nsp/2;

if(st_g<=left_g)
    dir = -1;
    nsp_right = st_nsp;
    max_w = left_w;
    max_g = left_g;
    max_nsp = st_nsp - gap_nsp;
else
    dir = 1;
    nsp_left = st_nsp - gap_nsp;
    max_w = st_w;
    max_g = st_g;
    max_nsp = st_nsp;
end

nsp = max_nsp;
pre_nsp = nsp;
pre_w = max_w;
pre_g = max_g;
gap = gap_nsp*5;

record_rad_w =  default_rad_w;
record_rad_g =  default_rad_g;

while(1)
    fprintf(f_log,'%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\r\n',yita0,nsp_left,nsp_right,max_nsp,max_w,max_g);
    
    nsp = nsp + dir * gap;
    iv_w = pre_w;
    iv_g = pre_g;
    para_update_s
    single_solve
    
    if(res_g<pre_g)
        if(dir>0)
            nsp_right = nsp;
        else
            nsp_left = nsp;
        end
        dir = -dir;
        gap = round(gap*0.6/gap_nsp)*gap_nsp;
        if(gap<=gap_nsp)
            gap = gap_nsp;
        end
    elseif(nsp+1.5*dir*gap_nsp>nsp_right||nsp+1.5*dir*gap_nsp<nsp_left)
        dir = -dir;
        gap = round(gap*0.6/gap_nsp)*gap_nsp;
        if(gap<=gap_nsp)
            gap = gap_nsp;
        end
    end
    
    if(res_g>max_g)
        if(nsp>max_nsp)
            nsp_left = max_nsp;
        end
        if(nsp<max_nsp)
            nsp_right = max_nsp;
        end
        max_w = res_w;
        max_g = res_g;
        max_nsp = nsp;
    else
        if(nsp<nsp_right&&nsp>max_nsp)
            nsp_right = nsp;
        end
        if(nsp>nsp_left&&nsp<max_nsp)
            nsp_left = nsp;
        end
    end
    
    pre_w = res_w;
    pre_g = res_g;
    pre_nsp = nsp;
    
    if((max_nsp-nsp_left)<1.5*gap_nsp&&(nsp_right-max_nsp)<1.5*gap_nsp)
        break;
    end
    
    default_rad_w = 2*abs(pre_w-st_w);
    default_rag_g = 2*abs(max_g-st_g);
end

nsp_max = max_nsp;
default_rad_w = record_rad_w;
default_rad_g = record_rad_g;