mid_i = center_i;
mid_w = center_w;
mid_g = center_g;
if(sign_j==1)
    for sj=center_j+1:1:max_j
        mid_j = sj;
        scan_peak
        [max_g,max_index] = max(data_g(:,mid_j));
        if(max_g<0.02)
            break
        end
        gap_i = max_index - mid_i;
        mid_i = max_index + gap_i;
        mid_w = data_w(mid_i,mid_j);
        mid_g = data_g(mid_i,mid_j);
    end
else
    for sj=center_j-1:-1:min_j
        mid_j = sj;
        scan_peak
        [max_g,max_index] = max(data_g(:,mid_j));
        if(max_g<0.02)
            break
        end
        gap_i = max_index - mid_i;
        mid_i = max_index + gap_i;
        mid_w = data_w(mid_i,mid_j);
        mid_g = data_g(mid_i,mid_j);
    end
end