function res = fw_par(ww)
w = ww(1)*1e-3+ww(2)*1e-4*1i;
res = -cond(PIM(w,2));
end
