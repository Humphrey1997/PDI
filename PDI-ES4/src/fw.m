function res = fw(ww)
w = ww(1)*1e-4+ww(2)*1e-6*1i;
res = -cond(PIM3(w,2));
end

