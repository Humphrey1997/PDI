%单点求解器
Tol = 1e-3;
solve_flag = false;
min_g = 0.01;
rad_w = default_rad_w;
rad_g = default_rad_g;
while(~solve_flag)
    if(iv_g<1.5*min_g)
        iv_g = 1.5*min_g;
        lb = [iv_w-rad_w,min_g];
        ub = [iv_w+rad_w,2*min_g];
    elseif(iv_g<rad_g)
        lb = [iv_w-rad_w,min_g];
        ub = [iv_w+rad_w,iv_g+rad_g];
    else
        lb = [iv_w-rad_w,iv_g-rad_g];
        ub = [iv_w+rad_w,iv_g+rad_g];
    end
    
    wiv = [iv_w,iv_g];
    options = optimoptions('patternsearch','TolMesh',Tol,'Display','off');
    [wx,fval,exitflag] = patternsearch(@fw,wiv,[],[],[],[],lb,ub,[],options);
    
    solve_flag = true;
    if(iv_g<2*min_g)&&(abs(wx(2)-2*min_g)<0.001)
        iv_g = 2*min_g;
        solve_flag = false;
    end
    if((abs(wx(1)-iv_w)>rad_w*0.95)||(abs(wx(2)-iv_g)>rad_g*0.95))
        if(abs(wx(1)-iv_w)>rad_w*0.95)
            iv_w = wx(1)+(wx(1)-iv_w)/abs(wx(1)-iv_w)*(rad_w*0.75);
        end
        if(abs(wx(2)-iv_g)>rad_g*0.95)
            rad_g = rad_g+0.1;
            iv_g = wx(2)+(wx(2)-iv_g)/abs(wx(2)-iv_g)*(rad_g*0.75);
        end
        solve_flag = false;
    end
    if(abs(wx(2)-iv_g)<Tol)
        Tol = Tol/10;
        iv_w = wx(1);
        rad_g = rad_g/10;
        rad_w = rad_w/10;
        solve_flag = false;
    end
    
    if(solve_flag)
        [Del nsz nsp n1p wx(1) wx(2)]
        fprintf(f_log,'%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\r\n',Del,nsz,nsp,n1p,wx(1),wx(2));
        res_w = wx(1);
        res_g = wx(2);
    end
end