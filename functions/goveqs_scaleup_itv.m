function out = goveqs_scaleup_itv(t, in, M0, M1, times, i, s,r, p, sel, agg,hivpoints)
    scale = min((t-times(1))/(times(2)-times(1)),1); 

    if (scale<0),scale=0;end
    Mt = M1;
    Mt.lin = M0.lin + scale*(M1.lin-M0.lin); 
    Mt.linHIV = M0.linHIV + scale*(M1.linHIV-M0.linHIV); 
    Mt.cfy_U=M0.cfy_U + scale*(M1.cfy_U-M0.cfy_U);
    Mt.cfy_Lf=M0.cfy_Lf + scale*(M1.cfy_Lf-M0.cfy_Lf);
    Mt.cfy_Ls=M0.cfy_Ls + scale*(M1.cfy_Ls-M0.cfy_Ls);
    Mt.cfy_I=M0.cfy_I + scale*(M1.cfy_I-M0.cfy_I);
    Mt.screen_U=M0.screen_U + scale*(M1.screen_U-M0.screen_U);
    Mt.screen_Lf=M0.screen_Lf + scale*(M1.screen_Lf-M0.screen_Lf);
    Mt.screen_Ls=M0.screen_Ls + scale*(M1.screen_Ls-M0.screen_Ls);
    Mt.screen_I=M0.screen_I + scale*(M1.screen_I-M0.screen_I);
    Mt.screen_hiv=M0.screen_hiv + scale*(M1.screen_hiv-M0.screen_hiv);

out = governing_equations_itv2(t, in, Mt, i, s,r,p, sel, agg, p.growth,hivpoints);
     



        