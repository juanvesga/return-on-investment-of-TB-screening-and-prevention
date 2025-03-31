function out = goveqs_scaleup(t, in, M0, M1, times, i, s,r, p, sel, agg,hivpoints)
    scale = min((t-times(1))/(times(2)-times(1)),1); 
    if (scale<0),scale=0;end
    Mt = M1;
    Mt.lin = M0.lin + scale*(M1.lin-M0.lin);
    % Mt.Dxlin = M0.Dxlin + scale*(M1.Dxlin-M0.Dxlin);

    
out = governing_equations2(t, in, Mt, i, s,r,p, sel, agg, p.growth,hivpoints);

        