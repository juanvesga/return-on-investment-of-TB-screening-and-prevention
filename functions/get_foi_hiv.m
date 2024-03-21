function FOI = get_foi_hiv(t,r,points)

hiv = 0;
yrstart = 1980;
hivIR = 0;
if (t >= yrstart)
    hiv = 1;
    it = round(t);% integer of t to used as index in pre-def vectors.
    ii = (it - yrstart)+1;
    % if (t > 2021)
    %     ii = (2020 - yrstart)+1;
    % end
    iii = ((it - yrstart) + 1)+1;
    % if (t > 2021)
    %     iii = ((2020 - yrstart) + 1)+1;
    % end
    x0 = it;
    x1 = it + 1;
    
    % Built in
    %hivIR = interp1([x0 x1], [points(ii), points(iii)],t)*r.turnoffHIV;
    
    % Faster performance
    hivIR = (points(ii)*(x1-t) + points(iii)*(t-x0))/(x1-x0)*r.turnoffHIV;
    
    % -- When year is > 2017 use decline as stated

%     if(t>2018)
%         f=9;
% 
%     end
% 
% if(t>2020)
%         f=9;
% 
%     end


    % 
    % 
    % if (t > 2022)
    %     incdecline = (1 + (t - 2022)) * r.hivdecline;
    %     hivIR = hivIR * (1 - incdecline);
    % end
    
end

FOI=hivIR*hiv;