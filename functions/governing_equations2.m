% This functions execute the governing equations in the model

function [out, lam_ds] = governing_equations2(t, in, M, i, s, r, p, sel, agg, growth,hivpoints)
%#codegen
invec = in(1:i.nstates);

% -- Get force of infection
lam_ds =  M.lambda_ds*invec;
lam_mdr =  M.lambda_dr*invec;

% -- If year is > 1980 get the HIV acquiaiton rate into the model
foi=0;
if t>=1980
    % foi = get_foi_fun(t,r,hivpoints);

    hiv = 0;
    yrstart = 1980;
    hivIR = 0;
    if (t >= yrstart)
        hiv = 1;
        it = round(t);% integer of t to used as index in pre-def vectors.
        ii = (it - yrstart)+1;
        iii = ((it - yrstart) + 1)+1;

        x0 = it;
        x1 = it + 1;


        % Faster performance
        hivIR = (hivpoints(ii)*(x1-t) + hivpoints(iii)*(t-x0))/(x1-x0)*r.turnoffHIV;

    end

    foi=hivIR*hiv;

end

% -- Get all model components together

allmat = M.lin + M.Dxlin + ...
    M.nlinHIV.*foi + ...
    (M.nlin.a0_4.ds.*lam_ds(1))+...
    (M.nlin.a5_9.ds.*lam_ds(2))+...
    (M.nlin.a10_15.ds.*lam_ds(3))+...
    (M.nlin.a15p.ds.*lam_ds(4))+...
    (M.nlin.a0_4.mdr.*lam_mdr(1))+...
    (M.nlin.a5_9.mdr.*lam_mdr(2))+...
    (M.nlin.a10_15.mdr.*lam_mdr(3))+...
    (M.nlin.a15p.mdr.*lam_mdr(4));


% -- Evaluate numbers transitioning in this time step
out = allmat*invec;

% -- Implement deaths
morts = sum(M.mortvec,2).*invec;
morts_slum = sum(M.mortvec(s.slum,:),2).*invec(s.slum);
out = out - morts;

% -- Implement births
births = sum(morts)*p.popn_turnover;
births_slum = sum(morts_slum)*p.popn_turnover;

if(p.bra+p.geo > 0)
    out(i.U.a15p.slum.hneg) = out(i.U.a15p.slum.hneg)+births_slum*p.birth_slum;
else
    out(i.U.a0_4.slum.hneg) = out(i.U.a0_4.slum.hneg)+births_slum*p.birth_slum;
end

out(i.U.a0_4.no_slum.hneg) = out(i.U.a0_4.no_slum.hneg) + ((births-births_slum)*p.birth_slum + births*(1-p.birth_slum)) ;


% -- Get the auxiliaries
auxil=zeros(i.nx,1);

auxil(i.aux.inc)      = 0;% (agg.inc*(sel.inc.*allmat))*invec;
auxil(i.aux.inc(2))   = 0;%(agg.inc(2,:)*(sel.inc.*allmat))*invec + sum((sel.acqu.*allmat)*invec);
auxil(i.aux.notif)    = 0;%(agg.notif*(sel.notif.*allmat))*invec;
auxil(i.aux.mort(1))  = 0;%sum(M.mortvec(:,2).*invec);
auxil(i.aux.mort(2))  = 0;%sum(M.mortvec(:,3).*invec);
auxil(i.aux.pt)       = 0;%(agg.ipt *(sel.ipt.*allmat))*invec;
auxil(i.aux.newart)   = 0;%(agg.art *(sel.art.*allmat))*invec;

% HIV inc
% if (t>1979)
% auxil(i.aux.hiv)  = 0;%(agg.hiv*(sel.hiv.*allmat))*invec;
% else
auxil(i.aux.hiv)  = 0;
% end


%%auxilput interventions set to 0
auxil(i.aux.inc_ds)= 0;
auxil(i.aux.inc_dr)= 0;
auxil(i.aux.screen_all)= 0;
auxil(i.aux.screen_plhiv)= 0;
auxil(i.aux.screen_hhc)= 0;
auxil(i.aux.screen_slum)= 0;
auxil(i.aux.tst_all)= 0;
auxil(i.aux.tst_plhiv)= 0;
auxil(i.aux.tst_hhc)= 0;
auxil(i.aux.tst_slum)= 0;
auxil(i.aux.tpt_all)= 0;
auxil(i.aux.tpt_plhiv)= 0;
auxil(i.aux.tpt_hhc)= 0;
auxil(i.aux.tpt_slum)= 0;
auxil(i.aux.tx_fl_all)=  0;
auxil(i.aux.tx_fl_plhiv)=  0;
auxil(i.aux.tx_fl_hhc)=  0;
auxil(i.aux.tx_fl_slum)=  0;
auxil(i.aux.tx_sl_all)=  0;
auxil(i.aux.tx_sl_plhiv)=  0;
auxil(i.aux.tx_sl_hhc)=  0;
auxil(i.aux.tx_sl_slum)=  0;
auxil(i.aux.daly_all)=   0;
auxil(i.aux.daly_plhiv)=   0;
auxil(i.aux.daly_slum)=   0;
auxil(i.aux.yll_all)=   0;
auxil(i.aux.yll_plhiv)=   0;
auxil(i.aux.yll_slum)=   0;
auxil(i.aux.yld_all)=   0 ;
auxil(i.aux.yld_plhiv)=   0 ;
auxil(i.aux.yld_slum)=   0 ;
auxil(i.aux.mu_ds)= 0;
auxil(i.aux.mu_dr)=  0;

out=[out;auxil((i.nstates+1):i.nx) ];

