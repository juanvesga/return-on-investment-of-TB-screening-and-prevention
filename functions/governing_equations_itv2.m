% This functions execute the governing equations in the model

function [out, lam_ds] = governing_equations_itv2(t, in, M, i, s, r, p, sel, agg, growth,hivpoints)

invec = in(1:i.nstates);
% N= sum(invec(s.nstates));

% -- Get force of transition for CYF
fot_u = M.fot*invec/0.485;  % U
fot_lf = M.fot*invec/0.3605;% Lf
fot_ls = M.fot*invec/0.1545;% Ls
fot_i = M.fot*invec/0.0310; % I



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

%SLum switch
slum_itv=0;
if t>=2024 && t <= 2027
    slum_itv = 1;
end

% -- Get all model components together


allmat = M.lin + M.slumlin.*slum_itv + M.linHIV+...
    M.nlinHIV.*foi + ...
    (M.nlin.a0_4.ds.*lam_ds(1))+...
    (M.nlin.a5_9.ds.*lam_ds(2))+...
    (M.nlin.a10_15.ds.*lam_ds(3))+...
    (M.nlin.a15p.ds.*lam_ds(4))+...
    (M.nlin.a0_4.mdr.*lam_mdr(1))+...
    (M.nlin.a5_9.mdr.*lam_mdr(2))+...
    (M.nlin.a10_15.mdr.*lam_mdr(3))+...
    (M.nlin.a15p.mdr.*lam_mdr(4))+...
    ((M.cfy_U).*fot_u) +...
    ((M.cfy_Lf).*fot_lf) +...
    ((M.cfy_Ls).*fot_ls) +...
    ((M.cfy_I).*fot_i);


% -- Evaluate numbers transitioning in this time step
out = allmat*invec;

% -- Implement deaths
morts =  M.mu.*invec;
out = out - morts;
morts_slum =  M.mu(s.slum).*invec(s.slum);
births = sum(morts)*p.popn_turnover;
births_slum = sum(morts_slum)*p.popn_turnover;
if(p.bra+p.geo > 0)
    out(i.U.a15p.slum.hneg) = out(i.U.a15p.slum.hneg)+births_slum*p.birth_slum;
else
    out(i.U.a0_4.slum.hneg) = out(i.U.a0_4.slum.hneg)+births_slum*p.birth_slum;
end
out(i.U.a0_4.no_slum.hneg) = out(i.U.a0_4.no_slum.hneg) + (births-births_slum)*p.birth_slum + (births*(1-p.birth_slum)) ;

% -- Get the auxiliaries
auxil=zeros(i.nx,1);
if (t<2022)
    auxil(i.aux.inc)      = (agg.inc*(sel.inc.*allmat))*invec;
    auxil(i.aux.inc(2))   = (agg.inc(2,:)*(sel.inc.*allmat))*invec + sum((sel.acqu.*allmat)*invec);
    auxil(i.aux.notif)    = (agg.notif*(sel.notif.*allmat))*invec;
    auxil(i.aux.mort(1))  = sum(M.mortvec(:,2).*invec);
    auxil(i.aux.mort(2))  = sum(M.mortvec(:,3).*invec);
    auxil(i.aux.pt)       = (agg.ipt *(sel.ipt.*allmat))*invec;
    auxil(i.aux.newart)   = (agg.art *(sel.art.*allmat))*invec;
else
    auxil(i.aux.inc)      = 0;
    auxil(i.aux.inc(2))   = 0;
    auxil(i.aux.notif)    = 0;
    auxil(i.aux.mort(1)) = 0;
    auxil(i.aux.mort(2)) = 0;
    auxil(i.aux.pt)       = 0;
    auxil(i.aux.newart)   = 0;


end

% auxil(i.aux.notif)    = (agg.notif*(sel.notif.*allmat))*invec;

auxil(i.aux.tpt_hhc)=    (agg.tpt_hhc*(sel.ipt.*(M.cfy_U.*fot_u + M.cfy_Lf.*fot_lf + M.cfy_Ls.*fot_ls)))*invec;


if (t<2021)

    % if (t>1979)
    %     auxil(i.aux.hiv)  = (agg.hiv*(sel.hiv.*allmat))*invec;
    % else
         auxil(i.aux.hiv)  = 0;
    % end
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

else
    auxil(i.aux.hiv)  = 0;
    screenMat= M.screen_hiv + M.screen_slum.*slum_itv + fot_u*M.screen_U + fot_lf*M.screen_Lf+...
        fot_ls*M.screen_Ls + fot_i*M.screen_I;

    %%auxilput interventions set to 0
    auxil(i.aux.inc_ds)= (agg.inc_ds*(sel.inc_ds.*allmat))*invec;
    auxil(i.aux.inc_dr)= (agg.inc_dr*(sel.inc_dr.*allmat))*invec;
    auxil(i.aux.screen_all)= (agg.screen_all*(sel.screen.*screenMat))*invec;
    auxil(i.aux.screen_plhiv)= (agg.screen_plhiv*(sel.screen.*M.screen_hiv))*invec;
    auxil(i.aux.screen_hhc)= (agg.screen_hhc*(sel.screen.*(M.screen_U.*fot_u +...
        M.screen_Lf.*fot_lf*+M.screen_Ls.*fot_ls + M.screen_I.*fot_i)))*invec;
    auxil(i.aux.screen_slum)= (agg.screen_slum*(sel.screen.*(M.screen_slum.*slum_itv)))*invec;

    auxil(i.aux.tst_all)= (agg.tst_all*(sel.screen.*screenMat))*invec;
    auxil(i.aux.tst_plhiv)= (agg.tst_plhiv*(sel.screen.*M.screen_hiv))*invec;
    auxil(i.aux.tst_hhc)= (agg.tst_hhc*(sel.screen.*(M.cfy_U.*fot_u +...
        M.cfy_Lf.*fot_lf + M.cfy_Ls.*fot_ls + M.cfy_I.*fot_i)))*invec;
    auxil(i.aux.tst_slum)= (agg.tst_slum*(sel.screen.*(M.screen_slum.*slum_itv)))*invec;

    auxil(i.aux.tpt_all)=    (agg.tpt_all*(sel.ipt.*allmat))*invec;
    auxil(i.aux.tpt_plhiv)=    (agg.tpt_all*(sel.ipt.*M.linHIV))*invec;
    auxil(i.aux.tpt_slum)=    (agg.tpt_slum*(sel.ipt.*M.slumlin.*slum_itv))*invec;


    auxil(i.aux.tx_fl_all)=  (agg.tx_fl_all*(sel.tx_fl.*allmat))*invec;
    auxil(i.aux.tx_fl_plhiv)=  (agg.tx_fl_plhiv*(sel.tx_fl.*allmat))*invec;
    auxil(i.aux.tx_fl_hhc)=  (agg.tx_fl_hhc*(sel.tx_fl.* M.cfy_I.*fot_i))*invec;
    auxil(i.aux.tx_fl_slum)=    (agg.tx_fl_slum*(sel.tx_fl.*allmat))*invec;


    auxil(i.aux.tx_sl_all)=  (agg.tx_sl_all*(sel.tx_sl.*allmat))*invec;
    auxil(i.aux.tx_sl_plhiv)=  (agg.tx_sl_plhiv*(sel.tx_sl.*allmat))*invec;
    auxil(i.aux.tx_sl_hhc)=  (agg.tx_sl_hhc*(sel.tx_sl.*M.cfy_I.*fot_i))*invec;
    auxil(i.aux.tx_sl_slum)=    (agg.tx_sl_slum*(sel.tx_sl.*allmat))*invec;

    auxil(i.aux.daly_all)=  0;% (agg.yld_all*invec) + (agg.yll_all*morts);
    auxil(i.aux.daly_plhiv)= 0;%    (agg.yld_plhiv*invec) + (agg.yll_plhiv*morts);
    auxil(i.aux.daly_slum)= 0;%  (agg.yld_slum*invec) + (agg.yll_slum*morts);
    auxil(i.aux.yll_all)=   agg.yll_all*morts;
    auxil(i.aux.yll_plhiv)=   agg.yll_plhiv*morts;
    auxil(i.aux.yll_slum)=   agg.yll_slum*morts;
    auxil(i.aux.yld_all)=   agg.yld_all*invec ;
    auxil(i.aux.yld_plhiv)=   agg.yld_plhiv*invec ;
    auxil(i.aux.yld_slum)=   agg.yld_slum*invec ;
    auxil(i.aux.mu_ds)=  agg.mu_ds*morts;
    auxil(i.aux.mu_dr)=  agg.mu_dr*morts;
end
out=[out;auxil((i.nstates+1):i.nx) ];