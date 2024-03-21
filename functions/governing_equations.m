% This functions execute the governing equations in the model 

function [out, lam_ds] = governing_equations(t, in, M, i, s, r, p, sel, agg, growth,hivpoints)
%#codegen
invec = in(1:i.nstates);

% -- Get force of infection
lam_ds =  M.lambda_ds*invec;
lam_mdr =  M.lambda_dr*invec;

% -- If year is > 1980 get the HIV acquiaiton rate into the model 
foi=0;
if t>=1980
    foi = get_foi_hiv(t,r,hivpoints);
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

out(i.U.a0_4.no_slum.hneg) = out(i.U.a0_4.no_slum.hneg);% + ((births-births_slum)*p.birth_slum + births*(1-p.birth_slum)) ;


% -- Get the auxiliaries
out(i.aux.inc)      =0;% (agg.inc*(sel.inc.*allmat))*invec;
out(i.aux.inc(2))   =0;% (agg.inc(2,:)*(sel.inc.*allmat))*invec + sum((sel.acqu.*allmat)*invec);
out(i.aux.notif)    =0;% (agg.notif*(sel.notif.*allmat)*invec);
out(i.aux.mort(1))  =0;% sum(M.mortvec(:,2).*invec);
out(i.aux.mort(2))  =0;% sum(M.mortvec(:,3).*invec);
out(i.aux.pt)       =0;% (agg.ipt *(sel.ipt.*allmat))*invec;
out(i.aux.newart)   =0;% (agg.art *(sel.art.*allmat))*invec;

% HIV inc
% if (t>1979)
% out(i.aux.hiv)  = 0;%(agg.hiv*(sel.hiv.*allmat))*invec;
% else
out(i.aux.hiv)  = 0;
% end  


%%Output interventions set to 0
    out(i.aux.inc_ds)= 0;
    out(i.aux.inc_dr)= 0;
    out(i.aux.screen_all)= 0;
    out(i.aux.screen_plhiv)= 0;
    out(i.aux.screen_hhc)= 0;
    out(i.aux.screen_slum)= 0;
    out(i.aux.tst_all)= 0;
    out(i.aux.tst_plhiv)= 0;
    out(i.aux.tst_hhc)= 0;
    out(i.aux.tst_slum)= 0;
    out(i.aux.tpt_all)= 0;
    out(i.aux.tpt_plhiv)= 0;
    out(i.aux.tpt_hhc)= 0;
    out(i.aux.tpt_slum)= 0;
    out(i.aux.tx_fl_all)=  0;
    out(i.aux.tx_fl_plhiv)=  0;
    out(i.aux.tx_fl_hhc)=  0;
    out(i.aux.tx_fl_slum)=  0;
    out(i.aux.tx_sl_all)=  0;
    out(i.aux.tx_sl_plhiv)=  0;
    out(i.aux.tx_sl_hhc)=  0;
    out(i.aux.tx_sl_slum)=  0;
    out(i.aux.daly_all)=   0;
    out(i.aux.daly_plhiv)=   0;
    out(i.aux.daly_slum)=   0;
    out(i.aux.yll_all)=   0;
    out(i.aux.yll_plhiv)=   0;
    out(i.aux.yll_slum)=   0;
    out(i.aux.yld_all)=   0 ;
    out(i.aux.yld_plhiv)=   0 ;
    out(i.aux.yld_slum)=   0 ;
    out(i.aux.mu_ds)= 0;
    out(i.aux.mu_dr)=  0;

