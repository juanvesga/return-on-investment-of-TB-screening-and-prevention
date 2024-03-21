
function [r,p] = allocate_parameters(x,r,p,xi)

r.beta          =  x(xi.r_beta)/p.scale(xi.r_beta);
rf_beta_mdr      = x(xi.rf_beta_mdr)/p.scale(xi.rf_beta_mdr);
rf_mort_TB      =   x(xi.rf_mort_TB)/p.scale(xi.rf_mort_TB);
rf_mort_TBh     =x(xi.rf_mort_TBh)/p.scale(xi.rf_mort_TBh);
r.RRslum        = x(xi.r_RRslum)/p.scale(xi.r_RRslum);
r.RRhiv         = x(xi.r_RRhiv )/p.scale(xi.r_RRhiv );
rf_careseeking_mo= 12;% x(xi.rf_careseeking_mo )/p.scale(xi.rf_careseeking_mo);
p.Tx_init       = x(xi.p_Tx_init)/p.scale(xi.p_Tx_init);
p.Tx_init2       = x(xi.p_Tx_init)/p.scale(xi.p_Tx_init);
r.ART_init      =  x(xi.r_ART_init)/p.scale(xi.r_ART_init);
p.IPThiv        =  x(xi.p_IPThiv)/p.scale(xi.p_IPThiv);
p.ntpcov       =  x(xi.p_ntpcov_dr)/p.scale(xi.p_ntpcov_dr);
p.ntpcov_dr    =  x(xi.p_ntpcov_dr)/p.scale(xi.p_ntpcov_dr);
p.hhcovu5      =  x(xi.p_hhcovu5)/p.scale(xi.p_hhcovu5);
p.hhcov        =  x(xi.p_hhcov)/p.scale(xi.p_hhcov);


r.beta_mdr = r.beta *rf_beta_mdr;

%Fixies
r.f_self_cure   = 1;%x(xi.rf_self_cure)/p.scale(xi.rf_self_cure);

r.f_reactivation= 1;%x(xi.rf_reactivation)/p.scale(xi.rf_reactivation);
p.HIVlam        = 1;% x(xi.p_HIVlam)/p.scale(xi.p_HIVlam);
r.HIV_mort      = 1;%x(xi.r_HIV_mort)/p.scale(xi.r_HIV_mort);



%% Construct morta;ity
r.mort_TB(1)    = r.mort_TB0*rf_mort_TB;
r.mort_TB(2)    = r.mort_TB0*rf_mort_TBh;
r.mort_TB(3)    = r.mort_TB0*rf_mort_TBh*r.ARTred;

%Careseking
r.careseeking=12/rf_careseeking_mo;

% % % %% pass manual fitting values
% x(xi.r_beta) = r.beta* p.scale(xi.r_beta)
% x(xi.rf_beta_mdr) = rf_beta_mdr* p.scale(xi.rf_beta_mdr)
% x(xi.rf_mort_TB) = rf_mort_TB* p.scale(xi.rf_mort_TB)
% x(xi.rf_mort_TBh) = rf_mort_TBh* p.scale(xi.rf_mort_TBh)
% x(xi.r_RRslum) = r.RRslum * p.scale(xi.r_RRslum)
% x(xi.r_RRhiv) = r.RRhiv * p.scale(xi.r_RRhiv)
% x(xi.rf_careseeking_mo) = rf_careseeking_mo * p.scale(xi.rf_careseeking_mo)
% x(xi.p_Tx_init) = p.Tx_init * p.scale(xi.p_Tx_init)
% x(xi.p_Tx_init2) = p.Tx_init * p.scale(xi.p_Tx_init2)
% x(xi.r_ART_init) = r.ART_init* p.scale(xi.r_ART_init)
% x(xi.p_IPThiv) = p.IPThiv * p.scale(xi.p_IPThiv)
% x(xi.p_ntpcov) = p.ntpcov * p.scale(xi.p_ntpcov)
% x(xi.p_ntpcov_dr) = p.ntpcov_dr * p.scale(xi.p_ntpcov_dr)
% x(xi.p_hhcovu5) = p.hhcovu5 * p.scale(xi.p_hhcovu5)
% x(xi.p_hhcov) = p.hhcov * p.scale(xi.p_hhcov)
% % % % % %
% save_results(x,'bestset','ZAF','mle');

% Construct vector for progression from 'fast' latent infection, each
% element representing a different HIV status: 1.HIV-ve, 2.HIV+ve, 3.On ART
tmp = r.progression0(1)*[1 1 1;1 1 1];
tmp(:,[2,3]) = tmp(:,[2,3])*r.RRhiv .*[1 r.ARTred ;1 r.ARTred ];
tmp(2,:)=tmp(2,:).*r.RRslum;
r.progression = tmp;

% Construct vector for progression from 'slow' latent infection, each
% element representing a different HIV status: 1.HIV-ve, 2.HIV+ve, 3.On ART
tmp = r.reactivation0(1)*[1 1 1;1 1 1];
tmp(:,[2,3]) = tmp(:,[2,3])*r.RRhiv .*[1 r.ARTred ;1 r.ARTred] ;
tmp(2,:)=tmp(2,:).*r.RRslum;
r.reactivation = tmp;

r.self_cure([1,3])    = r.self_cure0([1,3]);
% p.imm                 = p_imm * [1 0 1];

