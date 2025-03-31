
function [r,p] = allocate_parameters(x,r,p,xi)

r.beta            =  exp(x(xi.r_beta));%(xi.r_beta);
rf_beta_mdr       =  exp(x(xi.rf_beta_mdr));%(xi.rf_beta_mdr);
rf_mort_TB        =  exp(x(xi.rf_mort_TB));%(xi.rf_mort_TB);
rf_mort_TBh       =  exp(x(xi.rf_mort_TBh));%(xi.rf_mort_TBh);
r.RRslum          =  exp(x(xi.r_RRslum));%(xi.r_RRslum);
r.RRhiv           =  exp(x(xi.r_RRhiv ));%(xi.r_RRhiv );
rf_careseeking_mo =  exp(x(xi.rf_careseeking_mo));%(xi.rf_careseeking_mo);
p.Tx_init         =  exp(x(xi.p_Tx_init));%(xi.p_Tx_init);
p.Tx_init2        =  exp(x(xi.p_Tx_init2));%(xi.p_Tx_init);
r.ART_init        =  exp(x(xi.r_ART_init));%(xi.r_ART_init);
p.IPThiv          =  exp(x(xi.p_IPThiv));%(xi.p_IPThiv);
p.ntpcov          =  exp(x(xi.p_ntpcov));%(xi.p_ntpcov_dr);
p.ntpcov_dr       =  exp(x(xi.p_ntpcov_dr));%(xi.p_ntpcov_dr);
p.hhcscaleu5      =  exp(x(xi.p_hhcscaleu5));%(xi.p_hhcscaleu5);
p.hhcscale        =  exp(x(xi.p_hhcscale));%(xi.p_hhcscale);


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

p.tx_init       =  p.Tx_init.*[1 1];
p.tx_init2       =  p.Tx_init2.*[1 1];

%% 
r.hhc_scale = [p.hhcscaleu5 ,p.hhcscale ,p.hhcscale ,p.hhcscale ];

% Construct vector for progression from 'fast' latent infection, each
% element representing a different HIV status: 1.HIV-ve, 2.HIV+ve, 3.On ART
tmp = r.progression0(1)*[1 1 1;1 1 1];
tmp(:,[2,3]) = tmp(:,[2,3])*r.RRhiv .*[1 r.ARTred ;1 r.ARTred ];
tmp(2,:)=tmp(2,:).*r.RRslum;
r.progression = tmp;

% Construct vector for progression from 'slow' latent infection, each
% element representing a different HIV status: 1.HIV-ve, 2.HIV+ve, 3.On ART
tmp = r.reactivation0(1)*[1 1 1;1 1 1];
tmp(:,[2,3]) = tmp(:,[2,3])*r.RRhiv .*[1 r.ARTred ;1 r.ARTred];
tmp(2,:)=tmp(2,:).*r.RRslum;
r.reactivation = tmp;

r.self_cure([1,3])    = r.self_cure0([1,3]);
% p.imm                 = p_imm * [1 0 1];

