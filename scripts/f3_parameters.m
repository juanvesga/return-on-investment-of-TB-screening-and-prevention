
%% Free parameters

xi = [];
xnames ={...
    'r_beta',...
    'rf_beta_mdr',...
    'rf_mort_TB',...
    'rf_mort_TBh',...
    'r_RRslum',...
    'r_RRhiv',...
    'rf_careseeking_mo',...
    'p_Tx_init',...
    'p_Tx_init2',...
    'r_ART_init',...
    'p_IPThiv',...
    'p_ntpcov',...
    'p_ntpcov_dr',...
    'p_hhcovu5', ...
    'p_hhcov'};


xnums  = ones(1,size(xnames,2));
lim = 0;
for ii = 1:length(xnames)
    inds = lim + (1:xnums(ii));
    xi.(xnames{ii}) = inds;
    lim = inds(end);
end
bds=zeros(length(xnames),2);
xi.nx = lim;
bds(xi.r_beta,:)          = [1 25];
bds(xi.rf_beta_mdr,:)     = [0.8 1.2];
bds(xi.rf_mort_TB,:)      = [0.4 1.6];
bds(xi.rf_mort_TBh,:)     = [1 10];
bds(xi.r_RRslum,:)        = [1, 70];
bds(xi.r_RRhiv,:)         = [1 40];  
bds(xi.rf_careseeking_mo,:)= [1 12];  
bds(xi.p_Tx_init,:)       = [0 1];
bds(xi.p_Tx_init2,:)       = [0 1];
bds(xi.r_ART_init,:)      = [0 5];
bds(xi.p_IPThiv,:)        = [0 1];
bds(xi.p_ntpcov,:)        = [0 1];
bds(xi.p_ntpcov_dr,:)     = [0 1];
bds(xi.p_hhcovu5,:)       = [0 0.2];
bds(xi.p_hhcov,:)         = [0 0.2];

prm.bds_pure=bds';
p.scale=max(diff(prm.bds_pure)/2)./(diff(prm.bds_pure)/2);
prm.bds = bds'.*p.scale;



% Haario algo params

params_haario = {
    %      name,  init,        min, max, mu,  sig, target?, local?
    {'r_beta',           0 , prm.bds(1,1), prm.bds(2,1)}
    {'rf_beta_mdr',       0 , prm.bds(1,2), prm.bds(2,2)}
    {'rf_mort_TB',       0 , prm.bds(1,3), prm.bds(2,3)}
    {'rf_mort_TBh',      0 , prm.bds(1,4), prm.bds(2,4)}
    {'r_RRslum',         0 , prm.bds(1,5), prm.bds(2,5)}
    {'r_RRhiv',          0 , prm.bds(1,6), prm.bds(2,6)}
    {'rf_careseeking_mo',0 , prm.bds(1,7), prm.bds(2,7)}
    {'p_Tx_init',        0 , prm.bds(1,8), prm.bds(2,8)}
    {'p_Tx_init2',        0 , prm.bds(1,9), prm.bds(2,9)}
    {'r_ART_init',       0 , prm.bds(1,10), prm.bds(2,10)}
    {'p_IPThiv',         0 , prm.bds(1,11), prm.bds(2,11)}
    {'p_ntpcov',      0 , prm.bds(1,12), prm.bds(2,12)}
    {'p_ntpcov_dr',      0 , prm.bds(1,13), prm.bds(2,13)}
    {'p_hhcovu5',        0 , prm.bds(1,14), prm.bds(2,14)}
    {'p_hhcov',          0 , prm.bds(1,15), prm.bds(2,15)}};

%% Fixed Parameters

%load parameters
pars  =load_data('data/params.xlsx');
reg_profile=load_data('data/reg_profiles.xlsx');
regimen="6IPT";  %for calibration purposes

% --- Define baseline parameter values, organising all parameters as
% per-capita rates (r) or proportions (p)

% Natural history (fast, slow and fast to slow from Emery et al 2021)
r.progression0  = [ 0.086*[1 10 10*0.4] ;0.0826*[1 10 10*0.4]];
r.reactivation0 = [6.2228e-04*[1 10 10*0.4];0.0006*[1 10 10*0.4]];
r.self_cure0    = 1/6*[1 0 1];
r.mort_TB0      = 1/6;
r.LTBI_stabil  = [0.92*[1 0 1];0.872*[1 0 1];];
r.relapse      = [0.14 0.0015]; % High and stable only
r.mort          = 1/66;
p.imm           = [0.77 0 0.77];

% Treatment stage
p.covid_disrupt=0.3;
p.fl_fail = pars{country,'fl_fail'};
p.fl_lost = pars{country,'fl_lost'};
p.fl_die  = pars{country,'fl_die'};
p.fl_suc  =  1-(p.fl_fail+p.fl_lost+p.fl_die);%

p.sl_short = pars{country,'sl_short'};
r.MDR_acqu = 0.02;
p.SL_trans = 0.8; %Transfer from FL


% demographics and mixing
p.popn_turnover = 1;
p.pop_vec=pop_vec;
p.pop =sum(pop_vec);
p.pop1970 = pars{country,'pop1970'};
p.house_size = pars{country,'house_size'};

p.lex    = pars{country,'lex'};% Life expectancy
r.ageing = [1/5,1/5,1/6];          % Rate of ageing into >15 (tweaked to match simulated population fraction)
p.frac_paed=0.4089;       % Population fraction <15 2016;
p.resp_symptomatic=0.045; % Fraction Respiratoty symnptomatic (From kenya survey (eligible only by symptoms)(TB survey)
p.growth            = pars{country,'growth'};
p.frac_pop          = pop_vec/sum(pop_vec);        % Population fraction <15 2016;
if strcmp(country,'GEO')
    p.frac_slum =      17556/pop_vec(4);
elseif strcmp(country,'BRA')
    p.frac_slum =      0.003;
elseif strcmp(country,'KEN')
    n_slum=686177;
    n_slum15=416243;

    p.frac_slum = n_slum/sum(pop_vec);

    dist=p.frac_pop(1:3)./(1-p.frac_pop(4));

    p.frac_pop_slum=[...
        dist(1)*(1-( n_slum15/n_slum)),...
        dist(2)*(1-( n_slum15/n_slum)),...
        dist(3)*(1-( n_slum15/n_slum)),...
        ( n_slum15/n_slum)];

elseif strcmp(country,'ZAF')

    n_slum=1554870;
    n_slum15=1109090;

    p.frac_slum = n_slum/sum(pop_vec);

    dist=p.frac_pop(1:3)./(1-p.frac_pop(4));

    p.frac_pop_slum=[...
        dist(1)*(1-( n_slum15/n_slum)),...
        dist(2)*(1-( n_slum15/n_slum)),...
        dist(3)*(1-( n_slum15/n_slum)),...
        ( n_slum15/n_slum)];

end

p.bra=0;
p.geo=0;
if strcmp(country,'BRA')
    p.entry_highrisk =      200/1e5;
    p.out_highrisk =      12/24;% AVerage prison time of 2 years
    p.geo_mix= [1,0;0,1];
    p.bra=1;
elseif strcmp(country,'GEO')
    p.entry_highrisk =      0;
    p.out_highrisk =      0;
    p.geo_mix= [1,0.3;0.3,1];
    p.geo=1;
else
    p.entry_highrisk =      0;
    p.out_highrisk =      0;
    p.geo_mix= [1,0.3;0.3,1];
end

% if strcmp(country,'BRA') || strcmp(country,'ZAF')
%     p.birth_slum=1;
% else
p.birth_slum=1;
% end

p.cross = 0;              % Average contact pattern withlike for <15 and >15 In Kenya (ChapaKiti et al 2014,PlosOne)
r.RRsize=0;               % Factor reducing careseeking amongs resp symptomatic
r.RRage=[1 1];            % If needed to vary careseeking by age

% Natural history
r.beta = 0;
r.beta_mdr=0;
p.kappa     = 1;
r.RRslum    = 1;
r.self_clearance = 0;% 0.024;%Emery et al 2012 https://doi.org/10.1098/rspb.2020.1635
r.selfcure      =  1/3*0.5; % Tiemersma
r.muTB          =  1/3*0.5; % Tiemersma
r.RRmortTBhiv   = 3;
r.mort_ht= 1/36;          % Non-TB-related mortality in treated HIV+ (assumption: ex6tebeds live as far as observed, 2017-1982)
r.RRTBtxmort = (1-0.92);        % 92% red in TTx Based on: Borgdorff MW, Floyd K, Broekmans JF. Interventions to reduce tuberculosis mortality and transmission in low- and middle-income countries: effectiveness,                          % cost-effectiveness, and constraints to scaling up. (CMH Working Paper Series, Paper No. WG5: 8. Available at: URL: www.cmhealth.org/wg5_paper8.pdf
r.mort   = [1/p.lex  , 1/(p.lex-5), 1/(p.lex-10), 1/(p.lex-(10+30)) ];    % Non-disease-related mortality
r.mort_h = 1/15;          % Non-TB-related mortality in untreated HIV+

r.RRprog_age = [0.2,0.3,0.4,1];
r.fast_react      =0.0826.*r.RRprog_age;%Menxies [2.4090 ,  0.9855 , 0.0985 , 0.0985 ];      % Ragonnet 2018_Epidemics
r.slow_react      =0.000594.*r.RRprog_age;%Menzies [6.9350e-09 ,   0.0023 , 0.0012 , 0.0012 ]; % Ragonnet 2018_Epidemics
r.slow            =0.8728;%Menzies   [4.38, 4.38 , 1.9710 , 1.9710];             % Ragonnet 2018_Epidemics
p.crossg    = 0.3;
r.symp_del = [12/4 12/1.3 0.35*12/1.3]; %Bossen et al lancet Rmed 2023 https://www.thelancet.com/action/showPdf?pii=S2213-2600%2823%2900097-8
r.bRRpaed  =0;


% Pre care & Health System
r.careseek=0;
r.careseekRR_h=0;
r.careseeking2=12;
r.cs2           = 2;
r.csRRhiv       = 3;% Relative in crease in careseeking for HIV + on treatemt
r.Dx = 52;        %1 week between seek care and set on treatmkent
p.Dx = 0.8;%Probability of Dx attempt (availability of services at point of service) by provider and HIv status
p.smear_sens=0.8; %Smear test speccificity (Swai F 2011, BMC resreach notes)
p.smear_spec=0.94;%Smear test speccificity (Swai F 2011, BMC resreach notes)
p.xpert_sens=0.9; %Smear test speccificity (Swai F 2011, BMC resreach notes)
p.xpert_spec=0.99;%Smear test speccificity (Swai F 2011, BMC resreach notes)
p.smear=0.5;%prob of microscopy by sector
p.xpert=pars{country,'xpert'};%prob of microscopy by sector
p.yrntp = pars{country,'yrntp'};% care seeking 0.18pry as estimated by Olney et al 2016
p.ART_start = pars{country,'yrhaart'};



%First line

% r.TxOut   = r.Tx/(p.fl_suc+p.fl_fail); % total rate of treatment outcomes
% r.default = r.TxOut*p.fl_lost;
% r.fl_mort = r.TxOut*p.fl_die;
% p.cure    = p.fl_suc/(p.fl_suc+p.fl_fail);
%
% %Second line
%
% r.TxOut2   = r.Tx2/(p.sl_suc+p.sl_fail); % total rate of treatment outcomes
% r.default2 = r.TxOut2*p.sl_lost;
% r.sl_mort = r.TxOut2*p.sl_die;
% p.cure2    = p.sl_suc/(p.sl_suc+p.sl_fail);

p.IPT=zeros(1,3) ;%ia,ig
p.LTBItest  = 0;% set to 1 to set on PT only latent infections
p.housedist = [0.7 0.2 0.1];
p.rif_free = 1; % Fraction f TPT that is rifampn free
r.ptdur           =  reg_profile{regimen,'ptdur'};    % Duration of protection post regimen
r.ptregdur        =  reg_profile{regimen,'ptregdur'};     % Regimen duration (halved to show first and 2nd half)
p.pteffi          =  reg_profile{regimen,'pteffi'};     % efficacy
p.ptdrbarrier     =  reg_profile{regimen,'ptdrbarr'};     % DR Barrier
p.ptforg          =  reg_profile{regimen,'ptforg'};     % Regimen Forgiveness
p.ptcomp          =  reg_profile{regimen,'ptcomp'};     % Regimen completion (ease of adherence
p.nomido          =  reg_profile{regimen,'ptnomido'};    %
p.midoforg        =  reg_profile{regimen,'ptmidoforg'};    %
p.potency         =  reg_profile{regimen,'potency'};    %
p.reinfprotection =  1;
r.corr_hiv  = 1;%  0.75; % ratio of efficacy HIV+/HIV-

r.regdur= r.ptregdur;
p.tpt_short = pars{country,'tpt_short'};
p.tpt_base_comp = 0.63;
p.tpt_target_comp = 0.81;
p.tpt_target_compa0 = 0.81; % change to 0.71 in any intervention
p.tpt_target_compa15p = 0.81; % set to 90% in enhanced in >15

p.effondr       =  0; % Effect of PT on Dr strains
p.prevComm=1;

p.notified = p.pop*notifiedrate(end)/1e5;



% HIV cascades
r.bRRhiv    = 0.7; % realtive TB infectiousness of HIV+
r.progRRhiv = 26; % Getahun , Selwyn
p.hivtest_notb = pars{country,'hivtest_notb'};% care seeking 0.18pry as estimated by Olney et al 2016
p.hivtest_tb   = pars{country,'hivtest_tb'}  ;% WHO report % known status.. care seeking inside TB care system
p.tbtest_hiv   = pars{country,'tbtest_hiv'}  ;% Tb test among HIV pos
p.art_notb     = pars{country,'art_notb'}   ;% Coverage expected in 2015 * proportion estimated to adhreing and getting suppressed
p.art_tb       = pars{country,'art_tb'}   ;% Coverage expected in TB+ 2015 * proportion estimated to adhreing and getting suppressed
r.art_dropout  = pars{country,'art_dropout'};% Average UNGAS2016 + AMPATH/
r.hivdecline   = hivdecline;%pars{country,'hivdecline'};    % Annual Decline in HIV incidence after 2017
r.reinfhiv=1; % Susceptibility to reinfection on HIV+
r.ARTred    = 0.35;        % Reduction in progression given ART (Suthar2010PlosMed)
r.ARTrec    = 0;           % Recruitment probability

if strcmp(country,'GEO')
    p.hivcq_ad      = [0.1 0.1 0.1 2];%[0.05 0.13 1 1]; % Ratio to waight HIV rate in children and adults
elseif strcmp(country,'BRA')
    p.hivcq_ad      = [0.1 0.2 0.5 2.2];%[0.05 0.13 1 1]; % Ratio to waight HIV rate in children and adults
elseif strcmp(country,'KEN')
    p.hivcq_ad      = [0.5 0.5 0.5 1.5];%[0.05 0.13 1 1]; % Ratio to waight HIV rate in children and adults
elseif strcmp(country,'ZAF')
    p.hivcq_ad      = [0.1 0.2 0.5 4.5];
end


% Interventions
r.OptimART = 1;           % Parameter to optimize HIV cascade in intervention
r.tent=0;                 % Control parameter for tent-shape scale up of IPT catch up
p.ptcure=0;               % Fraction curative efficacy of pt
r.acf_sym = [0 0 0;0 0 0;0 0 0;0 0 0];% ACF by Age, HIV status
r.acf_asym =0;            % ACF by Age, HIV status in asymptomatics
p.demand = 0.0;           % Demand generation
p.cfy_all=0;              % Switch parametyer to do contact trace
p.cfy    =0;              % Switch parametyer to do contact trace
p.prevComm=1;             % TB prevalence in the community for CFY
r.demand1=1;              % Demand generation
r.demand2=1;              % Demand generation
p.redM=1;                 % TB Mortality rate reduction

% DALYs
p.yll=[68.7-10 68.7-30]; % Life Years lost in chidrenand aults
p.weights=[0.33 1];      % Utility weights: 1. Untreated TB 2.Deceased



% Switch Model Control
p.priv_off   = 1; %(==0 cuts pass to private)
r.turnoffHIV = 1; %(==0 cuts HIV)
p.ARToff     = 1; %(==0 turns ARToff)


%% Get transmission matrix ready
tmp= zeros(numel(gps.age),numel(gps.age),numel(gps.georisk),numel(gps.georisk));
tmp2=zeros(numel(gps.age),1);
for ii = 1:length(gps.age)
    for jj = 1:length(gps.age)
        for kk = 1:length(gps.georisk)
            for ll = 1:length(gps.georisk)
                tmp(ii,jj,kk,ll) = c_mat(ii,jj) * p.geo_mix(kk,ll);
            end
        end
    end
end


p.transmission=squeeze(sum(tmp,[2,3,4]));

%% Interventions

p.sl_short_itv = pars{country,'sl_short_fut'};
p.tpt_short_itv = pars{country,'tpt_short_fut'};

p.case_finding_U=zeros(numel(gps.age),numel(gps.georisk),numel(gps.hiv));
p.case_finding_Lf=zeros(numel(gps.age),numel(gps.georisk),numel(gps.hiv));
p.case_finding_Ls=zeros(numel(gps.age),numel(gps.georisk),numel(gps.hiv));
p.case_finding_I=zeros(numel(gps.age),numel(gps.georisk),numel(gps.hiv));
p.case_findingHHC_U=zeros(numel(gps.age),numel(gps.georisk),numel(gps.hiv));
p.case_findingHHC_Lf=zeros(numel(gps.age),numel(gps.georisk),numel(gps.hiv));
p.case_findingHHC_Ls=zeros(numel(gps.age),numel(gps.georisk),numel(gps.hiv));
p.case_findingHHC_I=zeros(numel(gps.age),numel(gps.georisk),numel(gps.hiv));
p.case_findingPLHIV_U=zeros(1,numel(gps.age));
p.case_findingPLHIV_Lf=zeros(1,numel(gps.age));
p.case_findingPLHIV_Ls=zeros(1,numel(gps.age));
p.case_findingPLHIV_I=zeros(1,numel(gps.age));
p.case_findingslum_U=zeros(numel(gps.age),numel(gps.georisk),numel(gps.hiv));
p.case_findingslum_Lf=zeros(numel(gps.age),numel(gps.georisk),numel(gps.hiv));
p.case_findingslum_Ls=zeros(numel(gps.age),numel(gps.georisk),numel(gps.hiv));
p.case_findingslum_I=zeros(numel(gps.age),numel(gps.georisk),numel(gps.hiv));

%baseline algorithms
if strcmp(country,'KEN')
    p.hiv_a0_baseline_noart_sens=0.7918;
    p.hiv_a5_baseline_noart_sens=0.7918;
    p.hiv_a10_baseline_noart_sens=0.7918;
    p.hiv_a15_baseline_noart_sens=0.7918;
    p.hiv_a0_baseline_noart_spec = 0.9467;
    p.hiv_a5_baseline_noart_spec = 0.9467;
    p.hiv_a10_baseline_noart_spec = 0.9467;
    p.hiv_a15_baseline_noart_spec = 0.9467;

    p.hiv_a0_baseline_art_sens=0.4985;
    p.hiv_a5_baseline_art_sens=0.4985;
    p.hiv_a10_baseline_art_sens=0.4985;
    p.hiv_a15_baseline_art_sens=0.4985;
    p.hiv_a0_baseline_art_spec = 0.9745;
    p.hiv_a5_baseline_art_spec = 0.9745;
    p.hiv_a10_baseline_art_spec = 0.9745;
    p.hiv_a15_baseline_art_spec = 0.9745;

    p.hhc_a0_baseline_sens = 0.5461;
    p.hhc_a5_baseline_sens = 0.5461;
    p.hhc_a10_baseline_sens = 0.5461;
    p.hhc_a15_baseline_sens = 0.5461;
    p.hhc_a0_baseline_spec = 0.9992;
    p.hhc_a5_baseline_spec = 0.9992;
    p.hhc_a10_baseline_spec = 0.9992;
    p.hhc_a15_baseline_spec = 0.9992;

    p.slum_a15_baseline_sens = 0.5693;
    p.slum_a15_baseline_spec = 0.9964 ;

elseif strcmp(country,'ZAF')

    p.hiv_a0_baseline_noart_sens=0.69;
    p.hiv_a5_baseline_noart_sens=0.69;
    p.hiv_a10_baseline_noart_sens=0.69;
    p.hiv_a15_baseline_noart_sens=0.69;
    p.hiv_a0_baseline_noart_spec = 0.98;
    p.hiv_a5_baseline_noart_spec = 0.98;
    p.hiv_a10_baseline_noart_spec = 0.98;
    p.hiv_a15_baseline_noart_spec = 0.98;

    p.hiv_a0_baseline_art_sens=0.69;
    p.hiv_a5_baseline_art_sens=0.69;
    p.hiv_a10_baseline_art_sens=0.69;
    p.hiv_a15_baseline_art_sens=0.69;
    p.hiv_a0_baseline_art_spec = 0.98;
    p.hiv_a5_baseline_art_spec = 0.98;
    p.hiv_a10_baseline_art_spec = 0.98;
    p.hiv_a15_baseline_art_spec = 0.98;

    p.hhc_a0_baseline_sens = 0.4957;
    p.hhc_a5_baseline_sens = 0.4957;
    p.hhc_a10_baseline_sens = 0.4957;
    p.hhc_a15_baseline_sens = 0.6991;
    p.hhc_a0_baseline_spec = 0.9901;
    p.hhc_a5_baseline_spec = 0.9901;
    p.hhc_a10_baseline_spec = 0.9901;
    p.hhc_a15_baseline_spec = 0.9896;

    p.slum_a15_baseline_sens = 0.6396;
    p.slum_a15_baseline_spec = 0.9856 ;

elseif strcmp(country,'GEO')

    p.hiv_a0_baseline_noart_sens=0.7378;
    p.hiv_a5_baseline_noart_sens=0.7378;
    p.hiv_a10_baseline_noart_sens=0.8741;
    p.hiv_a15_baseline_noart_sens=0.8741;
    p.hiv_a0_baseline_noart_spec = 0.9555;
    p.hiv_a5_baseline_noart_spec = 0.9555;
    p.hiv_a10_baseline_noart_spec = 0.9406;
    p.hiv_a15_baseline_noart_spec = 0.9406;

    p.hiv_a0_baseline_art_sens=0.4645;
    p.hiv_a5_baseline_art_sens=0.4645;
    p.hiv_a10_baseline_art_sens=0.7882;
    p.hiv_a15_baseline_art_sens=0.7882;
    p.hiv_a0_baseline_art_spec = 0.9787;
    p.hiv_a5_baseline_art_spec = 0.9787;
    p.hiv_a10_baseline_art_spec = 0.9604;
    p.hiv_a15_baseline_art_spec = 0.9604;

    p.hhc_a0_baseline_sens = 0.7181;
    p.hhc_a5_baseline_sens = 0.7181;
    p.hhc_a10_baseline_sens = 0.7181;
    p.hhc_a15_baseline_sens = 0.8613;
    p.hhc_a0_baseline_spec = 0.9888;
    p.hhc_a5_baseline_spec = 0.9888;
    p.hhc_a10_baseline_spec = 0.9888;
    p.hhc_a15_baseline_spec = 0.9846;

    p.slum_a15_baseline_sens = 0.5693;
    p.slum_a15_baseline_spec = 0.9964;


elseif strcmp(country,'BRA')
    p.hiv_a0_baseline_noart_sens=0.7918;
    p.hiv_a5_baseline_noart_sens=0.7918;
    p.hiv_a10_baseline_noart_sens=0.7918;
    p.hiv_a15_baseline_noart_sens=0.7918;
    p.hiv_a0_baseline_noart_spec = 0.9467;
    p.hiv_a5_baseline_noart_spec = 0.9467;
    p.hiv_a10_baseline_noart_spec = 0.9467;
    p.hiv_a15_baseline_noart_spec = 0.9467;

    p.hiv_a0_baseline_art_sens=0.4985;
    p.hiv_a5_baseline_art_sens=0.4985;
    p.hiv_a10_baseline_art_sens=0.4985;
    p.hiv_a15_baseline_art_sens=0.4985;
    p.hiv_a0_baseline_art_spec = 0.9745;
    p.hiv_a5_baseline_art_spec = 0.9745;
    p.hiv_a10_baseline_art_spec = 0.9745;
    p.hiv_a15_baseline_art_spec = 0.9745;

    p.hhc_a0_baseline_sens = 0.5989;
    p.hhc_a5_baseline_sens = 0.5989;
    p.hhc_a10_baseline_sens = 0.5989;
    p.hhc_a15_baseline_sens = 0.7153;
    p.hhc_a0_baseline_spec = 0.9984;
    p.hhc_a5_baseline_spec = 0.9984;
    p.hhc_a10_baseline_spec = 0.9984;
    p.hhc_a15_baseline_spec = 0.999;

    p.slum_a15_baseline_sens = 0.5693;
    p.slum_a15_baseline_spec = 0.9964;


end

%% Scenario algorithms
p.hiv_a0_basic_noart_sens=0.7378;
p.hiv_a0_basic_art_sens = 0.4645;
p.hiv_a0_basic_noart_spec = 0.9555;
p.hiv_a0_basic_art_spec=0.9787;

p.hiv_a0_enhanced_noart_sens=0.7378;
p.hiv_a0_enhanced_art_sens = 0.4645;
p.hiv_a0_enhanced_noart_spec = 0.9555;
p.hiv_a0_enhanced_art_spec=0.9787;

p.hiv_a5_basic_noart_sens=0.7378;
p.hiv_a5_basic_art_sens = 0.4645;
p.hiv_a5_basic_noart_spec = 0.9555;
p.hiv_a5_basic_art_spec=0.9787;

p.hiv_a5_enhanced_noart_sens=0.8274;
p.hiv_a5_enhanced_art_sens = 0.4645;
p.hiv_a5_enhanced_noart_spec = 0.9499;
p.hiv_a5_enhanced_art_spec=0.9787;

p.hiv_a10_basic_noart_sens=0.6565;
p.hiv_a10_basic_art_sens = 0.4645;
p.hiv_a10_basic_noart_spec = 0.9794;
p.hiv_a10_basic_art_spec=0.9787;

p.hiv_a10_enhanced_noart_sens=0.6565;
p.hiv_a10_enhanced_art_sens = 0.7265;
p.hiv_a10_enhanced_noart_spec = 0.9794;
p.hiv_a10_enhanced_art_spec=0.968;

p.hiv_a15_basic_noart_sens=0.6565;
p.hiv_a15_basic_art_sens = 0.4645;
p.hiv_a15_basic_noart_spec = 0.9794;
p.hiv_a15_basic_art_spec=0.9787;

p.hiv_a15_enhanced_noart_sens=0.6565;
p.hiv_a15_enhanced_art_sens = 0.7265;
p.hiv_a15_enhanced_noart_spec = 0.9794;
p.hiv_a15_enhanced_art_spec=0.968;

p.hhc_a0_basic_sens = 0.7181;
p.hhc_a0_basic_spec = 0.9888;
p.hhc_a0_enhanced_sens = 0.7181;
p.hhc_a0_enhanced_spec = 0.9888;
p.hhc_a5_basic_sens = 0.7181;
p.hhc_a5_basic_spec = 0.9888;
p.hhc_a5_enhanced_sens = 0.7181;
p.hhc_a5_enhanced_spec = 0.9888;
p.hhc_a10_basic_sens = 0.7181;
p.hhc_a10_basic_spec = 0.9888;
p.hhc_a10_enhanced_sens = 0.7181;
p.hhc_a10_enhanced_spec = 0.9888;
p.hhc_a15_basic_sens = 0.8613;
p.hhc_a15_basic_spec = 0.9846;
p.hhc_a15_enhanced_sens = 0.8613;
p.hhc_a15_enhanced_spec = 0.9846;

p.slum_a15_basic_sens = 0.8613;
p.slum_a15_basic_spec = 0.9846 ;
p.slum_a15_enhanced_sens = 0.8613;
p.slum_a15_enhanced_spec = 0.9846;

p.cfy = 0;
p.hhcovu5 = 0.2;
p.hhcov = 0.2;
p.hhc_a_cov_tpt=[0,0,0,0];
p.hhc_a_cov_screen=[0,0,0,0];
p.hhc_distr=[0,0,0,0];
p.hhc_distr(4)=0.031;% I Gregory J. Fox, European Respiratory Journal 2013 41: 140-156; DOI: 10.1183/09031936.00070812
p.hhc_distr(3)=0.515*0.3; % Ls
p.hhc_distr(2)=0.515*0.7; % Lf
p.hhc_distr(1)= 1- sum(p.hhc_distr(1:3)) ; % U

p.tst_sens=0.77;
p.tst_spec=0.97;
p.tbst_sens=0.78;
p.tbst_spec=0.98;
p.hiv_cov_tpt=0;
p.hiv_cov_screen=0;
p.slum_enhance=0;
p.slum_cov_tpt=0;
p.slum_cov_screen=0;
%% parameter sets

%Pass all info to structure
prm.p = p; prm.r = r;prm.endy=endyear;prm.fity=fityear;
ref.xi = xi;

%Produce a Random set
xran=zeros(1,size(bds,1));
for u=1:size(bds,1)
    a = bds(u,1);
    b = bds(u,2);
    xran(u) = (b-a).*rand(1,1) + a;
end
prm.random_set=xran.*p.scale;

% Common sense parameter set
% xgood=[
%     0.004,... beta
%     0.002,... r.beta_mdr
%     0.5*1/3,...r.muTB
%     3,...RRmortTBhiv
%     15,... progRRhiv
%     3,... careseek
%     0.9,...Dx
%     2,... ARTrec
%     0.5,...art_cov
%     2010,...art_year
%     0.7,...IPThiv
%     0.5,...nptcovds
%     0.5].*p.scale; %ntp


clear xi lim bds xnums