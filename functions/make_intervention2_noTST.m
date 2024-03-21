%enhanced intervention

function M = make_intervention2_noTST(p, r, i,s, gps,intvn)



% --- Baseline state
M0 = make_model3(p, r, i, s, gps);
Mset=[];
% for iseq = 1:length(intvn_seq)
%     intvn = intvn_seq(iseq);
p.tst_spec=0;
p.tst_sens=1;

if intvn == 0
    p0=p;
    r0=r;
    p0.sl_short=p.sl_short_itv;
    p0.tpt_short=p.tpt_short_itv;

    p0.cfy =(p.house_size-1);
    p0.hhc_a_cov_tpt=[p.hhcovu5,p.hhcov,p.hhcov,p.hhcov];
    p0.hhc_a_cov_screen=[p.hhcovu5,p.hhcov,p.hhcov,p.hhcov];

    %% PLHIV
    % Newly starting ART
    p0.case_findingPLHIV_U(1) =p.hiv_a0_baseline_noart_spec;
    p0.case_findingPLHIV_U(2) =p.hiv_a5_baseline_noart_spec;
    p0.case_findingPLHIV_U(3) =p.hiv_a10_baseline_noart_spec;
    p0.case_findingPLHIV_U(4) =p.hiv_a15_baseline_noart_spec;

    p0.case_findingPLHIV_Lf(1) =p.hiv_a0_baseline_noart_spec;
    p0.case_findingPLHIV_Lf(2) =p.hiv_a5_baseline_noart_spec;
    p0.case_findingPLHIV_Lf(3) =p.hiv_a10_baseline_noart_spec;
    p0.case_findingPLHIV_Lf(4) =p.hiv_a15_baseline_noart_spec;

    p0.case_findingPLHIV_Ls(1) =p.hiv_a0_baseline_noart_spec;
    p0.case_findingPLHIV_Ls(2) =p.hiv_a5_baseline_noart_spec;
    p0.case_findingPLHIV_Ls(3) =p.hiv_a10_baseline_noart_spec;
    p0.case_findingPLHIV_Ls(4) =p.hiv_a15_baseline_noart_spec;

    % HHC
    p0.case_findingHHC_U(1,1,1) = p.hhc_a0_baseline_spec;
    p0.case_findingHHC_U(2,1,1) = p.hhc_a5_baseline_spec*(1-p.tst_spec);
    p0.case_findingHHC_U(3,1,1) = p.hhc_a10_baseline_spec*(1-p.tst_spec);
    p0.case_findingHHC_U(4,1,1) = p.hhc_a15_baseline_spec*(1-p.tst_spec);

    p0.case_findingHHC_Lf(1,1,1) = p.hhc_a0_baseline_spec;
    p0.case_findingHHC_Lf(2,1,1) = p.hhc_a5_baseline_spec*p.tst_sens;
    p0.case_findingHHC_Lf(3,1,1) = p.hhc_a10_baseline_spec*p.tst_sens;
    p0.case_findingHHC_Lf(4,1,1) = p.hhc_a15_baseline_spec*p.tst_sens;

    p0.case_findingHHC_Ls(1,1,1) = p.hhc_a0_baseline_spec;
    p0.case_findingHHC_Ls(2,1,1) = p.hhc_a5_baseline_spec*p.tst_sens;
    p0.case_findingHHC_Ls(3,1,1) = p.hhc_a10_baseline_spec*p.tst_sens;
    p0.case_findingHHC_Ls(4,1,1) = p.hhc_a15_baseline_spec*p.tst_sens;

    p0.case_findingHHC_I(1,1,1) = p.hhc_a0_baseline_sens;
    p0.case_findingHHC_I(2,1,1) = p.hhc_a5_baseline_sens;
    p0.case_findingHHC_I(3,1,1) = p.hhc_a10_baseline_sens;
    p0.case_findingHHC_I(4,1,1) = p.hhc_a15_baseline_sens;

    Mset = make_model3(p0, r0, i, s, gps);

elseif intvn == 1
    %% PLHIV
    % Newly starting ART
    p1=p;
    r1=r;
    p1.hiv_cov_tpt=1;
    p1.hiv_cov_screen=1;
    r1.outReg = 12/3; % 3HP
    p1.potency = 1; % 3HP
    p1.sl_short=1;
    p1.tpt_short=1;
    p1.tpt_target_compa0=0.71;
    p1.tpt_target_compa15p=0.9;
    p1.cfy =(p.house_size-1);
    p1.hhc_a_cov_tpt=[p.hhcovu5,p.hhcov,p.hhcov,p.hhcov];
    p1.hhc_a_cov_screen=[p.hhcovu5,p.hhcov,p.hhcov,p.hhcov];


    %baseline HHC
     % HHC
     p1.case_findingHHC_U(1,1,1) = p.hhc_a0_baseline_spec;
     p1.case_findingHHC_U(2,1,1) = p.hhc_a5_baseline_spec*(1-p.tst_spec);
     p1.case_findingHHC_U(3,1,1) = p.hhc_a10_baseline_spec*(1-p.tst_spec);
     p1.case_findingHHC_U(4,1,1) = p.hhc_a15_baseline_spec*(1-p.tst_spec);
 
     p1.case_findingHHC_Lf(1,1,1) = p.hhc_a0_baseline_spec;
     p1.case_findingHHC_Lf(2,1,1) = p.hhc_a5_baseline_spec*p.tst_sens;
     p1.case_findingHHC_Lf(3,1,1) = p.hhc_a10_baseline_spec*p.tst_sens;
     p1.case_findingHHC_Lf(4,1,1) = p.hhc_a15_baseline_spec*p.tst_sens;
 
     p1.case_findingHHC_Ls(1,1,1) = p.hhc_a0_baseline_spec;
     p1.case_findingHHC_Ls(2,1,1) = p.hhc_a5_baseline_spec*p.tst_sens;
     p1.case_findingHHC_Ls(3,1,1) = p.hhc_a10_baseline_spec*p.tst_sens;
     p1.case_findingHHC_Ls(4,1,1) = p.hhc_a15_baseline_spec*p.tst_sens;
 
     p1.case_findingHHC_I(1,1,1) = p.hhc_a0_baseline_sens;
     p1.case_findingHHC_I(2,1,1) = p.hhc_a5_baseline_sens;
     p1.case_findingHHC_I(3,1,1) = p.hhc_a10_baseline_sens;
     p1.case_findingHHC_I(4,1,1) = p.hhc_a15_baseline_sens;

     %PLHIV
   
    p1.case_findingPLHIV_U(1) =p.hiv_a0_enhanced_noart_spec;
    p1.case_findingPLHIV_U(2) =p.hiv_a5_enhanced_noart_spec;
    p1.case_findingPLHIV_U(3) =p.hiv_a10_enhanced_noart_spec;
    p1.case_findingPLHIV_U(4) =p.hiv_a15_enhanced_noart_spec;

    p1.case_findingPLHIV_Lf(1) =p.hiv_a0_enhanced_noart_spec;
    p1.case_findingPLHIV_Lf(2) =p.hiv_a5_enhanced_noart_spec;
    p1.case_findingPLHIV_Lf(3) =p.hiv_a10_enhanced_noart_spec;
    p1.case_findingPLHIV_Lf(4) =p.hiv_a15_enhanced_noart_spec;

    p1.case_findingPLHIV_Ls(1) =p.hiv_a0_enhanced_noart_spec;
    p1.case_findingPLHIV_Ls(2) =p.hiv_a5_enhanced_noart_spec;
    p1.case_findingPLHIV_Ls(3) =p.hiv_a10_enhanced_noart_spec;
    p1.case_findingPLHIV_Ls(4) =p.hiv_a15_enhanced_noart_spec;

    % those on ART yearly screen
  
    p1.case_finding_U(1,1,3) = p.hiv_a0_enhanced_art_spec;
    p1.case_finding_U(2,1,3) = p.hiv_a5_enhanced_art_spec;
    p1.case_finding_U(3,1,3) = p.hiv_a10_enhanced_art_spec;
    p1.case_finding_U(4,1,3) = p.hiv_a15_enhanced_art_spec;

    p1.case_finding_Lf(1,1,3) = p.hiv_a0_enhanced_art_spec;
    p1.case_finding_Lf(2,1,3) = p.hiv_a5_enhanced_art_spec;
    p1.case_finding_Lf(3,1,3) = p.hiv_a10_enhanced_art_spec;
    p1.case_finding_Lf(4,1,3) = p.hiv_a15_enhanced_art_spec;

    p1.case_finding_Ls(1,1,3) = p.hiv_a0_enhanced_art_spec;
    p1.case_finding_Ls(2,1,3) = p.hiv_a5_enhanced_art_spec;
    p1.case_finding_Ls(3,1,3) = p.hiv_a10_enhanced_art_spec;
    p1.case_finding_Ls(4,1,3) = p.hiv_a15_enhanced_art_spec;

    p1.case_finding_I(1,1,3) = p.hiv_a0_enhanced_art_sens;
    p1.case_finding_I(2,1,3) = p.hiv_a5_enhanced_art_sens;
    p1.case_finding_I(3,1,3) = p.hiv_a10_enhanced_art_sens;
    p1.case_finding_I(4,1,3) = p.hiv_a15_enhanced_art_sens;


    Mset = make_model3(p1, r1, i, s, gps);

elseif intvn == 2
    %% HHC
  
    p2=p;
    r2=r;
    r2.outReg = 12/3; % 3HP
    p2.potency = 1; % 3HP
    p2.sl_short=1;
    p2.tpt_short=1;
    p2.tpt_target_compa0=0.71;
    p2.tpt_target_compa15p=0.9;
    p2.cfy =(p.house_size-1);
    p2.hhc_a_cov_tpt=[1,1,1,1];
    p2.hhc_a_cov_screen=[1,1,1,1];

%Baseline PLHIV

    p2.case_findingPLHIV_U(1) =p.hiv_a0_baseline_noart_spec;
    p2.case_findingPLHIV_U(2) =p.hiv_a5_baseline_noart_spec;
    p2.case_findingPLHIV_U(3) =p.hiv_a10_baseline_noart_spec;
    p2.case_findingPLHIV_U(4) =p.hiv_a15_baseline_noart_spec;

    p2.case_findingPLHIV_Lf(1) =p.hiv_a0_baseline_noart_spec;
    p2.case_findingPLHIV_Lf(2) =p.hiv_a5_baseline_noart_spec;
    p2.case_findingPLHIV_Lf(3) =p.hiv_a10_baseline_noart_spec;
    p2.case_findingPLHIV_Lf(4) =p.hiv_a15_baseline_noart_spec;

    p2.case_findingPLHIV_Ls(1) =p.hiv_a0_baseline_noart_spec;
    p2.case_findingPLHIV_Ls(2) =p.hiv_a5_baseline_noart_spec;
    p2.case_findingPLHIV_Ls(3) =p.hiv_a10_baseline_noart_spec;
    p2.case_findingPLHIV_Ls(4) =p.hiv_a15_baseline_noart_spec;

% HHC
    p2.case_findingHHC_U(1,1,1) = p.hhc_a0_enhanced_spec;
    p2.case_findingHHC_U(2,1,1) = p.hhc_a5_enhanced_spec*(1-p.tbst_spec);
    p2.case_findingHHC_U(3,1,1) = p.hhc_a10_enhanced_spec*(1-p.tbst_spec);
    p2.case_findingHHC_U(4,1,1) = p.hhc_a15_enhanced_spec*(1-p.tbst_spec);

    p2.case_findingHHC_Lf(1,1,1) = p.hhc_a0_enhanced_spec;
    p2.case_findingHHC_Lf(2,1,1) = p.hhc_a5_enhanced_spec*p.tbst_sens;
    p2.case_findingHHC_Lf(3,1,1) = p.hhc_a10_enhanced_spec*p.tbst_sens;
    p2.case_findingHHC_Lf(4,1,1) = p.hhc_a15_enhanced_spec*p.tbst_sens;

    p2.case_findingHHC_Ls(1,1,1) = p.hhc_a0_enhanced_spec;
    p2.case_findingHHC_Ls(2,1,1) = p.hhc_a5_enhanced_spec*p.tbst_sens;
    p2.case_findingHHC_Ls(3,1,1) = p.hhc_a10_enhanced_spec*p.tbst_sens;
    p2.case_findingHHC_Ls(4,1,1) = p.hhc_a15_enhanced_spec*p.tbst_sens;

    p2.case_findingHHC_I(1,1,1) = p.hhc_a0_enhanced_sens;
    p2.case_findingHHC_I(2,1,1) = p.hhc_a5_enhanced_sens;
    p2.case_findingHHC_I(3,1,1) = p.hhc_a10_enhanced_sens;
    p2.case_findingHHC_I(4,1,1) = p.hhc_a15_enhanced_sens;


    Mset = make_model3(p2, r2, i, s, gps);


elseif intvn == 3
    %% PLHIV + HHC
    p3=p;
    r3=r;
    r3.outReg = 12/3; % 3HP
    p3.potency = 1; % 3HP
    p3.sl_short=1;
    p3.tpt_short=1;
    p3.tpt_target_compa0=0.71;
    p3.tpt_target_compa15p=0.9;
    p3.cfy =(p.house_size-1);
    p3.hhc_a_cov_tpt=[1,1,1,1];
    p3.hhc_a_cov_screen=[1,1,1,1];
    p3.hiv_cov_tpt=1;
    p3.hiv_cov_screen=1;

    %% PLHIV
    % Newly starting ART
    p3.case_findingPLHIV_U(1) =p.hiv_a0_enhanced_noart_spec;
    p3.case_findingPLHIV_U(2) =p.hiv_a5_enhanced_noart_spec;
    p3.case_findingPLHIV_U(3) =p.hiv_a10_enhanced_noart_spec;
    p3.case_findingPLHIV_U(4) =p.hiv_a15_enhanced_noart_spec;

    p3.case_findingPLHIV_Lf(1) =p.hiv_a0_enhanced_noart_spec;
    p3.case_findingPLHIV_Lf(2) =p.hiv_a5_enhanced_noart_spec;
    p3.case_findingPLHIV_Lf(3) =p.hiv_a10_enhanced_noart_spec;
    p3.case_findingPLHIV_Lf(4) =p.hiv_a15_enhanced_noart_spec;

    p3.case_findingPLHIV_Ls(1) =p.hiv_a0_enhanced_noart_spec;
    p3.case_findingPLHIV_Ls(2) =p.hiv_a5_enhanced_noart_spec;
    p3.case_findingPLHIV_Ls(3) =p.hiv_a10_enhanced_noart_spec;
    p3.case_findingPLHIV_Ls(4) =p.hiv_a15_enhanced_noart_spec;

    % those on ART yearly screen
    p3.case_finding_U(1,1,3) = p.hiv_a0_enhanced_art_spec;
    p3.case_finding_U(2,1,3) = p.hiv_a5_enhanced_art_spec;
    p3.case_finding_U(3,1,3) = p.hiv_a10_enhanced_art_spec;
    p3.case_finding_U(4,1,3) = p.hiv_a15_enhanced_art_spec;

    p3.case_finding_Lf(1,1,3) = p.hiv_a0_enhanced_art_spec;
    p3.case_finding_Lf(2,1,3) = p.hiv_a5_enhanced_art_spec;
    p3.case_finding_Lf(3,1,3) = p.hiv_a10_enhanced_art_spec;
    p3.case_finding_Lf(4,1,3) = p.hiv_a15_enhanced_art_spec;

    p3.case_finding_Ls(1,1,3) = p.hiv_a0_enhanced_art_spec;
    p3.case_finding_Ls(2,1,3) = p.hiv_a5_enhanced_art_spec;
    p3.case_finding_Ls(3,1,3) = p.hiv_a10_enhanced_art_spec;
    p3.case_finding_Ls(4,1,3) = p.hiv_a15_enhanced_art_spec;

    p3.case_finding_I(1,1,3) = p.hiv_a0_enhanced_art_sens;
    p3.case_finding_I(2,1,3) = p.hiv_a5_enhanced_art_sens;
    p3.case_finding_I(3,1,3) = p.hiv_a10_enhanced_art_sens;
    p3.case_finding_I(4,1,3) = p.hiv_a15_enhanced_art_sens;

    % HHC
    p3.case_findingHHC_U(1,1,1) = p.hhc_a0_enhanced_spec;
    p3.case_findingHHC_U(2,1,1) = p.hhc_a5_enhanced_spec*(1-p.tbst_spec);
    p3.case_findingHHC_U(3,1,1) = p.hhc_a10_enhanced_spec*(1-p.tbst_spec);
    p3.case_findingHHC_U(4,1,1) = p.hhc_a15_enhanced_spec*(1-p.tbst_spec);

    p3.case_findingHHC_Lf(1,1,1) = p.hhc_a0_enhanced_spec;
    p3.case_findingHHC_Lf(2,1,1) = p.hhc_a5_enhanced_spec*p.tbst_sens;
    p3.case_findingHHC_Lf(3,1,1) = p.hhc_a10_enhanced_spec*p.tbst_sens;
    p3.case_findingHHC_Lf(4,1,1) = p.hhc_a15_enhanced_spec*p.tbst_sens;

    p3.case_findingHHC_Ls(1,1,1) = p.hhc_a0_enhanced_spec;
    p3.case_findingHHC_Ls(2,1,1) = p.hhc_a5_enhanced_spec*p.tbst_sens;
    p3.case_findingHHC_Ls(3,1,1) = p.hhc_a10_enhanced_spec*p.tbst_sens;
    p3.case_findingHHC_Ls(4,1,1) = p.hhc_a15_enhanced_spec*p.tbst_sens;

    p3.case_findingHHC_I(1,1,1) = p.hhc_a0_enhanced_sens;
    p3.case_findingHHC_I(2,1,1) = p.hhc_a5_enhanced_sens;
    p3.case_findingHHC_I(3,1,1) = p.hhc_a10_enhanced_sens;
    p3.case_findingHHC_I(4,1,1) = p.hhc_a15_enhanced_sens;


    Mset = make_model3(p3, r3, i, s, gps);


elseif intvn == 4
    %% PLHIV+HHC+High risk
    %% PLHIV + HHC
    p4=p;
    r4=r;
    r4.outReg = 12/3; % 3HP
    p4.potency = 1; % 3HP
    p4.sl_short=1;
    p4.tpt_short=1;
    p4.tpt_target_compa0=0.71;
    p4.tpt_target_compa15p=0.9;
    p4.cfy =(p.house_size-1);
    p4.hhc_a_cov_tpt=[1,1,1,1];
    p4.hhc_a_cov_screen=[1,1,1,1];
    p4.hiv_cov_tpt=1;
    p4.hiv_cov_screen=1;
    p4.slum_enhance=1;
    p4.slum_cov_tpt=0.9;
    p4.slum_cov_screen=0.9;


        %% PLHIV
    % Newly starting ART
    p4.case_findingPLHIV_U(1) =p.hiv_a0_enhanced_noart_spec;
    p4.case_findingPLHIV_U(2) =p.hiv_a5_enhanced_noart_spec;
    p4.case_findingPLHIV_U(3) =p.hiv_a10_enhanced_noart_spec;
    p4.case_findingPLHIV_U(4) =p.hiv_a15_enhanced_noart_spec;

    p4.case_findingPLHIV_Lf(1) =p.hiv_a0_enhanced_noart_spec;
    p4.case_findingPLHIV_Lf(2) =p.hiv_a5_enhanced_noart_spec;
    p4.case_findingPLHIV_Lf(3) =p.hiv_a10_enhanced_noart_spec;
    p4.case_findingPLHIV_Lf(4) =p.hiv_a15_enhanced_noart_spec;

    p4.case_findingPLHIV_Ls(1) =p.hiv_a0_enhanced_noart_spec;
    p4.case_findingPLHIV_Ls(2) =p.hiv_a5_enhanced_noart_spec;
    p4.case_findingPLHIV_Ls(3) =p.hiv_a10_enhanced_noart_spec;
    p4.case_findingPLHIV_Ls(4) =p.hiv_a15_enhanced_noart_spec;

    % those on ART yearly screen
    p4.case_finding_U(1,1,3) = p.hiv_a0_enhanced_art_spec;
    p4.case_finding_U(2,1,3) = p.hiv_a5_enhanced_art_spec;
    p4.case_finding_U(3,1,3) = p.hiv_a10_enhanced_art_spec;
    p4.case_finding_U(4,1,3) = p.hiv_a15_enhanced_art_spec;

    p4.case_finding_Lf(1,1,3) = p.hiv_a0_enhanced_art_spec;
    p4.case_finding_Lf(2,1,3) = p.hiv_a5_enhanced_art_spec;
    p4.case_finding_Lf(3,1,3) = p.hiv_a10_enhanced_art_spec;
    p4.case_finding_Lf(4,1,3) = p.hiv_a15_enhanced_art_spec;

    p4.case_finding_Ls(1,1,3) = p.hiv_a0_enhanced_art_spec;
    p4.case_finding_Ls(2,1,3) = p.hiv_a5_enhanced_art_spec;
    p4.case_finding_Ls(3,1,3) = p.hiv_a10_enhanced_art_spec;
    p4.case_finding_Ls(4,1,3) = p.hiv_a15_enhanced_art_spec;

    p4.case_finding_I(1,1,3) = p.hiv_a0_enhanced_art_sens;
    p4.case_finding_I(2,1,3) = p.hiv_a5_enhanced_art_sens;
    p4.case_finding_I(3,1,3) = p.hiv_a10_enhanced_art_sens;
    p4.case_finding_I(4,1,3) = p.hiv_a15_enhanced_art_sens;

    % HHC
    p4.case_findingHHC_U(1,1,1) = p.hhc_a0_enhanced_spec;
    p4.case_findingHHC_U(2,1,1) = p.hhc_a5_enhanced_spec*(1-p.tbst_spec);
    p4.case_findingHHC_U(3,1,1) = p.hhc_a10_enhanced_spec*(1-p.tbst_spec);
    p4.case_findingHHC_U(4,1,1) = p.hhc_a15_enhanced_spec*(1-p.tbst_spec);

    p4.case_findingHHC_Lf(1,1,1) = p.hhc_a0_enhanced_spec;
    p4.case_findingHHC_Lf(2,1,1) = p.hhc_a5_enhanced_spec*p.tbst_sens;
    p4.case_findingHHC_Lf(3,1,1) = p.hhc_a10_enhanced_spec*p.tbst_sens;
    p4.case_findingHHC_Lf(4,1,1) = p.hhc_a15_enhanced_spec*p.tbst_sens;

    p4.case_findingHHC_Ls(1,1,1) = p.hhc_a0_enhanced_spec;
    p4.case_findingHHC_Ls(2,1,1) = p.hhc_a5_enhanced_spec*p.tbst_sens;
    p4.case_findingHHC_Ls(3,1,1) = p.hhc_a10_enhanced_spec*p.tbst_sens;
    p4.case_findingHHC_Ls(4,1,1) = p.hhc_a15_enhanced_spec*p.tbst_sens;

    p4.case_findingHHC_I(1,1,1) = p.hhc_a0_enhanced_sens;
    p4.case_findingHHC_I(2,1,1) = p.hhc_a5_enhanced_sens;
    p4.case_findingHHC_I(3,1,1) = p.hhc_a10_enhanced_sens;
    p4.case_findingHHC_I(4,1,1) = p.hhc_a15_enhanced_sens;

    %High risk
    p4.case_findingslum_U(4,2,1) = p.slum_a15_enhanced_spec*(1-p.tbst_spec);
    p4.case_findingslum_Lf(4,2,1) = p.slum_a15_enhanced_spec*p.tbst_sens;
    p4.case_findingslum_Ls(4,2,1) = p.slum_a15_enhanced_spec*p.tbst_sens;
    p4.case_findingslum_I(4,2,1) = p.slum_a15_enhanced_sens;


    Mset = make_model3(p4, r4, i, s, gps);



end
M = Mset;



