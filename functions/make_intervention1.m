%Basic intervention

function M = make_intervention1(p, r, i,s, gps,intvn)



% --- Baseline state

Mset=[];



if intvn == 0

     p0=p;
    r0=r;
    p0.sl_short=p.sl_short_itv;
    p0.tpt_short(:,:)=p.tpt_short_itv;

    p0.cfy =(p.house_size-1);
    p0.hhc_a_cov_tpt   =[p.hhcovu5,p.hhcov,p.hhcov,p.hhcov];
    p0.hhc_a_cov_screen=[p.hhcovu5,p.hhcov,p.hhcov,p.hhcov];

     [r0, p0]=get_country_pars(r0,p0,p.country,'Baseline','all');


     if(strcmp(p.country,"GEO")||strcmp(p.country,"BRA"))

         hhc_slum=1;
     else

         hhc_slum=[1 2];
     end


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
    p0.case_findingHHC_U(1,:,hhc_slum) = p.hhc_a0_baseline_spec;
    p0.case_findingHHC_U(2,:,hhc_slum) = p.hhc_a5_baseline_spec*(1-p.tst_spec);
    p0.case_findingHHC_U(3,:,hhc_slum) = p.hhc_a10_baseline_spec*(1-p.tst_spec);
    p0.case_findingHHC_U(4,:,hhc_slum) = p.hhc_a15_baseline_spec*(1-p.tst_spec);

    p0.case_findingHHC_Lf(1,:,hhc_slum) = p.hhc_a0_baseline_spec;
    p0.case_findingHHC_Lf(2,:,hhc_slum) = p.hhc_a5_baseline_spec*p.tst_sens;
    p0.case_findingHHC_Lf(3,:,hhc_slum) = p.hhc_a10_baseline_spec*p.tst_sens;
    p0.case_findingHHC_Lf(4,:,hhc_slum) = p.hhc_a15_baseline_spec*p.tst_sens;

    p0.case_findingHHC_Ls(1,:,hhc_slum) = p.hhc_a0_baseline_spec;
    p0.case_findingHHC_Ls(2,:,hhc_slum) = p.hhc_a5_baseline_spec*p.tst_sens;
    p0.case_findingHHC_Ls(3,:,hhc_slum) = p.hhc_a10_baseline_spec*p.tst_sens;
    p0.case_findingHHC_Ls(4,:,hhc_slum) = p.hhc_a15_baseline_spec*p.tst_sens;

    p0.case_findingHHC_I(1,:,hhc_slum) = p.hhc_a0_baseline_sens;
    p0.case_findingHHC_I(2,:,hhc_slum) = p.hhc_a5_baseline_sens;
    p0.case_findingHHC_I(3,:,hhc_slum) = p.hhc_a10_baseline_sens;
    p0.case_findingHHC_I(4,:,hhc_slum) = p.hhc_a15_baseline_sens;
    
    Mset = make_model3(p0, r0, i, s, gps);
elseif intvn == 1
    %% PLHIV
    % Newly starting ART
    p1=p;
    r1=r;
    p1.hiv_cov_tpt=0.9;
    p1.hiv_cov_screen=0.9;
    p1.sl_short=1;
    p1.tpt_short(:,:)=p.tpt_short_itv;
    p1.tpt_short(1,3)=1;

    p1.tpt_target_compa0=0.71;
    p1.cfy =(p.house_size-1);
    p1.hhc_a_cov_tpt   =[p.hhcovu5,p.hhcov,p.hhcov,p.hhcov];
    p1.hhc_a_cov_screen=[p.hhcovu5,p.hhcov,p.hhcov,p.hhcov];

    if(strcmp(p.country,"GEO")||strcmp(p.country,"BRA"))

         hhc_slum=1;
     else

         hhc_slum=[1 2];
     end

     [r1, p1]=get_country_pars(r1,p1,p.country,'TPT_basic','plhiv');


    %baseline HHC
    % HHC
    p1.case_findingHHC_U(1,:,hhc_slum) = p.hhc_a0_baseline_spec;
    p1.case_findingHHC_U(2,:,hhc_slum) = p.hhc_a5_baseline_spec*(1-p.tst_spec);
    p1.case_findingHHC_U(3,:,hhc_slum) = p.hhc_a10_baseline_spec*(1-p.tst_spec);
    p1.case_findingHHC_U(4,:,hhc_slum) = p.hhc_a15_baseline_spec*(1-p.tst_spec);

    p1.case_findingHHC_Lf(1,:,hhc_slum) = p.hhc_a0_baseline_spec;
    p1.case_findingHHC_Lf(2,:,hhc_slum) = p.hhc_a5_baseline_spec*p.tst_sens;
    p1.case_findingHHC_Lf(3,:,hhc_slum) = p.hhc_a10_baseline_spec*p.tst_sens;
    p1.case_findingHHC_Lf(4,:,hhc_slum) = p.hhc_a15_baseline_spec*p.tst_sens;

    p1.case_findingHHC_Ls(1,:,hhc_slum) = p.hhc_a0_baseline_spec;
    p1.case_findingHHC_Ls(2,:,hhc_slum) = p.hhc_a5_baseline_spec*p.tst_sens;
    p1.case_findingHHC_Ls(3,:,hhc_slum) = p.hhc_a10_baseline_spec*p.tst_sens;
    p1.case_findingHHC_Ls(4,:,hhc_slum) = p.hhc_a15_baseline_spec*p.tst_sens;

    p1.case_findingHHC_I(1,:,hhc_slum) = p.hhc_a0_baseline_sens;
    p1.case_findingHHC_I(2,:,hhc_slum) = p.hhc_a5_baseline_sens;
    p1.case_findingHHC_I(3,:,hhc_slum) = p.hhc_a10_baseline_sens;
    p1.case_findingHHC_I(4,:,hhc_slum) = p.hhc_a15_baseline_sens;

    %PLHIV
    p1.case_findingPLHIV_U(1) =p.hiv_a0_basic_noart_spec;
    p1.case_findingPLHIV_U(2) =p.hiv_a5_basic_noart_spec;
    p1.case_findingPLHIV_U(3) =p.hiv_a10_basic_noart_spec;
    p1.case_findingPLHIV_U(4) =p.hiv_a15_basic_noart_spec;

    p1.case_findingPLHIV_Lf(1) =p.hiv_a0_basic_noart_spec;
    p1.case_findingPLHIV_Lf(2) =p.hiv_a5_basic_noart_spec;
    p1.case_findingPLHIV_Lf(3) =p.hiv_a10_basic_noart_spec;
    p1.case_findingPLHIV_Lf(4) =p.hiv_a15_basic_noart_spec;

    p1.case_findingPLHIV_Ls(1) =p.hiv_a0_basic_noart_spec;
    p1.case_findingPLHIV_Ls(2) =p.hiv_a5_basic_noart_spec;
    p1.case_findingPLHIV_Ls(3) =p.hiv_a10_basic_noart_spec;
    p1.case_findingPLHIV_Ls(4) =p.hiv_a15_basic_noart_spec;

    p1.case_findingPLHIV_I(1) =p.hiv_a0_basic_noart_sens;
    p1.case_findingPLHIV_I(2) =p.hiv_a5_basic_noart_sens;
    p1.case_findingPLHIV_I(3) =p.hiv_a10_basic_noart_sens;
    p1.case_findingPLHIV_I(4) =p.hiv_a15_basic_noart_sens;




    % those on ART yearly screen

    p1.case_finding_U(1,1,3) = p.hiv_a0_basic_art_spec;
    p1.case_finding_U(2,1,3) = p.hiv_a5_basic_art_spec;
    p1.case_finding_U(3,1,3) = p.hiv_a10_basic_art_spec;
    p1.case_finding_U(4,1,3) = p.hiv_a15_basic_art_spec;

    p1.case_finding_Lf(1,1,3) = p.hiv_a0_basic_art_spec;
    p1.case_finding_Lf(2,1,3) = p.hiv_a5_basic_art_spec;
    p1.case_finding_Lf(3,1,3) = p.hiv_a10_basic_art_spec;
    p1.case_finding_Lf(4,1,3) = p.hiv_a15_basic_art_spec;

    p1.case_finding_Ls(1,1,3) = p.hiv_a0_basic_art_spec;
    p1.case_finding_Ls(2,1,3) = p.hiv_a5_basic_art_spec;
    p1.case_finding_Ls(3,1,3) = p.hiv_a10_basic_art_spec;
    p1.case_finding_Ls(4,1,3) = p.hiv_a15_basic_art_spec;

    p1.case_finding_I(1,1,3) = p.hiv_a0_basic_art_sens;
    p1.case_finding_I(2,1,3) = p.hiv_a5_basic_art_sens;
    p1.case_finding_I(3,1,3) = p.hiv_a10_basic_art_sens;
    p1.case_finding_I(4,1,3) = p.hiv_a15_basic_art_sens;


    Mset = make_model3(p1, r1, i, s, gps);

elseif intvn == 2
    %% HHC

    p2=p;
    r2=r;
    p2.sl_short=1;
    p2.tpt_short(:,:)=p.tpt_short_itv;

    % if(strcmp(p.country,"GEO")||strcmp(p.country,"BRA")||strcmp(p.country,"KEN"))
    p2.tpt_short(1,:)=1;
    % else 
    % p2.tpt_short(:,:)=1;
    % end
    

     if(strcmp(p.country,"GEO")||strcmp(p.country,"BRA"))

         hhc_slum=1;
     else

         hhc_slum=[1 2];
     end


    p2.tpt_target_compa0=0.71;
    p2.cfy =(p.house_size-1);
    p2.hhc_a_cov_tpt   =[1,0.5,0.5,0.5];
    p2.hhc_a_cov_screen=[1,0.5,0.5,0.5];

   [r2, p2]=get_country_pars(r2,p2,p.country,'TPT_basic','hhc');


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

    p2.case_findingHHC_U(1,:,hhc_slum) = p.hhc_a0_basic_spec;
    p2.case_findingHHC_U(2,:,hhc_slum) = p.hhc_a5_basic_spec*(1-p.tst_spec);
    p2.case_findingHHC_U(3,:,hhc_slum) = p.hhc_a10_basic_spec*(1-p.tst_spec);
    p2.case_findingHHC_U(4,:,hhc_slum) = p.hhc_a15_basic_spec*(1-p.tst_spec);

    p2.case_findingHHC_Lf(1,:,hhc_slum) = p.hhc_a0_basic_spec;
    p2.case_findingHHC_Lf(2,:,hhc_slum) = p.hhc_a5_basic_spec*p.tst_sens;
    p2.case_findingHHC_Lf(3,:,hhc_slum) = p.hhc_a10_basic_spec*p.tst_sens;
    p2.case_findingHHC_Lf(4,:,hhc_slum) = p.hhc_a15_basic_spec*p.tst_sens;

    p2.case_findingHHC_Ls(1,:,hhc_slum) = p.hhc_a0_basic_spec;
    p2.case_findingHHC_Ls(2,:,hhc_slum) = p.hhc_a5_basic_spec*p.tst_sens;
    p2.case_findingHHC_Ls(3,:,hhc_slum) = p.hhc_a10_basic_spec*p.tst_sens;
    p2.case_findingHHC_Ls(4,:,hhc_slum) = p.hhc_a15_basic_spec*p.tst_sens;

    p2.case_findingHHC_I(1,:,hhc_slum) = p.hhc_a0_basic_sens;
    p2.case_findingHHC_I(2,:,hhc_slum) = p.hhc_a5_basic_sens;
    p2.case_findingHHC_I(3,:,hhc_slum) = p.hhc_a10_basic_sens;
    p2.case_findingHHC_I(4,:,hhc_slum) = p.hhc_a15_basic_sens;


    Mset = make_model3(p2, r2, i, s, gps);


elseif intvn == 3
    %% PLHIV + HHC
    p3=p;
    r3=r;
    p3.sl_short=1;
    p3.tpt_short(:,:)=p.tpt_short_itv;


    % if(strcmp(p.country,"GEO")||strcmp(p.country,"BRA")||strcmp(p.country,"KEN"))
    p3.tpt_short(1,:)=1;
    % else 
    % p2.tpt_short(:,:)=1;
    % end

    p3.tpt_target_compa0=0.71;
    p3.cfy =(p.house_size-1);
    p3.hhc_a_cov_tpt   =[1,0.5,0.5,0.5];
    p3.hhc_a_cov_screen=[1,0.5,0.5,0.5];
    p3.hiv_cov_tpt=0.9;
    p3.hiv_cov_screen=0.9;

     if(strcmp(p.country,"GEO")||strcmp(p.country,"BRA"))

         hhc_slum=1;
     else

         hhc_slum=[1 2];
     end

    [r3, p3]=get_country_pars(r3,p3,p.country,'TPT_basic','hhc_plhiv');


    %% PLHIV
    % Newly starting ART
    p3.case_findingPLHIV_U(1) =p.hiv_a0_basic_noart_spec;
    p3.case_findingPLHIV_U(2) =p.hiv_a5_basic_noart_spec;
    p3.case_findingPLHIV_U(3) =p.hiv_a10_basic_noart_spec;
    p3.case_findingPLHIV_U(4) =p.hiv_a15_basic_noart_spec;

    p3.case_findingPLHIV_Lf(1) =p.hiv_a0_basic_noart_spec;
    p3.case_findingPLHIV_Lf(2) =p.hiv_a5_basic_noart_spec;
    p3.case_findingPLHIV_Lf(3) =p.hiv_a10_basic_noart_spec;
    p3.case_findingPLHIV_Lf(4) =p.hiv_a15_basic_noart_spec;

    p3.case_findingPLHIV_Ls(1) =p.hiv_a0_basic_noart_spec;
    p3.case_findingPLHIV_Ls(2) =p.hiv_a5_basic_noart_spec;
    p3.case_findingPLHIV_Ls(3) =p.hiv_a10_basic_noart_spec;
    p3.case_findingPLHIV_Ls(4) =p.hiv_a15_basic_noart_spec;

    p3.case_findingPLHIV_I(1) =p.hiv_a0_basic_noart_sens;
    p3.case_findingPLHIV_I(2) =p.hiv_a5_basic_noart_sens;
    p3.case_findingPLHIV_I(3) =p.hiv_a10_basic_noart_sens;
    p3.case_findingPLHIV_I(4) =p.hiv_a15_basic_noart_sens;

    % those on ART yearly screen
    p3.case_finding_U(1,1,3) = p.hiv_a0_basic_art_spec;
    p3.case_finding_U(2,1,3) = p.hiv_a5_basic_art_spec;
    p3.case_finding_U(3,1,3) = p.hiv_a10_basic_art_spec;
    p3.case_finding_U(4,1,3) = p.hiv_a15_basic_art_spec;

    p3.case_finding_Lf(1,1,3) = p.hiv_a0_basic_art_spec;
    p3.case_finding_Lf(2,1,3) = p.hiv_a5_basic_art_spec;
    p3.case_finding_Lf(3,1,3) = p.hiv_a10_basic_art_spec;
    p3.case_finding_Lf(4,1,3) = p.hiv_a15_basic_art_spec;

    p3.case_finding_Ls(1,1,3) = p.hiv_a0_basic_art_spec;
    p3.case_finding_Ls(2,1,3) = p.hiv_a5_basic_art_spec;
    p3.case_finding_Ls(3,1,3) = p.hiv_a10_basic_art_spec;
    p3.case_finding_Ls(4,1,3) = p.hiv_a15_basic_art_spec;

    p3.case_finding_I(1,1,3) = p.hiv_a0_basic_art_sens;
    p3.case_finding_I(2,1,3) = p.hiv_a5_basic_art_sens;
    p3.case_finding_I(3,1,3) = p.hiv_a10_basic_art_sens;
    p3.case_finding_I(4,1,3) = p.hiv_a15_basic_art_sens;

    % HHC
    p3.case_findingHHC_U(1,:,hhc_slum) = p.hhc_a0_basic_spec;
    p3.case_findingHHC_U(2,:,hhc_slum) = p.hhc_a5_basic_spec*(1-p.tst_spec);
    p3.case_findingHHC_U(3,:,hhc_slum) = p.hhc_a10_basic_spec*(1-p.tst_spec);
    p3.case_findingHHC_U(4,:,hhc_slum) = p.hhc_a15_basic_spec*(1-p.tst_spec);

    p3.case_findingHHC_Lf(1,:,hhc_slum) = p.hhc_a0_basic_spec;
    p3.case_findingHHC_Lf(2,:,hhc_slum) = p.hhc_a5_basic_spec*p.tst_sens;
    p3.case_findingHHC_Lf(3,:,hhc_slum) = p.hhc_a10_basic_spec*p.tst_sens;
    p3.case_findingHHC_Lf(4,:,hhc_slum) = p.hhc_a15_basic_spec*p.tst_sens;

    p3.case_findingHHC_Ls(1,:,hhc_slum) = p.hhc_a0_basic_spec;
    p3.case_findingHHC_Ls(2,:,hhc_slum) = p.hhc_a5_basic_spec*p.tst_sens;
    p3.case_findingHHC_Ls(3,:,hhc_slum) = p.hhc_a10_basic_spec*p.tst_sens;
    p3.case_findingHHC_Ls(4,:,hhc_slum) = p.hhc_a15_basic_spec*p.tst_sens;

    p3.case_findingHHC_I(1,:,hhc_slum) = p.hhc_a0_basic_sens;
    p3.case_findingHHC_I(2,:,hhc_slum) = p.hhc_a5_basic_sens;
    p3.case_findingHHC_I(3,:,hhc_slum) = p.hhc_a10_basic_sens;
    p3.case_findingHHC_I(4,:,hhc_slum) = p.hhc_a15_basic_sens;


    Mset = make_model3(p3, r3, i, s, gps);


elseif intvn == 4
    %% PLHIV+HHC+High risk
    %% PLHIV + HHC
    p4=p;
    r4=r;
    p4.sl_short=1;
    p4.tpt_short(:,:)=p.tpt_short_itv;
    p4.tpt_short(:,:)=1;
    p4.tpt_target_compa0=0.71;
    p4.cfy =(p.house_size-1);
    p4.hhc_a_cov_tpt   =[1,0.5,0.5,0.5];
    p4.hhc_a_cov_screen=[1,0.5,0.5,0.5];
    p4.hiv_cov_tpt=0.9;
    p4.hiv_cov_screen=0.9;
    p4.slum_enhance=0;
    p4.slum_cov_tpt=0.6;
    p4.slum_cov_screen=0.6;

    if(strcmp(p.country,"GEO")||strcmp(p.country,"BRA"))

         hhc_slum=1;
     else

         hhc_slum=[1 2];
     end

   
   [r4, p4]=get_country_pars(r4,p4,p.country,'TPT_basic','slum');

 
    %% PLHIV
    % Newly starting ART
    p4.case_findingPLHIV_U(1) =p.hiv_a0_basic_noart_spec;
    p4.case_findingPLHIV_U(2) =p.hiv_a5_basic_noart_spec;
    p4.case_findingPLHIV_U(3) =p.hiv_a10_basic_noart_spec;
    p4.case_findingPLHIV_U(4) =p.hiv_a15_basic_noart_spec;

    p4.case_findingPLHIV_Lf(1) =p.hiv_a0_basic_noart_spec;
    p4.case_findingPLHIV_Lf(2) =p.hiv_a5_basic_noart_spec;
    p4.case_findingPLHIV_Lf(3) =p.hiv_a10_basic_noart_spec;
    p4.case_findingPLHIV_Lf(4) =p.hiv_a15_basic_noart_spec;

    p4.case_findingPLHIV_Ls(1) =p.hiv_a0_basic_noart_spec;
    p4.case_findingPLHIV_Ls(2) =p.hiv_a5_basic_noart_spec;
    p4.case_findingPLHIV_Ls(3) =p.hiv_a10_basic_noart_spec;
    p4.case_findingPLHIV_Ls(4) =p.hiv_a15_basic_noart_spec;

    p4.case_findingPLHIV_I(1) =p.hiv_a0_basic_noart_sens;
    p4.case_findingPLHIV_I(2) =p.hiv_a5_basic_noart_sens;
    p4.case_findingPLHIV_I(3) =p.hiv_a10_basic_noart_sens;
    p4.case_findingPLHIV_I(4) =p.hiv_a15_basic_noart_sens;

    % those on ART yearly screen
    p4.case_finding_U(1,1,3) = p.hiv_a0_basic_art_spec;
    p4.case_finding_U(2,1,3) = p.hiv_a5_basic_art_spec;
    p4.case_finding_U(3,1,3) = p.hiv_a10_basic_art_spec;
    p4.case_finding_U(4,1,3) = p.hiv_a15_basic_art_spec;

    p4.case_finding_Lf(1,1,3) = p.hiv_a0_basic_art_spec;
    p4.case_finding_Lf(2,1,3) = p.hiv_a5_basic_art_spec;
    p4.case_finding_Lf(3,1,3) = p.hiv_a10_basic_art_spec;
    p4.case_finding_Lf(4,1,3) = p.hiv_a15_basic_art_spec;

    p4.case_finding_Ls(1,1,3) = p.hiv_a0_basic_art_spec;
    p4.case_finding_Ls(2,1,3) = p.hiv_a5_basic_art_spec;
    p4.case_finding_Ls(3,1,3) = p.hiv_a10_basic_art_spec;
    p4.case_finding_Ls(4,1,3) = p.hiv_a15_basic_art_spec;

    p4.case_finding_I(1,1,3) = p.hiv_a0_basic_art_sens;
    p4.case_finding_I(2,1,3) = p.hiv_a5_basic_art_sens;
    p4.case_finding_I(3,1,3) = p.hiv_a10_basic_art_sens;
    p4.case_finding_I(4,1,3) = p.hiv_a15_basic_art_sens;

    % HHC
    p4.case_findingHHC_U(1,:,hhc_slum) = p.hhc_a0_basic_spec;
    p4.case_findingHHC_U(2,:,hhc_slum) = p.hhc_a5_basic_spec*(1-p.tst_spec);
    p4.case_findingHHC_U(3,:,hhc_slum) = p.hhc_a10_basic_spec*(1-p.tst_spec);
    p4.case_findingHHC_U(4,:,hhc_slum) = p.hhc_a15_basic_spec*(1-p.tst_spec);

    p4.case_findingHHC_Lf(1,:,hhc_slum) = p.hhc_a0_basic_spec;
    p4.case_findingHHC_Lf(2,:,hhc_slum) = p.hhc_a5_basic_spec*p.tst_sens;
    p4.case_findingHHC_Lf(3,:,hhc_slum) = p.hhc_a10_basic_spec*p.tst_sens;
    p4.case_findingHHC_Lf(4,:,hhc_slum) = p.hhc_a15_basic_spec*p.tst_sens;

    p4.case_findingHHC_Ls(1,:,hhc_slum) = p.hhc_a0_basic_spec;
    p4.case_findingHHC_Ls(2,:,hhc_slum) = p.hhc_a5_basic_spec*p.tst_sens;
    p4.case_findingHHC_Ls(3,:,hhc_slum) = p.hhc_a10_basic_spec*p.tst_sens;
    p4.case_findingHHC_Ls(4,:,hhc_slum) = p.hhc_a15_basic_spec*p.tst_sens;

    p4.case_findingHHC_I(1,:,hhc_slum) = p.hhc_a0_basic_sens;
    p4.case_findingHHC_I(2,:,hhc_slum) = p.hhc_a5_basic_sens;
    p4.case_findingHHC_I(3,:,hhc_slum) = p.hhc_a10_basic_sens;
    p4.case_findingHHC_I(4,:,hhc_slum) = p.hhc_a15_basic_sens;

    %High risk
    p4.case_findingslum_U(4,2,1) = p.slum_a15_basic_spec*(1-p.tst_spec);
    p4.case_findingslum_Lf(4,2,1) = p.slum_a15_basic_spec*p.tst_sens;
    p4.case_findingslum_Ls(4,2,1) = p.slum_a15_basic_spec*p.tst_sens;
    p4.case_findingslum_I(4,2,1) = p.slum_a15_basic_sens;


    Mset = make_model3(p4, r4, i, s, gps);



end
M = Mset;



