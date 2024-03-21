%Baseline run

function M = make_intervention0(p, r, i,s, gps,intvn_seq)



% --- Baseline state
M0 = make_model3(p, r, i, s, gps);

for iseq = 1:length(intvn_seq)
    intvn = intvn_seq(iseq);


    if intvn == 0

        Mset =M0;

    elseif intvn == 1
        %% PLHIV
        % Newly starting ART
        p1=p;
        r1=r;
        p1.case_findingPLHIV_U(1) =(1-p.IPThiv)*p.hiv_a0_basic_noart_spec*0;
        p1.case_findingPLHIV_U(2) =(1-p.IPThiv)*p.hiv_a5_basic_noart_spec*0;
        p1.case_findingPLHIV_U(3) =(1-p.IPThiv)*p.hiv_a10_basic_noart_spec*0;
        p1.case_findingPLHIV_U(4) =(1-p.IPThiv)*p.hiv_a15_basic_noart_spec*0;

        p1.case_findingPLHIV_Lf(1) =(1-p.IPThiv)*p.hiv_a0_basic_noart_spec*0;
        p1.case_findingPLHIV_Lf(2) =(1-p.IPThiv)*p.hiv_a5_basic_noart_spec*0;
        p1.case_findingPLHIV_Lf(3) =(1-p.IPThiv)*p.hiv_a10_basic_noart_spec*0;
        p1.case_findingPLHIV_Lf(4) =(1-p.IPThiv)*p.hiv_a15_basic_noart_spec*0;

        p1.case_findingPLHIV_Ls(1) =(1-p.IPThiv)*p.hiv_a0_basic_noart_spec*0;
        p1.case_findingPLHIV_Ls(2) =(1-p.IPThiv)*p.hiv_a5_basic_noart_spec*0;
        p1.case_findingPLHIV_Ls(3) =(1-p.IPThiv)*p.hiv_a10_basic_noart_spec*0;
        p1.case_findingPLHIV_Ls(4) =(1-p.IPThiv)*p.hiv_a15_basic_noart_spec*0;

        p1.case_findingPLHIV_I(1) =p.hiv_a0_basic_noart_sens*0;
        p1.case_findingPLHIV_I(2) =p.hiv_a5_basic_noart_sens*0;
        p1.case_findingPLHIV_I(3) =p.hiv_a10_basic_noart_sens*0;
        p1.case_findingPLHIV_I(4) =p.hiv_a15_basic_noart_sens*0;

        % those on ART yearly screen
        p1.case_finding_U(1,1,3) = p.hiv_a0_basic_noart_spec*0;
        p1.case_finding_U(1,1,3) = p.hiv_a0_basic_art_spec*0;
        p1.case_finding_U(2,1,3) = p.hiv_a5_basic_noart_spec*0;
        p1.case_finding_U(2,1,3) = p.hiv_a5_basic_art_spec*0;
        p1.case_finding_U(3,1,3) = p.hiv_a10_basic_noart_spec*0;
        p1.case_finding_U(3,1,3) = p.hiv_a10_basic_art_spec*0;
        p1.case_finding_U(4,1,3) = p.hiv_a15_basic_noart_spec*0;
        p1.case_finding_U(4,1,3) = p.hiv_a15_basic_art_spec*0;

        p1.case_finding_Lf(1,1,3) = p.hiv_a0_basic_noart_spec*0;
        p1.case_finding_Lf(1,1,3) = p.hiv_a0_basic_art_spec*0;
        p1.case_finding_Lf(2,1,3) = p.hiv_a5_basic_noart_spec*0;
        p1.case_finding_Lf(2,1,3) = p.hiv_a5_basic_art_spec*0;
        p1.case_finding_Lf(3,1,3) = p.hiv_a10_basic_noart_spec*0;
        p1.case_finding_Lf(3,1,3) = p.hiv_a10_basic_art_spec*0;
        p1.case_finding_Lf(4,1,3) = p.hiv_a15_basic_noart_spec*0;
        p1.case_finding_Lf(4,1,3) = p.hiv_a15_basic_art_spec*0;

        p1.case_finding_Ls(1,1,3) = p.hiv_a0_basic_noart_spec*0;
        p1.case_finding_Ls(1,1,3) = p.hiv_a0_basic_art_spec*0;
        p1.case_finding_Ls(2,1,3) = p.hiv_a5_basic_noart_spec*0;
        p1.case_finding_Ls(2,1,3) = p.hiv_a5_basic_art_spec*0;
        p1.case_finding_Ls(3,1,3) = p.hiv_a10_basic_noart_spec*0;
        p1.case_finding_Ls(3,1,3) = p.hiv_a10_basic_art_spec*0;
        p1.case_finding_Ls(4,1,3) = p.hiv_a15_basic_noart_spec*0;
        p1.case_finding_Ls(4,1,3) = p.hiv_a15_basic_art_spec*0;

        p1.case_finding_I(1,1,3) = p.hiv_a0_basic_noart_sens*0;
        p1.case_finding_I(1,1,3) = p.hiv_a0_basic_art_sens*0;
        p1.case_finding_I(2,1,3) = p.hiv_a5_basic_noart_sens*0;
        p1.case_finding_I(2,1,3) = p.hiv_a5_basic_art_sens*0;
        p1.case_finding_I(3,1,3) = p.hiv_a10_basic_noart_sens*0;
        p1.case_finding_I(3,1,3) = p.hiv_a10_basic_art_sens*0;
        p1.case_finding_I(4,1,3) = p.hiv_a15_basic_noart_sens*0;
        p1.case_finding_I(4,1,3) = p.hiv_a15_basic_art_sens*0;


        Mset = make_model3(p1, r1, i, s, gps);

    elseif intvn == 2
        %% HHC
        % Newly starting ART
        p2=p;
        r2=r;
        p2.cfy =(p.house_size-1)*0;
        p2.hhc_a_cov=[0.9,0.5,0.5].*0;

        p2.case_findingHHC_U(1,1,1) = p.hhc_a0_basic_spec*p.hhc_distr(1)*0;
        p2.case_findingHHC_U(2,1,1) = p.hhc_a5_basic_spec*p.tst_spec*p.hhc_distr(1)*0;
        p2.case_findingHHC_U(3,1,1) = p.hhc_a10_basic_spec*p.tst_spec*p.hhc_distr(1)*0;
        p2.case_findingHHC_U(1,1,1) = p.hhc_a15_basic_spec*p.tst_spec*p.hhc_distr(1)*0;

        p2.case_findingHHC_Lf(1,1,1) = p.hhc_a0_basic_spec*p.hhc_distr(2)*0;
        p2.case_findingHHC_Lf(2,1,1) = p.hhc_a5_basic_spec*p.tst_sens*p.hhc_distr(2)*0;
        p2.case_findingHHC_Lf(3,1,1) = p.hhc_a10_basic_spec*p.tst_sens*p.hhc_distr(2)*0;
        p2.case_findingHHC_Lf(1,1,1) = p.hhc_a15_basic_spec*p.tst_sens*p.hhc_distr(2)*0;

        p2.case_findingHHC_Ls(1,1,1) = p.hhc_a0_basic_spec*p.hhc_distr(3)*0;
        p2.case_findingHHC_Ls(2,1,1) = p.hhc_a5_basic_spec*p.tst_sens*p.hhc_distr(3)*0;
        p2.case_findingHHC_Ls(3,1,1) = p.hhc_a10_basic_spec*p.tst_sens*p.hhc_distr(3)*0;
        p2.case_findingHHC_Ls(1,1,1) = p.hhc_a15_basic_spec*p.tst_sens*p.hhc_distr(3)*0;

        p2.case_findingHHC_I(1,1,1) = p.hhc_a0_basic_sens*p.hhc_distr(4)*0;
        p2.case_findingHHC_I(2,1,1) = p.hhc_a5_basic_sens*p.hhc_distr(4)*0;
        p2.case_findingHHC_I(3,1,1) = p.hhc_a10_basic_sens*p.hhc_distr(4)*0;
        p2.case_findingHHC_I(1,1,1) = p.hhc_a15_basic_sens*p.hhc_distr(4)*0;
        
   
        Mset = make_model3(p2, r2, i, s, gps);


    elseif intvn == 3
        %% PLHIV + HHC
        p3=p;
        r3=r;
        p3.cfy =(p.house_size-1)*0;
        p3.hhc_a_cov=[0.9,0.5,0.5].*0;

              %% PLHIV
        % Newly starting ART
        p3.case_findingPLHIV_U(1) =(1-p.IPThiv)*p.hiv_a0_basic_noart_spec*0;
        p3.case_findingPLHIV_U(2) =(1-p.IPThiv)*p.hiv_a5_basic_noart_spec*0;
        p3.case_findingPLHIV_U(3) =(1-p.IPThiv)*p.hiv_a10_basic_noart_spec*0;
        p3.case_findingPLHIV_U(4) =(1-p.IPThiv)*p.hiv_a15_basic_noart_spec*0;

        p3.case_findingPLHIV_Lf(1) =(1-p.IPThiv)*p.hiv_a0_basic_noart_spec*0;
        p3.case_findingPLHIV_Lf(2) =(1-p.IPThiv)*p.hiv_a5_basic_noart_spec*0;
        p3.case_findingPLHIV_Lf(3) =(1-p.IPThiv)*p.hiv_a10_basic_noart_spec*0;
        p3.case_findingPLHIV_Lf(4) =(1-p.IPThiv)*p.hiv_a15_basic_noart_spec*0;

        p3.case_findingPLHIV_Ls(1) =(1-p.IPThiv)*p.hiv_a0_basic_noart_spec*0;
        p3.case_findingPLHIV_Ls(2) =(1-p.IPThiv)*p.hiv_a5_basic_noart_spec*0;
        p3.case_findingPLHIV_Ls(3) =(1-p.IPThiv)*p.hiv_a10_basic_noart_spec*0;
        p3.case_findingPLHIV_Ls(4) =(1-p.IPThiv)*p.hiv_a15_basic_noart_spec*0;

        p3.case_findingPLHIV_I(1) =p.hiv_a0_basic_noart_sens*0;
        p3.case_findingPLHIV_I(2) =p.hiv_a5_basic_noart_sens*0;
        p3.case_findingPLHIV_I(3) =p.hiv_a10_basic_noart_sens*0;
        p3.case_findingPLHIV_I(4) =p.hiv_a15_basic_noart_sens*0;

        % those on ART yearly screen
        p3.case_finding_U(1,1,3) = p.hiv_a0_basic_noart_spec*0;
        p3.case_finding_U(1,1,3) = p.hiv_a0_basic_art_spec*0;
        p3.case_finding_U(2,1,3) = p.hiv_a5_basic_noart_spec*0;
        p3.case_finding_U(2,1,3) = p.hiv_a5_basic_art_spec*0;
        p3.case_finding_U(3,1,3) = p.hiv_a10_basic_noart_spec*0;
        p3.case_finding_U(3,1,3) = p.hiv_a10_basic_art_spec*0;
        p3.case_finding_U(4,1,3) = p.hiv_a15_basic_noart_spec*0;
        p3.case_finding_U(4,1,3) = p.hiv_a15_basic_art_spec*0;

        p3.case_finding_Lf(1,1,3) = p.hiv_a0_basic_noart_spec*0;
        p3.case_finding_Lf(1,1,3) = p.hiv_a0_basic_art_spec*0;
        p3.case_finding_Lf(2,1,3) = p.hiv_a5_basic_noart_spec*0;
        p3.case_finding_Lf(2,1,3) = p.hiv_a5_basic_art_spec*0;
        p3.case_finding_Lf(3,1,3) = p.hiv_a10_basic_noart_spec*0;
        p3.case_finding_Lf(3,1,3) = p.hiv_a10_basic_art_spec*0;
        p3.case_finding_Lf(4,1,3) = p.hiv_a15_basic_noart_spec*0;
        p3.case_finding_Lf(4,1,3) = p.hiv_a15_basic_art_spec*0;

        p3.case_finding_Ls(1,1,3) = p.hiv_a0_basic_noart_spec*0;
        p3.case_finding_Ls(1,1,3) = p.hiv_a0_basic_art_spec*0;
        p3.case_finding_Ls(2,1,3) = p.hiv_a5_basic_noart_spec*0;
        p3.case_finding_Ls(2,1,3) = p.hiv_a5_basic_art_spec*0;
        p3.case_finding_Ls(3,1,3) = p.hiv_a10_basic_noart_spec*0;
        p3.case_finding_Ls(3,1,3) = p.hiv_a10_basic_art_spec*0;
        p3.case_finding_Ls(4,1,3) = p.hiv_a15_basic_noart_spec*0;
        p3.case_finding_Ls(4,1,3) = p.hiv_a15_basic_art_spec*0;

        p3.case_finding_I(1,1,3) = p.hiv_a0_basic_noart_sens*0;
        p3.case_finding_I(1,1,3) = p.hiv_a0_basic_art_sens*0;
        p3.case_finding_I(2,1,3) = p.hiv_a5_basic_noart_sens*0;
        p3.case_finding_I(2,1,3) = p.hiv_a5_basic_art_sens*0;
        p3.case_finding_I(3,1,3) = p.hiv_a10_basic_noart_sens*0;
        p3.case_finding_I(3,1,3) = p.hiv_a10_basic_art_sens*0;
        p3.case_finding_I(4,1,3) = p.hiv_a15_basic_noart_sens*0;
        p3.case_finding_I(4,1,3) = p.hiv_a15_basic_art_sens*0;

        % HHC
        p3.case_findingHHC_U(1,1,1) = p.hhc_a0_basic_spec*p.hhc_distr(1)*0;
        p3.case_findingHHC_U(2,1,1) = p.hhc_a5_basic_spec*p.tst_spec*p.hhc_distr(1)*0;
        p3.case_findingHHC_U(3,1,1) = p.hhc_a10_basic_spec*p.tst_spec*p.hhc_distr(1)*0;
        p3.case_findingHHC_U(1,1,1) = p.hhc_a15_basic_spec*p.tst_spec*p.hhc_distr(1)*0;

        p3.case_findingHHC_Lf(1,1,1) = p.hhc_a0_basic_spec*p.hhc_distr(2)*0;
        p3.case_findingHHC_Lf(2,1,1) = p.hhc_a5_basic_spec*p.tst_sens*p.hhc_distr(2)*0;
        p3.case_findingHHC_Lf(3,1,1) = p.hhc_a10_basic_spec*p.tst_sens*p.hhc_distr(2)*0;
        p3.case_findingHHC_Lf(1,1,1) = p.hhc_a15_basic_spec*p.tst_sens*p.hhc_distr(2)*0;

        p3.case_findingHHC_Ls(1,1,1) = p.hhc_a0_basic_spec*p.hhc_distr(3)*0;
        p3.case_findingHHC_Ls(2,1,1) = p.hhc_a5_basic_spec*p.tst_sens*p.hhc_distr(3)*0;
        p3.case_findingHHC_Ls(3,1,1) = p.hhc_a10_basic_spec*p.tst_sens*p.hhc_distr(3)*0;
        p3.case_findingHHC_Ls(1,1,1) = p.hhc_a15_basic_spec*p.tst_sens*p.hhc_distr(3)*0;

        p3.case_findingHHC_I(1,1,1) = p.hhc_a0_basic_sens*p.hhc_distr(4)*0;
        p3.case_findingHHC_I(2,1,1) = p.hhc_a5_basic_sens*p.hhc_distr(4)*0;
        p3.case_findingHHC_I(3,1,1) = p.hhc_a10_basic_sens*p.hhc_distr(4)*0;
        p3.case_findingHHC_I(1,1,1) = p.hhc_a15_basic_sens*p.hhc_distr(4)*0;


        Mset = make_model3(p3, r3, i, s, gps);


    elseif intvn == 4
        %% PLHIV+HHC+High risk
                %% PLHIV + HHC
        p4=p;
        r4=r;
        p4.cfy =(p.house_size-1)*0;
        p4.hhc_a_cov=[0.9,0.5,0.5].*0;

              %% PLHIV
        % Newly starting ART
        p4.case_findingPLHIV_U(1) =(1-p.IPThiv)*p.hiv_a0_basic_noart_spec*0;
        p4.case_findingPLHIV_U(2) =(1-p.IPThiv)*p.hiv_a5_basic_noart_spec*0;
        p4.case_findingPLHIV_U(3) =(1-p.IPThiv)*p.hiv_a10_basic_noart_spec*0;
        p4.case_findingPLHIV_U(4) =(1-p.IPThiv)*p.hiv_a15_basic_noart_spec*0;

        p4.case_findingPLHIV_Lf(1) =(1-p.IPThiv)*p.hiv_a0_basic_noart_spec*0;
        p4.case_findingPLHIV_Lf(2) =(1-p.IPThiv)*p.hiv_a5_basic_noart_spec*0;
        p4.case_findingPLHIV_Lf(3) =(1-p.IPThiv)*p.hiv_a10_basic_noart_spec*0;
        p4.case_findingPLHIV_Lf(4) =(1-p.IPThiv)*p.hiv_a15_basic_noart_spec*0;

        p4.case_findingPLHIV_Ls(1) =(1-p.IPThiv)*p.hiv_a0_basic_noart_spec*0;
        p4.case_findingPLHIV_Ls(2) =(1-p.IPThiv)*p.hiv_a5_basic_noart_spec*0;
        p4.case_findingPLHIV_Ls(3) =(1-p.IPThiv)*p.hiv_a10_basic_noart_spec*0;
        p4.case_findingPLHIV_Ls(4) =(1-p.IPThiv)*p.hiv_a15_basic_noart_spec*0;

        p4.case_findingPLHIV_I(1) =p.hiv_a0_basic_noart_sens*0;
        p4.case_findingPLHIV_I(2) =p.hiv_a5_basic_noart_sens*0;
        p4.case_findingPLHIV_I(3) =p.hiv_a10_basic_noart_sens*0;
        p4.case_findingPLHIV_I(4) =p.hiv_a15_basic_noart_sens*0;

        % those on ART yearly screen
        p4.case_finding_U(1,1,3) = p.hiv_a0_basic_noart_spec*0;
        p4.case_finding_U(1,1,3) = p.hiv_a0_basic_art_spec*0;
        p4.case_finding_U(2,1,3) = p.hiv_a5_basic_noart_spec*0;
        p4.case_finding_U(2,1,3) = p.hiv_a5_basic_art_spec*0;
        p4.case_finding_U(3,1,3) = p.hiv_a10_basic_noart_spec*0;
        p4.case_finding_U(3,1,3) = p.hiv_a10_basic_art_spec*0;
        p4.case_finding_U(4,1,3) = p.hiv_a15_basic_noart_spec*0;
        p4.case_finding_U(4,1,3) = p.hiv_a15_basic_art_spec*0;

        p4.case_finding_Lf(1,1,3) = p.hiv_a0_basic_noart_spec*0;
        p4.case_finding_Lf(1,1,3) = p.hiv_a0_basic_art_spec*0;
        p4.case_finding_Lf(2,1,3) = p.hiv_a5_basic_noart_spec*0;
        p4.case_finding_Lf(2,1,3) = p.hiv_a5_basic_art_spec*0;
        p4.case_finding_Lf(3,1,3) = p.hiv_a10_basic_noart_spec*0;
        p4.case_finding_Lf(3,1,3) = p.hiv_a10_basic_art_spec*0;
        p4.case_finding_Lf(4,1,3) = p.hiv_a15_basic_noart_spec*0;
        p4.case_finding_Lf(4,1,3) = p.hiv_a15_basic_art_spec*0;

        p4.case_finding_Ls(1,1,3) = p.hiv_a0_basic_noart_spec*0;
        p4.case_finding_Ls(1,1,3) = p.hiv_a0_basic_art_spec*0;
        p4.case_finding_Ls(2,1,3) = p.hiv_a5_basic_noart_spec*0;
        p4.case_finding_Ls(2,1,3) = p.hiv_a5_basic_art_spec*0;
        p4.case_finding_Ls(3,1,3) = p.hiv_a10_basic_noart_spec*0;
        p4.case_finding_Ls(3,1,3) = p.hiv_a10_basic_art_spec*0;
        p4.case_finding_Ls(4,1,3) = p.hiv_a15_basic_noart_spec*0;
        p4.case_finding_Ls(4,1,3) = p.hiv_a15_basic_art_spec*0;

        p4.case_finding_I(1,1,3) = p.hiv_a0_basic_noart_sens*0;
        p4.case_finding_I(1,1,3) = p.hiv_a0_basic_art_sens*0;
        p4.case_finding_I(2,1,3) = p.hiv_a5_basic_noart_sens*0;
        p4.case_finding_I(2,1,3) = p.hiv_a5_basic_art_sens*0;
        p4.case_finding_I(3,1,3) = p.hiv_a10_basic_noart_sens*0;
        p4.case_finding_I(3,1,3) = p.hiv_a10_basic_art_sens*0;
        p4.case_finding_I(4,1,3) = p.hiv_a15_basic_noart_sens*0;
        p4.case_finding_I(4,1,3) = p.hiv_a15_basic_art_sens*0;

        % HHC
        p4.case_findingHHC_U(1,1,1) = p.hhc_a0_basic_spec*p.hhc_distr(1)*0;
        p4.case_findingHHC_U(2,1,1) = p.hhc_a5_basic_spec*(1-p.tst_spec)*p.hhc_distr(1)*0;
        p4.case_findingHHC_U(3,1,1) = p.hhc_a10_basic_spec*(1-p.tst_spec)*p.hhc_distr(1)*0;
        p4.case_findingHHC_U(1,1,1) = p.hhc_a15_basic_spec*(1-p.tst_spec)*p.hhc_distr(1)*0;

        p4.case_findingHHC_Lf(1,1,1) = p.hhc_a0_basic_spec*p.hhc_distr(2)*0;
        p4.case_findingHHC_Lf(2,1,1) = p.hhc_a5_basic_spec*p.tst_sens*p.hhc_distr(2)*0;
        p4.case_findingHHC_Lf(3,1,1) = p.hhc_a10_basic_spec*p.tst_sens*p.hhc_distr(2)*0;
        p4.case_findingHHC_Lf(1,1,1) = p.hhc_a15_basic_spec*p.tst_sens*p.hhc_distr(2)*0;

        p4.case_findingHHC_Ls(1,1,1) = p.hhc_a0_basic_spec*p.hhc_distr(3)*0;
        p4.case_findingHHC_Ls(2,1,1) = p.hhc_a5_basic_spec*p.tst_sens*p.hhc_distr(3)*0;
        p4.case_findingHHC_Ls(3,1,1) = p.hhc_a10_basic_spec*p.tst_sens*p.hhc_distr(3)*0;
        p4.case_findingHHC_Ls(1,1,1) = p.hhc_a15_basic_spec*p.tst_sens*p.hhc_distr(3)*0;

        p4.case_findingHHC_I(1,1,1) = p.hhc_a0_basic_sens*p.hhc_distr(4)*0;
        p4.case_findingHHC_I(2,1,1) = p.hhc_a5_basic_sens*p.hhc_distr(4)*0;
        p4.case_findingHHC_I(3,1,1) = p.hhc_a10_basic_sens*p.hhc_distr(4)*0;
        p4.case_findingHHC_I(1,1,1) = p.hhc_a15_basic_sens*p.hhc_distr(4)*0;

        %High risk
        p4.case_findingslum_U(4,2,1) = p.slum_a15_basic_spec*0.6*0;
        p4.case_findingslum_Lf(4,2,1) = p.slum_a15_basic_spec*0.6*0;
        p4.case_findingslum_Ls(4,2,1) = p.slum_a15_basic_spec*0.6*0;
        p4.case_findingslum_I(4,2,1) = p.slum_a15_basic_sens*0.6*0;


        Mset = make_model3(p4, r4, i, s, gps);


    end
end
M = Mset;



