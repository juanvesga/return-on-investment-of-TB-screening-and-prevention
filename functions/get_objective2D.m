%This function constructs the Models , rund simulations and give objective indicators back

function [out, aux] = get_objective2D(data,theta, prm, ref, sel, agg, gps, calfn, hivpoints,endyr)

if(isempty(endyr))
    endyr=2019;
end

r = prm.r; p = prm.p; i = ref.i; s = ref.s; xi = ref.xi;

X=theta;

[r,p] = allocate_parameters(X,r,p,xi);

% Check if parameters are within bounds
tmp = [prm.bds(:,1:length(X)); X]; tmp = diff(tmp([1,3,2],:));

if min(min(tmp)) < 0 %|| p.Dx(2) > p.Dx(1) || r.default(1) > r.default(2)
    out.llk = -Inf; aux = nan;
else

    % --- Set up the necessary models -----------------------------------------
    slum_birth=p.geo+p.bra;
    % Equilibrium model
    init = zeros(1,i.nx); seed = 1e-6;

    if (p.geo+p.bra>0)
        init(i.U.a0_4.slum.hneg) = 0;
        init(i.U.a0_4.no_slum.hneg) = (1-seed)*p.frac_pop(1);
        init(i.U.a5_9.slum.hneg) = 0;
        init(i.U.a5_9.no_slum.hneg) = (1-seed)*p.frac_pop(2);
        init(i.U.a10_15.slum.hneg) = 0;
        init(i.U.a10_15.no_slum.hneg) = (1-seed)*p.frac_pop(3);
        init(i.U.a15p.slum.hneg) = (1-seed)*p.frac_pop(4)*p.frac_slum;
        init(i.U.a15p.no_slum.hneg) = (1-seed)*p.frac_pop(4)*(1-p.frac_slum);
        init(i.I.a15p.slum.hneg.ds) = seed*0.3;
        init(i.I.a15p.no_slum.hneg.ds) = seed*0.3;
        init(i.I.a15p.slum.hneg.mdr) = seed*0.2;
        init(i.I.a15p.no_slum.hneg.mdr) = seed*0.2;

    else
        init(i.U.a0_4.slum.hneg) = (1-seed)*p.frac_pop_slum(1)*p.frac_slum*(1-slum_birth);
        init(i.U.a0_4.no_slum.hneg) = (1-seed)*p.frac_pop(1)*(1-p.frac_slum*(1-slum_birth));
        init(i.U.a5_9.slum.hneg) = (1-seed)*p.frac_pop_slum(2)*p.frac_slum*(1-slum_birth);
        init(i.U.a5_9.no_slum.hneg) = (1-seed)*p.frac_pop(2)*(1-p.frac_slum*(1-slum_birth));
        init(i.U.a10_15.slum.hneg) = (1-seed)*p.frac_pop_slum(3)*p.frac_slum*(1-slum_birth);
        init(i.U.a10_15.no_slum.hneg) = (1-seed)*p.frac_pop(3)*(1-p.frac_slum*(1-slum_birth));
        init(i.U.a15p.slum.hneg) = (1-seed)*p.frac_pop_slum(4)*p.frac_slum;
        init(i.U.a15p.no_slum.hneg) = (1-seed)*p.frac_pop(4)*(1-p.frac_slum);
        init(i.I.a15p.slum.hneg.ds) = seed*0.3;
        init(i.I.a15p.no_slum.hneg.ds) = seed*0.3;
        init(i.I.a15p.slum.hneg.mdr) = seed*0.2;
        init(i.I.a15p.no_slum.hneg.mdr) = seed*0.2;


    end

    %% Equilibrium model
    p0 = p; r0=r;
    r0.MDR_acqu=0;
    r0.ART_init = 0;
    p0.IPThiv = 0;
    p0.sl_short =0;
    p0.tpt_short =p.tpt_short*0;
    p0.xpert = 0;
    p0.Tx_init = p0.Tx_init*0;
    p0.Tx_init2 = p0.Tx_init2*0;
    M0 = make_model3(p0, r0, i, s, gps);

    options = odeset('NonNegative',1:i.nstates,'RelTol',1e-7,'AbsTol',1e-7);


    geq = @(t,in)governing_equations_itv2(t, in, M0, i, s, r0, p0, sel, agg,0,hivpoints);
    [t0, soln0] = ode15s(geq, [p.yrntp-1500 p.yrntp], init, options);
    init = soln0(end,:);

    %% Public sector, no ART 1970 ~ 2000
    p1 = p; r1 = r;
    r1.MDR_acqu=r.MDR_acqu*5;
    r1.ART_init = 0;
    p1.IPThiv = 0;
    p1.sl_short =0;
    p1.tpt_short =p.tpt_short*0;
    p1.xpert = 0;


    p1.Tx_init = p.Tx_init;
    p1.Tx_init2 = p.Tx_init2;

    M1 = make_model3(p1, r1, i, s, gps);

    [t1, soln1] = ode15s(@(t,in) goveqs_scaleup_itv(t, in, M0, M1,[p.yrntp p.ART_start],...
        i, s,r1, p1, sel, agg,hivpoints), [p.yrntp p.ART_start] , init, options);

    init = soln1(end,:);


    %% ART ~2000 2010
    p2 = p; r2 = r;

    p2.IPThiv = 0;
    p2.sl_short =0;
    p2.tpt_short =p.tpt_short*0;
    p2.xpert = p.xpert*p.ntpcov_dr;
    r2.ART_init = r.ART_init * 0.1;
    p2.Tx_init = p.Tx_init*p.ntpcov;
    p2.Tx_init2 = p.Tx_init2*p.ntpcov;



    M2 = make_model3(p2, r2, i, s, gps);



    [t2, soln2] = ode15s(@(t,in) goveqs_scaleup_itv(t, in, M1, M2,[p.ART_start p.ART_start+1],...
        i, s,r2, p2, sel, agg,hivpoints), [p.ART_start 2010] , init, options);

    init = soln2(end,:);



    %% Xpert expansion & IPT 2010-2022
    hhc_n=get_screening_pop(soln2(end,:),s);

    p3 = p; r3 = r;

    p3.cfy =(p.house_size-1);


    r3.cseek_fac(1:2)=0.35;%


    tmp=diff(soln2(:,i.aux.lam_ds));
    p3.lam_ds=tmp(end,:);

    tmp=diff(soln2(:,i.aux.lam_dr));
    p3.lam_dr=tmp(end,:);


    p3.U=hhc_n.U;
    p3.Lfds=hhc_n.Lfds;
    p3.Lsds=hhc_n.Lsds;
    p3.Ids=hhc_n.Ids;
    p3.Lfdr=hhc_n.Lfdr;
    p3.Lsdr=hhc_n.Lsdr;
    p3.Idr=hhc_n.Idr;


    p3.hhc_a_cov_tpt   =[p.hhcovu5,p.hhcov,p.hhcov,p.hhcov];
    p3.hhc_a_cov_screen=[p.hhcovu5,p.hhcov,p.hhcov,p.hhcov];

    % Newly starting ART
    p3.case_findingPLHIV_U(1) =p.hiv_a0_baseline_noart_spec;
    p3.case_findingPLHIV_U(2) =p.hiv_a5_baseline_noart_spec;
    p3.case_findingPLHIV_U(3) =p.hiv_a10_baseline_noart_spec;
    p3.case_findingPLHIV_U(4) =p.hiv_a15_baseline_noart_spec;

    p3.case_findingPLHIV_Lf(1) =p.hiv_a0_baseline_noart_spec;
    p3.case_findingPLHIV_Lf(2) =p.hiv_a5_baseline_noart_spec;
    p3.case_findingPLHIV_Lf(3) =p.hiv_a10_baseline_noart_spec;
    p3.case_findingPLHIV_Lf(4) =p.hiv_a15_baseline_noart_spec;

    p3.case_findingPLHIV_Ls(1) =p.hiv_a0_baseline_noart_spec;
    p3.case_findingPLHIV_Ls(2) =p.hiv_a5_baseline_noart_spec;
    p3.case_findingPLHIV_Ls(3) =p.hiv_a10_baseline_noart_spec;
    p3.case_findingPLHIV_Ls(4) =p.hiv_a15_baseline_noart_spec;

    % HHC
    p3.case_findingHHC_U(1,:,1) = p.hhc_a0_baseline_spec;
    p3.case_findingHHC_U(2,:,1) = p.hhc_a5_baseline_spec*(1-p.tst_spec);
    p3.case_findingHHC_U(3,:,1) = p.hhc_a10_baseline_spec*(1-p.tst_spec);
    p3.case_findingHHC_U(4,:,1) = p.hhc_a15_baseline_spec*(1-p.tst_spec);

    p3.case_findingHHC_Lf(1,:,1) = p.hhc_a0_baseline_spec;
    p3.case_findingHHC_Lf(2,:,1) = p.hhc_a5_baseline_spec*p.tst_sens;
    p3.case_findingHHC_Lf(3,:,1) = p.hhc_a10_baseline_spec*p.tst_sens;
    p3.case_findingHHC_Lf(4,:,1) = p.hhc_a15_baseline_spec*p.tst_sens;

    p3.case_findingHHC_Ls(1,:,1) = p.hhc_a0_baseline_spec;
    p3.case_findingHHC_Ls(2,:,1) = p.hhc_a5_baseline_spec*p.tst_sens;
    p3.case_findingHHC_Ls(3,:,1) = p.hhc_a10_baseline_spec*p.tst_sens;
    p3.case_findingHHC_Ls(4,:,1) = p.hhc_a15_baseline_spec*p.tst_sens;

    p3.case_findingHHC_I(1,:,1) = p.hhc_a0_baseline_sens;
    p3.case_findingHHC_I(2,:,1) = p.hhc_a5_baseline_sens;
    p3.case_findingHHC_I(3,:,1) = p.hhc_a10_baseline_sens;
    p3.case_findingHHC_I(4,:,1) = p.hhc_a15_baseline_sens;


    M3 = make_model3(p3, r3, i, s, gps);





    [t3, soln3] = ode15s(@(t,in) goveqs_scaleup_itv(t, in, M2, M3,[2012 2015],...
        i, s,r3, p3, sel, agg,hivpoints), [2010 endyr] , init, options);



    %% checks and Model output
    if (sum(any(isnan(soln2)))>0)
        out.llk = -Inf; aux = nan;
    else

        soln = soln3;
        t= t3;

        soln_long = cat(1,soln0,soln1(2:end,:),soln2(2:end,:),soln3(2:end,:));
        t_long =cat(1,t0,t1(2:end),t2(2:end),t3(2:end));


        % --- Get the objectives --------------------------------------------------
        sfin = soln(end,:);
        gs = @(cols)(sum(sfin(cols)));

        yrs=2012:endyr;
        %pop vector
        tmp = interp1(t,sum(soln(:,1:i.nstates),2),t(1):t(end));
        pop = (tmp(1:end-1)+tmp(2:end))/2;

        tmp = interp1(t,sum(soln(:,s.a15p),2),t(1):t(end));
        pop_ad = (tmp(1:end-1)+tmp(2:end))/2;



        tmp = interp1(t,sum(soln(:,s.slum),2),t(1):t(end));
        pop_slum = (tmp(1:end-1)+tmp(2:end))/2;

        tmp = interp1(t,sum(soln(:,intersect(s.slum,s.a15p)),2),t(1):t(end));
        pop_slum_a15 = (tmp(1:end-1)+tmp(2:end))/2;


        tmp = interp1(t,sum(soln(:,intersect(s.slum,s.ch )),2),t(1):t(end));
        pop_slum_a0_15 = (tmp(1:end-1)+tmp(2:end))/2;


        tmp = interp1(t,sum(soln(:,s.hivpositive),2),t(1):t(end));
        pop_plhiv = ((tmp(1:end-1)+tmp(2:end))/2).*p.pop;

        % plot(pop_plhiv)





        %Inc (All , MDR, TB/HIV)
        tmp = diff(interp1(t,soln(:,i.aux.inc),t(1):t(end)),1);
        inc=1e5*(tmp./[pop' pop' pop' pop_slum'])';
        model.inc_all.est= inc(1,:);
        model.inc_all.year= yrs;
        model.inc_mdr.est= inc(2,:);
        model.inc_mdr.year= yrs;
        model.inc_tbhiv.est= inc(3,:);
        model.inc_tbhiv.year= yrs;
        model.inc_slum.est= inc(4,:);
        model.inc_slum.year= yrs;


        %  Mort TB in HIV negative    Mort TB in HIV +
        tmp  = diff(interp1(t,soln(:,i.aux.mort),t(1):t(end)),1);
        mort= 1e5*[tmp(:,1)./pop' ,...
            tmp(:,2)./pop']';

        model.mort_tbhn.est= mort(1,:);
        model.mort_tbhn.year= yrs;
        model.mort_tbhiv.est= mort(2,:);
        model.mort_tbhiv.year= yrs;


        %Notif (All)
        tmp  = diff(interp1(t,soln(:,i.aux.notif),t(1):t(end)),1);
        notif=1e5*(tmp(:,1)./pop')';
        notif_cases=tmp(:,1).*p.pop;
        model.notif_cases.est= notif_cases(end);
        model.notif.est= notif;
        model.notif.year= yrs;


        %Prevalence
        tmp = interp1(t,sum(soln(:,s.prevalent),2),t(1):t(end));
        num = (tmp(1:end-1)+tmp(2:end))/2;



        model.prev.est= 1e5*(num./pop_ad);
        model.prev.year= yrs;



        tmp = interp1(t,sum(soln(:,intersect(s.prevalent,intersect(s.slum,s.ch))),2),t(1):t(end));
        num = (tmp(1:end-1)+tmp(2:end))/2;
        model.prev_hi0.est= 1e5*(num./pop_slum_a0_15);
        model.prev_hi0.year= yrs;


        if (endyr >2022)
            tend=  2022;
        else
            tend=t(end);
        end

        tmp = interp1(t,sum(soln(:,intersect(s.prevalent,intersect(s.slum,s.ad))),2),t(1):t(end));
        num = (tmp(1:end-1)+tmp(2:end))/2;
        model.prev_hi.est= 1e5*(num./pop_slum_a15);
        model.prev_hi.year= yrs;




        % Treatment coverage

        model.txcov.est= 1e5*notif_cases(end)/(model.inc_all.est(end)/1e5*p.pop);
        model.txcov.year=tend;


        %% Screening
        tmp  = diff(interp1(t,soln(:,i.aux.screen_hhc),t(1):t(end)),1);
        screenhhc=sum(tmp*p.pop,2);



        %% HHC on TPT




        tmp  = diff(interp1(t,soln(:,i.aux.tpt_hhc),t(1):t(end)),1);
        tpthhc=tmp(end-3,:)*p.pop;

        model.tpt_hhc_u5.est=tpthhc(1);
        model.tpt_hhc_u5.year= tend-3;
        model.tpt_hhc.est=sum(tpthhc);
        model.tpt_hhc.year= tend-3;


        %HIV cascade

        % HIV Prevalence
        tmp  = diff(interp1(t,soln(:,i.aux.pt),t(1):t(end)),1);
        pthiv = tmp(:,2);

        tmp  = diff(interp1(t,soln(:,i.aux.newart(1)),t(1):t(end)),1);

        newart = tmp';
        tmp=pthiv./newart;

        pr_onipt=1e5*tmp(end);
        pr_onart=1e5*gs([s.hart])/gs(s.hivpositive);


        model.pr_onart.est=pr_onart;
        model.pr_onart.year= tend-3;
        model.pr_onipt.est=pr_onipt;
        model.pr_onipt.year=tend-3;

        % Call LLk function
        llks=calfn(model);
        out.llk = sum(llks);

        % --- Get additional outputs and package ------------------------------

        out.solnL  =[soln_long, t_long];
        out.soln   =[soln, t];
        out.sfin  = soln3(end,:);
        out.model = model;
        %out.hhc_fac=hhc_fac;




    end
end
