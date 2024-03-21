%This function constructs the Models , rund simulations and give objective indicators back

function [out, aux] = get_objective(data,theta, prm, ref, sel, agg, gps, calfn, hivpoints,endyr)

if(isempty(endyr))
    endyr=2023;
end

X=theta;

r = prm.r; p = prm.p; i = ref.i; s = ref.s; xi = ref.xi;

[r,p] = allocate_parameters(X,r,p,xi);

% Check if parameters are within bounds
tmp = [prm.bds(:,1:length(X)); X]; tmp = diff(tmp([1,3,2],:));

if min(min(tmp)) < 0 %|| p.Dx(2) > p.Dx(1) || r.default(1) > r.default(2)
    out.llk = -Inf; aux = nan;
else
    
    
    % --- Set up the necessary models -----------------------------------------
    % COVID recovery
    p5 = p; r5=r;
    M5 = make_model(p5, r5, i, s, gps);

    % COVID disruption
    p4 = p; r4=r;
    p4.Dx = p.Dx * (1-p.covid_disrupt);
    p4.xpert = p.xpert* (1-p.covid_disrupt);
    M4 = make_model(p4, r4, i, s, gps);


    % Final conditions with full ART and IPT catch up
    p3 = p; r3=r;
    p3.IPThiv = 0;
    M3 = make_model(p3, r3, i, s, gps);
    
    % ART expansion
    p2 = p; r2=r;
    p2.IPThiv = 0;
    r2.ARTrec = r.ARTrec*p.art_cov;
    p2.xpert = p.xpert.*p.ntpcov_dr;
    M2 = make_model(p2, r2, i, s, gps);
    
    % NTP expansion
    p1 = p; r1=r;
    % p1.hivtest_notb = 0;
    % p1.hivtest_tb = 0;
    % p1.art_notb = 0;
    % p1.art_tb = 0;
    p1.IPThiv = 0;
    r1.ARTrec = 0;
    p1.Dx = p.Dx.*p.ntpcov_ds;
    p1.xpert = 0;
    M1 = make_model(p1, r1, i, s, gps);
    
    
    % Equil
    p0 = p; r0=r;
    %r0.beta_mdr = 0;
    r0.MDR_acqu = 0;
    p0.hivtest_notb = 0;
    p0.hivtest_tb = 0;
    p0.art_notb = 0;
    p0.art_tb = 0;
    p0.IPThiv = 0;
    r0.ARTrec = 0;
    p0.Dx = p.Dx.*0;
    p0.xpert = 0;
    M0 = make_model(p0, r0, i, s, gps);
    
    % Equilibrium model
    init = zeros(1,i.nx); seed = 1e-6;
    init(i.U.a0_4.slum.hneg) = (1-seed)*p.frac_pop(1)*p.frac_slum;
    init(i.U.a0_4.no_slum.hneg) = (1-seed)*p.frac_pop(1)*(1-p.frac_slum);
    init(i.U.a5_9.slum.hneg) = (1-seed)*p.frac_pop(2)*p.frac_slum;
    init(i.U.a5_9.no_slum.hneg) = (1-seed)*p.frac_pop(2)*(1-p.frac_slum);
    init(i.U.a10_15.slum.hneg) = (1-seed)*p.frac_pop(3)*p.frac_slum;
    init(i.U.a10_15.no_slum.hneg) = (1-seed)*p.frac_pop(3)*(1-p.frac_slum);
    init(i.U.a15p.slum.hneg) = (1-seed)*p.frac_pop(4)*p.frac_slum;
    init(i.U.a15p.no_slum.hneg) = (1-seed)*p.frac_pop(4)*(1-p.frac_slum);
    init(i.I.a15p.slum.hneg.ds) = seed*0.7*0.7;
    init(i.I.a15p.no_slum.hneg.ds) = seed*0.3*0.7;
    init(i.I.a15p.slum.hneg.mdr) = seed*0.7*0.3;
    init(i.I.a15p.no_slum.hneg.mdr) = seed*0.3*0.3;
  
    odeopts = odeset('NonNegative',[1:i.nstates],'Refine',64,'AbsTol',1e-10,'RelTol',1e-10);


    geq = @(t,in)governing_equations(t, in, M0, i, s, r0, p0, sel, agg,0,hivpoints);
    [t0, soln0] = ode15s(geq, [p.yrntp-600 p.yrntp], init, sodeset('NonNegative',[s.nstates]));

    % NTP
    init = soln0(end,:);
    geq = @(t,in)governing_equations(t, in, M1, i, s,r1, p1, sel, agg,p.growth,hivpoints);
    [t1, soln1] = ode15s(geq, [p.yrntp p.yrart], init,odeset('NonNegative',[s.nstates]));
    
    % ART
    init = soln1(end,:);
    [t2, soln2] = ode15s(@(t,in) goveqs_scaleup(t, in, M1, M2,[p.yrart+1 2007],...
        i, s,r2, p2, sel, agg,hivpoints), [p.yrart 2010] , init, odeset('NonNegative',[s.nstates]));
                                         
    % IPT to newly enrolles
    init = soln2(end,:);
    [t3, soln3] = ode15s(@(t,in) goveqs_scaleup(t, in, M2, M3,[p.art_year p.art_year+3],...
        i, s,r3, p3, sel, agg,hivpoints), [2010 2019] , init,odeset('NonNegative',[s.nstates]));

    % tic;
    % % Covid
    init = soln3(end,:);
    [t4, soln4] = ode15s(@(t,in) goveqs_scaleup(t, in, M3, M4,[2019 2020],...
        i, s,r4, p4, sel, agg,hivpoints), [2019 2021] , init, odeset('NonNegative',[s.nstates]));
    % geq = @(t,in)governing_equations(t, in, M4, i, s, r4, p4, sel, agg,0,hivpoints);
    % [t4, soln4] = ode15s(geq, [2020 2021], init, odeset('NonNegative',[s.nstates]));
  
    % Covid recovery
    init = soln4(end,:);
    [t5, soln5] = ode15s(@(t,in) goveqs_scaleup(t, in, M4, M5,[2021 2023],...
        i, s,r5, p5, sel, agg,hivpoints), [2021 endyr] , init, odeset('NonNegative',[s.nstates]));

    %checks
    if (sum(any(isnan(soln5)))>0)
        out.llk = -Inf; aux = nan;
    else
        
        soln = cat(1,soln3,soln4(2:end,:),soln5(2:end,:));
        t= cat(1,t3,t4(2:end),t5(2:end));
        
        
        soln_long = cat(1,soln0,soln1(2:end,:),soln2(2:end,:),soln3(2:end,:),...
            soln4(2:end,:), soln5(2:end,:));
        t_long =cat(1,t0,t1(2:end),t2(2:end),t3(2:end),t4(2:end),t5(2:end));
        
        
        % --- Get the objectives --------------------------------------------------
        sfin2021 = soln4(end,:); 
        gs = @(cols)(sum(sfin2021(cols)));
        N=gs(1:i.nstates);



        %id=find([t(1):t(end)]==2012)-1;
        yrs=2012:2021;
        %pop vector
        tmp = interp1(t,sum(soln(:,1:i.nstates),2),t(1):t(end));
        pop = (tmp(1:end-1)+tmp(2:end))/2;

        tmp = interp1(t,sum(soln(:,s.a15p),2),t(1):t(end));
        pop_ad = (tmp(1:end-1)+tmp(2:end))/2;
        
        %Inc (All , MDR, TB/HIV)
        tmp = diff(interp1(t,soln(:,i.aux.inc),t(1):t(end)),1);
        inc=1e5*(tmp./[pop' pop' pop'])';
        model.inc_all.est= inc(1,:);
        model.inc_all.year= yrs;
        model.inc_mdr.est= inc(2,:);
        model.inc_mdr.year= yrs;
        model.inc_tbhiv.est= inc(3,:);
        model.inc_tbhiv.year= yrs;
       
        
        %  Mort TB in HIV negative    Mort TB in HIV +
        tmp  = diff(interp1(t,soln(:,i.aux.mort),t(1):t(end)),1);
        mort= 1e5*[tmp(:,2)-tmp(:,3)./pop' ,...
              tmp(:,3)./pop']';

        model.mort_tbhn.est= mort(1,:);
        model.mort_tbhn.year= yrs;
        model.mort_tbhiv.est= mort(2,:);
        model.mort_tbhiv.year= yrs;
             
        
        %Notif (All)
        tmp  = diff(interp1(t,soln(:,i.aux.notif),t(1):t(end)),1);
        notif=1e5*(tmp(:,1)./pop')';
        model.notif.est= notif;
        model.notif.year= yrs;
        

        %Prevalence
        tmp = interp1(t,sum(soln(:,intersect(s.prevalent,s.a15p)),2),t(1):t(end));
        num = (tmp(1:end-1)+tmp(2:end))/2;
        model.prev.est= 1e5*(num./pop_ad);
        model.prev.year= yrs;

        % 
        % 
        % num= gs(intersect(s.prevalent,s.a15p));
        % denom= gs(s.a15p);
        % tbprevalence =1e5*num/denom;
        % 
        % 
        % 
        %HIV cascade
        % HIV Prevalence
        prev_hivad=(gs(s.hivpositive)/N)*1e5;
        tmp  = diff(interp1(t4,soln4(:,i.aux.pt),t4(1):t4(end)),1);
        pthiv = tmp(:,2);
        
        tmp  = diff(interp1(t4,soln4(:,i.aux.newart(1)),t4(1):t4(end)),1);
        
        newart = tmp';
        tmp=pthiv./newart;
        
        pr_onipt=tmp(end);
        idtmp=intersect(s.hivpositive,s.prevalent);
        pr_hivtb=(gs(idtmp)/gs(s.prevalent));
        pr_onart=gs([s.hart])/gs(s.hivpositive);
                    % pr_onart       pr_onipt
        model.pr_onart.est=pr_onart; 
        model.pr_onart.year=t4(end);
        model.pr_onipt.est=pr_onipt; 
        model.pr_onipt.year=t4(end);
        
        
        
        
        
        % Call LLk function
        out.llk = sum(calfn(data,model));
      
        % --- Get additional outputs and package ------------------------------
        
        out.solnL  =[soln_long, t_long];
        out.soln   =[soln, t];
        out.model = model;
     
        
        
    end
end
