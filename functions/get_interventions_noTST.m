function results=get_interventions_noTST(x,sfin,path,prm,ref,sel,agg,gps,location,scen,hivpoints,fityr,endyr)

rsto = prm.r; psto = prm.p;
i = ref.i; s = ref.s;
xi = ref.xi;

% handles
if (scen==1) % Basic
    fn_mk = @(p,r,intvn_seq)  make_intervention1_noTST(p, r, i, s, gps,intvn_seq);

elseif (scen==2) % Enhanced
    fn_mk = @(p,r,intvn_seq)  make_intervention2_noTST(p, r, i, s, gps,intvn_seq);
end


xvals=x;
nsims=size(x,1);
runtype='all';


itv_labs={'PLHIV','HHC','PLHIV + HHC','PLHIV + HHC + High risk'};

[ inc_ds, inc_dr, mu_ds, mu_dr]=...
    deal(zeros(nsims,(2050-fityr),3*numel(gps.age)+1, length(path)+1));

[screen_all, tst_all, tpt_all, tx_fl_all,tx_sl_all,...
    screen_plhiv, tst_plhiv, tpt_plhiv, tx_fl_plhiv,tx_sl_plhiv,...
    screen_hhc, tst_hhc, tpt_hhc, tx_fl_hhc,tx_sl_hhc,...
    screen_slum, tst_slum, tpt_slum, tx_fl_slum,tx_sl_slum,...
    daly_all, daly_plhiv, daly_slum, yll_all, yll_plhiv, yll_slum,yld_all, yld_plhiv, yld_slum]=...
    deal(zeros(nsims,(2050-fityr),numel(gps.age)+1, length(path)+1));

[prevtb, prevtb_hpos, prevtb_hart, prevtb_slum, prevltbi, prevltbi_hpos,...
    prevltbi_hart, prevltbi_slum]=...
    deal(zeros(nsims,(2050-fityr),numel(gps.age)+1, length(path)+1));

[alive_all, alive_hpos, alive_hart, alive_slum, alive_hhc]=...
    deal(zeros(nsims,(2050-fityr),numel(gps.age)+1, length(path)+1));


[avagemort,avagemort_plhiv, avagemort_slum, inc_all, inc_mdr,...
    mort_tb]=...
    deal(zeros(nsims,(2050-fityr),length(path)+1));

sfins                    =zeros(nsims,i.nx,length(path)+1);



%
parpool('local',8);
ppm = ParforProgMon(location, nsims);
parfor ii=1:nsims
    ppm.increment();

 % for ii=  1:nsims
 % 
 % disp(ii);
% profile on

    x0=xvals(ii,:);
    [r,p] = allocate_parameters(x0,rsto,psto,xi);
    init=sfin(ii,:);


    % Baseline
    Mset = cell(1,length(path)+1);


    for it=1:length(path)+1

        Mset{it}=fn_mk(p,r,it-1);%make_intervention(p, r, i, s, gps,intvn_seq);

    end
    M0=Mset{1};

    options = odeset('NonNegative',1:i.nstates,'RelTol',1e-7,'AbsTol',1e-7);

    getsol = @(Mfinal) ode45(@(t,in)...
        goveqs_scaleup_itv_noTST(t, in, M0, Mfinal,[2023 2030], i,s,r,p, sel, agg,hivpoints),...
        (fityr:endyr) , init, options);


    allsol=zeros(numel(fityr:endyr), i.nx, length(path)+1);
    for mi = 1:length(Mset)

        [t0, tmp] = getsol(Mset{mi});

        allsol(:,:,mi) = tmp;
    end


    %% Denominators for prevalence
    tmp = squeeze(sum(allsol(:,1:i.nstates,:),2));
    pop_all = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,s.a0_4,:),2));
    pop_a0 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,s.a5_9,:),2));
    pop_a5 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,s.a10_15,:),2));
    pop_a10 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,s.a15p,:),2));
    pop_a15p = (tmp(1:end-1,:)+tmp(2:end,:))/2;

    tmp = squeeze(sum(allsol(:,s.hpos,:),2));
    pop_hpos = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.hpos,s.a0_4),:),2));
    pop_hpos_a0 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.hpos,s.a5_9),:),2));
    pop_hpos_a5 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.hpos,s.a10_15),:),2));
    pop_hpos_a10 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.hpos,s.a15p),:),2));
    pop_hpos_a15p = (tmp(1:end-1,:)+tmp(2:end,:))/2;

    tmp = squeeze(sum(allsol(:,s.hart,:),2));
    pop_hart = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.hart,s.a0_4),:),2));
    pop_hart_a0 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.hart,s.a5_9),:),2));
    pop_hart_a5 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.hart,s.a10_15),:),2));
    pop_hart_a10 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.hart,s.a15p),:),2));
    pop_hart_a15p = (tmp(1:end-1,:)+tmp(2:end,:))/2;

    tmp = squeeze(sum(allsol(:,s.slum,:),2));
    pop_slum = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.slum,s.a0_4),:),2));
    pop_slum_a0 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.slum,s.a5_9),:),2));
    pop_slum_a5 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.slum,s.a10_15),:),2));
    pop_slum_a10 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.slum,s.a15p),:),2));
    pop_slum_a15p = (tmp(1:end-1,:)+tmp(2:end,:))/2;

    %% numerators for tb prevalence
    tmp = squeeze(sum(allsol(:,s.prevalent,:),2));
    ntb_all = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.prevalent,s.a0_4),:),2));
    ntb_a0 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.prevalent,s.a5_9),:),2));
    ntb_a5 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.prevalent,s.a10_15),:),2));
    ntb_a10 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.prevalent,s.a15p),:),2));
    ntb_a15p = (tmp(1:end-1,:)+tmp(2:end,:))/2;

    tmp = squeeze(sum(allsol(:,intersect(s.prevalent,s.hpos),:),2));
    ntb_hpos = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.prevalent,intersect(s.hpos,s.a0_4)),:),2));
    ntb_hpos_a0 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.prevalent,intersect(s.hpos,s.a5_9)),:),2));
    ntb_hpos_a5 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.prevalent,intersect(s.hpos,s.a10_15)),:),2));
    ntb_hpos_a10 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.prevalent,intersect(s.hpos,s.a15p)),:),2));
    ntb_hpos_a15p = (tmp(1:end-1,:)+tmp(2:end,:))/2;

    tmp = squeeze(sum(allsol(:,intersect(s.prevalent,s.hart),:),2));
    ntb_hart = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.prevalent,intersect(s.hart,s.a0_4)),:),2));
    ntb_hart_a0 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.prevalent,intersect(s.hart,s.a5_9)),:),2));
    ntb_hart_a5 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.prevalent,intersect(s.hart,s.a10_15)),:),2));
    ntb_hart_a10 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.prevalent,intersect(s.hart,s.a15p)),:),2));
    ntb_hart_a15p = (tmp(1:end-1,:)+tmp(2:end,:))/2;

    tmp = squeeze(sum(allsol(:,intersect(s.prevalent,s.slum),:),2));
    ntb_slum = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.prevalent,intersect(s.slum,s.a0_4)),:),2));
    ntb_slum_a0 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.prevalent,intersect(s.slum,s.a5_9)),:),2));
    ntb_slum_a5 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.prevalent,intersect(s.slum,s.a10_15)),:),2));
    ntb_slum_a10 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.prevalent,intersect(s.slum,s.a15p)),:),2));
    ntb_slum_a15p = (tmp(1:end-1,:)+tmp(2:end,:))/2;

    %% numerators for ltbi prevalence
    tmp = squeeze(sum(allsol(:,s.ltbiprevalent,:),2));
    nltb_all = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.ltbiprevalent,s.a0_4),:),2));
    nltb_a0 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.ltbiprevalent,s.a5_9),:),2));
    nltb_a5 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.ltbiprevalent,s.a10_15),:),2));
    nltb_a10 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.ltbiprevalent,s.a15p),:),2));
    nltb_a15p = (tmp(1:end-1,:)+tmp(2:end,:))/2;

    tmp = squeeze(sum(allsol(:,intersect(s.ltbiprevalent,s.hpos),:),2));
    nltb_hpos = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.ltbiprevalent,intersect(s.hpos,s.a0_4)),:),2));
    nltb_hpos_a0 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.ltbiprevalent,intersect(s.hpos,s.a5_9)),:),2));
    nltb_hpos_a5 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.ltbiprevalent,intersect(s.hpos,s.a10_15)),:),2));
    nltb_hpos_a10 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.ltbiprevalent,intersect(s.hpos,s.a15p)),:),2));
    nltb_hpos_a15p = (tmp(1:end-1,:)+tmp(2:end,:))/2;

    tmp = squeeze(sum(allsol(:,intersect(s.ltbiprevalent,s.hart),:),2));
    nltb_hart = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.ltbiprevalent,intersect(s.hart,s.a0_4)),:),2));
    nltb_hart_a0 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.ltbiprevalent,intersect(s.hart,s.a5_9)),:),2));
    nltb_hart_a5 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.ltbiprevalent,intersect(s.hart,s.a10_15)),:),2));
    nltb_hart_a10 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.ltbiprevalent,intersect(s.hart,s.a15p)),:),2));
    nltb_hart_a15p = (tmp(1:end-1,:)+tmp(2:end,:))/2;

    tmp = squeeze(sum(allsol(:,intersect(s.ltbiprevalent,s.slum),:),2));
    nltb_slum = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.ltbiprevalent,intersect(s.slum,s.a0_4)),:),2));
    nltb_slum_a0 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.ltbiprevalent,intersect(s.slum,s.a5_9)),:),2));
    nltb_slum_a5 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.ltbiprevalent,intersect(s.slum,s.a10_15)),:),2));
    nltb_slum_a10 = (tmp(1:end-1,:)+tmp(2:end,:))/2;
    tmp = squeeze(sum(allsol(:,intersect(s.ltbiprevalent,intersect(s.slum,s.a15p)),:),2));
    nltb_slum_a15p = (tmp(1:end-1,:)+tmp(2:end,:))/2;

    %% outputs (years  x age*cat x scenario)

    inc_ds(ii,:,:,:)  = cat(2,sum(diff(squeeze(allsol(:,i.aux.inc_ds(1:4),:)),1),2)*p.pop, diff(squeeze(allsol(:,i.aux.inc_ds,:)),1)*p.pop);
    inc_dr(ii,:,:,:)  = cat(2,sum(diff(squeeze(allsol(:,i.aux.inc_dr(1:4),:)),1),2)*p.pop, diff(squeeze(allsol(:,i.aux.inc_dr,:)),1)*p.pop);
    screen_all(ii,:,:,:)  = cat(2,sum(diff(squeeze(allsol(:,i.aux.screen_all,:)),1),2)*p.pop, diff(squeeze(allsol(:,i.aux.screen_all,:)),1)*p.pop);
    screen_plhiv(ii,:,:,:)  = cat(2,sum(diff(squeeze(allsol(:,i.aux.screen_plhiv,:)),1),2)*p.pop, diff(squeeze(allsol(:,i.aux.screen_plhiv,:)),1)*p.pop);
    screen_hhc(ii,:,:,:)  = cat(2,sum(diff(squeeze(allsol(:,i.aux.screen_hhc,:)),1),2)*p.pop, diff(squeeze(allsol(:,i.aux.screen_hhc,:)),1)*p.pop);
    screen_slum(ii,:,:,:)  = cat(2,sum(diff(squeeze(allsol(:,i.aux.screen_slum,:)),1),2)*p.pop, diff(squeeze(allsol(:,i.aux.screen_slum,:)),1)*p.pop);
    tst_all(ii,:,:,:)  = cat(2,sum(diff(squeeze(allsol(:,i.aux.tst_all,:)),1),2)*p.pop, diff(squeeze(allsol(:,i.aux.tst_all,:)),1)*p.pop);
    tst_plhiv(ii,:,:,:)  = cat(2,sum(diff(squeeze(allsol(:,i.aux.tst_plhiv,:)),1),2)*p.pop, diff(squeeze(allsol(:,i.aux.tst_plhiv,:)),1)*p.pop);
    tst_hhc(ii,:,:,:)  = cat(2,sum(diff(squeeze(allsol(:,i.aux.tst_hhc,:)),1),2)*p.pop, diff(squeeze(allsol(:,i.aux.tst_hhc,:)),1)*p.pop);
    tst_slum(ii,:,:,:)  = cat(2,sum(diff(squeeze(allsol(:,i.aux.tst_slum,:)),1),2)*p.pop, diff(squeeze(allsol(:,i.aux.tst_slum,:)),1)*p.pop);
    tpt_all(ii,:,:,:)  = cat(2,sum(diff(squeeze(allsol(:,i.aux.tpt_all,:)),1),2)*p.pop, diff(squeeze(allsol(:,i.aux.tpt_all,:)),1)*p.pop);
    tpt_plhiv(ii,:,:,:)  = cat(2,sum(diff(squeeze(allsol(:,i.aux.tpt_plhiv,:)),1),2)*p.pop, diff(squeeze(allsol(:,i.aux.tpt_plhiv,:)),1)*p.pop);
    tpt_hhc(ii,:,:,:)  = cat(2,sum(diff(squeeze(allsol(:,i.aux.tpt_hhc,:)),1),2)*p.pop, diff(squeeze(allsol(:,i.aux.tpt_hhc,:)),1)*p.pop);
    tpt_slum(ii,:,:,:)  = cat(2,sum(diff(squeeze(allsol(:,i.aux.tpt_slum,:)),1),2)*p.pop, diff(squeeze(allsol(:,i.aux.tpt_slum,:)),1)*p.pop);
    tx_fl_all(ii,:,:,:)  = cat(2,sum(diff(squeeze(allsol(:,i.aux.tx_fl_all,:)),1),2)*p.pop, diff(squeeze(allsol(:,i.aux.tx_fl_all,:)),1)*p.pop);
    tx_fl_plhiv(ii,:,:,:)  = cat(2,sum(diff(squeeze(allsol(:,i.aux.tx_fl_plhiv,:)),1),2)*p.pop, diff(squeeze(allsol(:,i.aux.tx_fl_plhiv,:)),1)*p.pop);
    tx_fl_hhc(ii,:,:,:)  = cat(2,sum(diff(squeeze(allsol(:,i.aux.tx_fl_hhc,:)),1),2)*p.pop, diff(squeeze(allsol(:,i.aux.tx_fl_hhc,:)),1)*p.pop);
    tx_fl_slum(ii,:,:,:)  = cat(2,sum(diff(squeeze(allsol(:,i.aux.tx_fl_slum,:)),1),2)*p.pop, diff(squeeze(allsol(:,i.aux.tx_fl_slum,:)),1)*p.pop);
    tx_sl_all(ii,:,:,:)  = cat(2,sum(diff(squeeze(allsol(:,i.aux.tx_sl_all,:)),1),2)*p.pop, diff(squeeze(allsol(:,i.aux.tx_sl_all,:)),1)*p.pop);
    tx_sl_plhiv(ii,:,:,:)  = cat(2,sum(diff(squeeze(allsol(:,i.aux.tx_sl_plhiv,:)),1),2)*p.pop, diff(squeeze(allsol(:,i.aux.tx_sl_plhiv,:)),1)*p.pop);
    tx_sl_hhc(ii,:,:,:)  = cat(2,sum(diff(squeeze(allsol(:,i.aux.tx_sl_hhc,:)),1),2)*p.pop, diff(squeeze(allsol(:,i.aux.tx_sl_hhc,:)),1)*p.pop);
    tx_sl_slum(ii,:,:,:)  = cat(2,sum(diff(squeeze(allsol(:,i.aux.tx_sl_slum,:)),1),2)*p.pop, diff(squeeze(allsol(:,i.aux.tx_sl_slum,:)),1)*p.pop);
    % daly_all(ii,:,:,:)    = cat(2,sum(diff(squeeze(allsol(:,i.aux.daly_all,:)),1),2)*p.pop,diff(squeeze(allsol(:,i.aux.daly_all,:)),1)*p.pop);
    % daly_plhiv(ii,:,:,:)    = cat(2,sum(diff(squeeze(allsol(:,i.aux.daly_plhiv,:)),1),2)*p.pop,diff(squeeze(allsol(:,i.aux.daly_plhiv,:)),1)*p.pop);
    % daly_slum(ii,:,:,:)    = cat(2,sum(diff(squeeze(allsol(:,i.aux.daly_slum,:)),1),2)*p.pop,diff(squeeze(allsol(:,i.aux.daly_slum,:)),1)*p.pop);
    yll_all(ii,:,:,:)     = cat(2,sum(diff(squeeze(allsol(:,i.aux.yll_all,:)),1),2)*p.pop,diff(squeeze(allsol(:,i.aux.yll_all,:)),1)*p.pop);
    yll_plhiv(ii,:,:,:)     = cat(2,sum(diff(squeeze(allsol(:,i.aux.yll_plhiv,:)),1),2)*p.pop,diff(squeeze(allsol(:,i.aux.yll_plhiv,:)),1)*p.pop);
    yll_slum(ii,:,:,:)     = cat(2,sum(diff(squeeze(allsol(:,i.aux.yll_slum,:)),1),2)*p.pop,diff(squeeze(allsol(:,i.aux.yll_slum,:)),1)*p.pop);
    yld_all(ii,:,:,:)     = cat(2,sum(diff(squeeze(allsol(:,i.aux.yld_all,:)),1),2)*p.pop,diff(squeeze(allsol(:,i.aux.yld_all,:)),1)*p.pop);
    yld_plhiv(ii,:,:,:)     = cat(2,sum(diff(squeeze(allsol(:,i.aux.yld_plhiv,:)),1),2)*p.pop,diff(squeeze(allsol(:,i.aux.yld_plhiv,:)),1)*p.pop);
    yld_slum(ii,:,:,:)     = cat(2,sum(diff(squeeze(allsol(:,i.aux.yld_slum,:)),1),2)*p.pop,diff(squeeze(allsol(:,i.aux.yld_slum,:)),1)*p.pop);
    mu_ds(ii,:,:,:)   = cat(2,sum(diff(squeeze(allsol(:,i.aux.mu_ds(1:4),:)),1),2)*p.pop,diff(squeeze(allsol(:,i.aux.mu_ds,:)),1)*p.pop);
    mu_dr(ii,:,:,:)   = cat(2,sum(diff(squeeze(allsol(:,i.aux.mu_dr(1:4),:)),1),2)*p.pop,diff(squeeze(allsol(:,i.aux.mu_dr,:)),1)*p.pop);

    daly_all(ii,:,:,:)    = yll_all(ii,:,:,:)+ yld_all(ii,:,:,:);%cat(2,sum(diff(squeeze(allsol(:,i.aux.daly_all,:)),1),2)*p.pop,diff(squeeze(allsol(:,i.aux.daly_all,:)),1)*p.pop);
    daly_plhiv(ii,:,:,:)  = yll_plhiv(ii,:,:,:)+ yld_plhiv(ii,:,:,:);% cat(2,sum(diff(squeeze(allsol(:,i.aux.daly_plhiv,:)),1),2)*p.pop,diff(squeeze(allsol(:,i.aux.daly_plhiv,:)),1)*p.pop);
    daly_slum(ii,:,:,:)   = yll_slum(ii,:,:,:)+ yld_slum(ii,:,:,:);% cat(2,sum(diff(squeeze(allsol(:,i.aux.daly_slum,:)),1),2)*p.pop,diff(squeeze(allsol(:,i.aux.daly_slum,:)),1)*p.pop);



    prevtb(ii,:,:,:)   = permute(cat(3,ntb_all./pop_all,ntb_a0./pop_a0,ntb_a5./pop_a5,ntb_a10./pop_a10, ntb_a15p./pop_a15p),[1 3 2]);
    prevtb_hpos(ii,:,:,:) = permute(cat(3,ntb_hpos./pop_hpos,ntb_hpos_a0./pop_hpos_a0,ntb_hpos_a5./pop_hpos_a5,ntb_hpos_a10./pop_hpos_a10, ntb_hpos_a15p./pop_hpos_a15p),[1 3 2]);
    prevtb_hart(ii,:,:,:) = permute(cat(3,ntb_hart./pop_hart,ntb_hart_a0./pop_hart_a0,ntb_hart_a5./pop_hart_a5,ntb_hart_a10./pop_hart_a10, ntb_hart_a15p./pop_hart_a15p),[1 3 2]);
    prevtb_slum(ii,:,:,:) = permute(cat(3,ntb_slum./pop_slum,ntb_slum_a0./pop_slum_a0,ntb_slum_a5./pop_slum_a5,ntb_slum_a10./pop_slum_a10, ntb_slum_a15p./pop_slum_a15p),[1 3 2]);
    prevltbi(ii,:,:,:)    = permute(cat(3,nltb_all./pop_all,nltb_a0./pop_a0,nltb_a5./pop_a5,nltb_a10./pop_a10, nltb_a15p./pop_a15p),[1 3 2]);
    prevltbi_hpos(ii,:,:,:) = permute(cat(3,nltb_hpos./pop_hpos,nltb_hpos_a0./pop_hpos_a0,nltb_hpos_a5./pop_hpos_a5,nltb_hpos_a10./pop_hpos_a10, nltb_hpos_a15p./pop_hpos_a15p),[1 3 2]);
    prevltbi_hart(ii,:,:,:) = permute(cat(3,nltb_hart./pop_hart,nltb_hart_a0./pop_hart_a0,nltb_hart_a5./pop_hart_a5,nltb_hart_a10./pop_hart_a10, nltb_hart_a15p./pop_hart_a15p),[1 3 2]);
    prevltbi_slum(ii,:,:,:)= permute(cat(3,nltb_slum./pop_slum,nltb_slum_a0./pop_slum_a0,nltb_slum_a5./pop_slum_a5,nltb_slum_a10./pop_slum_a10, nltb_slum_a15p./pop_slum_a15p),[1 3 2]);


    % Alive

    alive_all(ii,:,:,:) =  permute(cat(3,pop_all, pop_a0,pop_a5,pop_a10,pop_a15p)*p.pop,[1 3 2]);
    alive_hpos(ii,:,:,:) =  permute(cat(3,pop_hpos,pop_hpos_a0,pop_hpos_a5,pop_hpos_a10,pop_hpos_a15p)*p.pop,[1 3 2]);
    alive_hart(ii,:,:,:) =  permute(cat(3,pop_hart,pop_hart_a0,pop_hart_a5,pop_hart_a10, pop_hart_a15p)*p.pop,[1 3 2]);
    alive_slum(ii,:,:,:) =  permute(cat(3,pop_slum,pop_slum_a0,pop_slum_a5,pop_slum_a10, pop_slum_a15p)*p.pop,[1 3 2]);
    alive_hhc(ii,:,:,:) = permute(cat(3,ntb_all,ntb_a0,ntb_a5,ntb_a10, ntb_a15p).*p.pop*(p.house_size-1),[1 3 2]);

  

    %% Usual Epi trends

    tmp1=sum(diff(squeeze(allsol(:,i.aux.inc_ds(1:4),:)),1),2);
    tmp2=sum(diff(squeeze(allsol(:,i.aux.inc_dr(1:4),:)),1),2);
    inc_all(ii,:,:)= 1e5*squeeze(tmp1+tmp2)./pop_all;
    inc_mdr(ii,:,:)= 1e5*squeeze(tmp2)./pop_all;

    tmp1=sum(diff(squeeze(allsol(:,i.aux.mu_ds(1:4),:)),1),2);
    tmp2=sum(diff(squeeze(allsol(:,i.aux.mu_dr(1:4),:)),1),2);
    mort_tb(ii,:,:)= 1e5*squeeze(tmp1+tmp2)./pop_all;


    sfins(ii,:,:)=allsol(end,:,:);

disp(ii);
end

delete(gcp('nocreate'));

results.prevtb=prevtb;
results.prevtb_hpos=prevtb_hpos;
results.prevtb_hart=prevtb_hart;
results.prevtb_slum=prevtb_slum;
results.prevltbi=prevltbi;
results.prevltbi_hpos=prevltbi_hpos;
results.prevltbi_hart=prevltbi_hart;
results.prevltbi_slum=prevltbi_slum;
results.inc_ds=inc_ds;
results.inc_dr=inc_dr;
results.screen_all=screen_all;
results.screen_plhiv=screen_plhiv;
results.screen_hhc=screen_hhc;
results.screen_slum=screen_slum;

results.tst_all=tst_all;
results.tst_plhiv=tst_plhiv;
results.tst_hhc=tst_hhc;
results.tst_slum=tst_slum;

results.tpt_all=tpt_all;
results.tpt_plhiv=tpt_plhiv;
results.tpt_hhc=tpt_hhc;
results.tpt_slum=tpt_slum;

results.tx_fl_all=tx_fl_all;
results.tx_fl_plhiv=tx_fl_plhiv;
results.tx_fl_hhc=tx_fl_hhc;
results.tx_fl_slum=tx_fl_slum;

results.tx_sl_all=tx_sl_all;
results.tx_sl_plhiv=tx_sl_plhiv;
results.tx_sl_hhc=tx_sl_hhc;
results.tx_sl_slum=tx_sl_slum;

results.daly_all=daly_all;
results.daly_plhiv=daly_plhiv;
results.daly_slum=daly_slum;
results.yll_all=yll_all;
results.yll_plhiv=yll_plhiv;
results.yll_slum=yll_slum;
results.yld_all=yld_all;
results.yld_plhiv=yld_plhiv;
results.yld_slum=yld_slum;
results.mu_ds=mu_ds;
results.mu_dr=mu_dr;
results.avagemort=avagemort;
results.avagemort_plhiv=avagemort_plhiv;
results.avagemort_slum=avagemort_slum;
results.inc_all=inc_all;
results.inc_mdr=inc_mdr;
results.mort_tb=mort_tb;
results.alive_all=alive_all;
results.alive_hpos=alive_hpos;
results.alive_hart=alive_hart;
results.alive_slum=alive_slum;
results.alive_hhc=alive_hhc;


results.sfin=sfins;
results.s=s;

results.itv_labs=itv_labs(path);

txt=sprintf('%s','itvnoTST','_',num2str(scen));

save_results(results,txt,location,runtype);
