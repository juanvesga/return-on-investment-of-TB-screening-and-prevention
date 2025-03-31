%% This Funtions constructs Model structure with trasnitions (linear and non-linear)

function M = make_model3(p, r, i, s, gps)

% --- Structures to hold matrices
% --- Get the linear rates ------------------------------------------------
m  = zeros(i.nstates);     % Matrix for linear transitions in the model (i.e. except infection)
m2 = zeros(i.nstates);    % slum Linear
m3 = zeros(i.nstates);   % Linear components HIV cascade
m4 = zeros(i.nstates);   % force of transition for HHC case finding U
m8 = zeros(i.nstates);   % matrix for counting screening on plhiv
m9 = zeros(i.nstates);   % matrix for counting screening on hhcU
m13 = zeros(i.nstates);   % matrix for counting screening on slum
m14 = zeros(i.nstates);   % Linear components of HIV tx and prevention


fot = zeros(1,i.nstates);% Elements in the force of transition for CFY

% --- Apply diagnostic transitions to NTS and Preventive therapy for all model
for ia = 1:length(gps.age)
    age = gps.age{ia};

    for ir = 1:length(gps.georisk)
        risk = gps.georisk{ir};

        islum = (ir > 1);
        isnoslum = (ir == 1);

        for ig = 1:length(gps.hiv)
            hiv = gps.hiv{ig};
            getind1 = @(st) i.(st).(age).(risk).(hiv);

            ishneg = (ig ==1);

            U    = getind1('U');
            Pu   = getind1('Pu');
            Qu   = getind1('Qu');


            % Case finding (PLHIV)
            source = U; destin = Pu;
            ch=1;%p.lam_ds(ia)+p.lam_dr(ia);

            prop=p.case_finding_U(ia,ir,ig)*p.hiv_cov_tpt;
            rate=ch*prop/(1-prop)* isnoslum * p.plhiv_boost;

            % rate = p.case_finding_U(ia,ir,ig)*p.hiv_cov_tpt;
            m14(destin, source) = m14(destin, source) + rate;


            %% HHC Case finding

            r_u=get_HHC_rate(2,p,ia,ch,'U',p.U_age_hhc,'tpt');

            prop=p.case_findingHHC_U(ia,ir,ig) ;

            rate=r_u * prop* r.hhc_scale(ia)*p.hhc_boost*p.hhcUfac;

            m4(destin, source) = m4(destin, source) + rate;


            r_u=get_HHC_rate(2,p,ia,ch,'U',p.U_age_hhc,'screen');

            prop = 1 ;
            
            rate =  r_u * prop* r.hhc_scale(ia) * p.hhc_screen_dim *p.hhc_boost*p.hhcUfac;

            m9(destin, source) = m9(destin, source) + rate;

            % Case finding (slums)
            source = U; destin = Pu;
            prop=p.case_findingslum_U(ia,ir,ig)*p.slum_cov_tpt;
            rate= ishneg *p.slum_enhance*ch*prop/(1-prop)*p.slum_boost ;
            m2(destin, source) = m2(destin, source) + rate;

            prop=p.slum_cov_screen;
            
            rate=ishneg *p.slum_enhance*islum*ch*prop/(1-prop)*p.slum_boost*p.slum_screen_dim; 
            
            m13(destin, source) = m13(destin, source) + rate;

            % Regimen completion
            r.outReg        =  12/((r.tpt_long_reg_dur_mo)*(1-p.tpt_short(ir,ig)) +...
                (r.tpt_short_reg_dur_mo)*p.tpt_short(ir,ig));

           tpt_efficacy=p.tpt_eff_plhiv(ig)*(p.tpt_short(ir,ig)*p.pteffi_short + (1-p.tpt_short(ir,ig))*p.pteffi);


            p.tpt_complete(1) = (p.tpt_base_comp*(1-p.tpt_short(ir,ig)) + p.tpt_short(ir,ig)/2*p.tpt_target_comp + p.tpt_short(ir,ig)/2*p.tpt_target_compa0);
            p.tpt_complete(2) = (p.tpt_base_comp*(1-p.tpt_short(ir,ig)) + p.tpt_short(ir,ig)*p.tpt_target_comp);
            p.tpt_complete(3) = (p.tpt_base_comp*(1-p.tpt_short(ir,ig)) + p.tpt_short(ir,ig)*p.tpt_target_comp);
            p.tpt_complete(4) = (p.tpt_base_comp*(1-p.tpt_short(ir,ig)) + p.tpt_short(ir,ig)*p.tpt_target_compa15p);

            r.tpt_default  = r.outReg.*(1-p.tpt_complete)./p.tpt_complete;
            source = Pu; destin = Qu; rate = r.outReg;
            m(destin, source) = m(destin, source) + rate;

            % Regimen default
            source = Pu; destin = U; rate = r.tpt_default(ia);
            m(destin, source) = m(destin, source) + rate;

            for is = 1:length(gps.strain)
                strain = gps.strain{is};

                gi = @(st) i.(st).(age).(risk).(hiv).(strain);
                Lf    = gi('Lf');
                Ls    = gi('Ls');
                Pf    = gi('Pf');
                Ps    = gi('Ps');
                Qf    = gi('Qf');
                Qs    = gi('Qs');
                I     = gi('I');
                Dx    = gi('Dx');
                Tx    = gi('Tx');
                Tx2   = gi('Tx2');
                Rhi   = gi('Rhi');
                R     = gi('R');

                ismdr   = strcmp(strain, 'mdr');

                %% Pf to Ls
                source = Lf; destin = Pf;

                competing_hazard= r.progression(ir,ig)*r.RRprog_age(ia)+r.LTBI_stabil(ir,ig);
              
                %~~~~ Case finding (PLHIV)  
                prop=p.case_finding_Lf(ia,ir,ig)*p.hiv_cov_tpt;
                
                rate=isnoslum * competing_hazard*prop/(1-prop)* p.plhiv_boost;
                
                m14(destin, source) = m14(destin, source) + rate;

                %~~~~ HHC Case finding
                if(ismdr==0) ;st='Lfds'; else; st='Lfdr'; end
                r_lf=get_HHC_rate(ismdr,p,ia,competing_hazard,st,p.Lf_age_hhc,'tpt');

                prop= p.case_findingHHC_Lf(ia,ir,ig);

                rate=r_lf * prop* r.hhc_scale(ia)*p.hhc_boost*p.hhcLffac;
                
                
                m4(destin, source) = m4(destin, source) + rate;


                r_lf=get_HHC_rate(ismdr,p,ia,competing_hazard,st,p.Lf_age_hhc,'screen');

                prop= 1;
                
                rate= r_lf * prop* r.hhc_scale(ia)*p.hhc_boost*p.hhc_screen_dim*p.hhcLffac;

                m9(destin, source) = m9(destin, source) + rate;

                %~~~~ Slum Case finding
                prop= p.case_findingslum_Lf(ia,ir,ig)*p.slum_cov_tpt;

                rate=ishneg * p.slum_enhance*competing_hazard*prop/(1-prop)* p.slum_boost ;
                
                m2(destin, source) = m2(destin, source) + rate;


                prop=p.slum_cov_screen;
                
                rate=ishneg * islum*competing_hazard*prop/(1- prop)* p.slum_boost*p.slum_screen_dim; 
                
                m13(destin, source) = m13(destin, source) + rate;

                %% Ps to Ls
                source = Ls; destin = Ps;

                competing_hazard= r.reactivation(ir,ig);
                
                %~~~~ Case finding (PLHIV)  
                prop=p.case_finding_Ls(ia,ir,ig)*p.hiv_cov_tpt;
                
                rate = isnoslum * competing_hazard*prop/(1-prop)* p.plhiv_boost;
                
                m14(destin, source) = m14(destin, source) + rate;

                %~~~~ HHC Case finding    
                if(ismdr==0); st='Lsds'; else; st='Lsdr'; end
                
                r_ls=get_HHC_rate(ismdr,p,ia,competing_hazard,st,p.Ls_age_hhc,'tpt');

                prop= p.case_findingHHC_Ls(ia,ir,ig);

                rate=r_ls * prop* r.hhc_scale(ia)*p.hhc_boost*p.hhcLsfac;
                
                m4(destin, source) = m4(destin, source) + rate;


                r_ls=get_HHC_rate(ismdr,p,ia,competing_hazard,st,p.Ls_age_hhc,'screen');

                prop= 1;
                
                rate= r_ls * prop* r.hhc_scale(ia)*p.hhc_boost*p.hhc_screen_dim*p.hhcLsfac;

                m9(destin, source) = m9(destin, source) + rate;


                %~~~~ Slum Case finding
                prop= p.case_findingslum_Ls(ia,ir,ig)*p.slum_cov_tpt;

                rate= ishneg*p.slum_enhance*competing_hazard*prop/(1-prop)* p.slum_boost ;
                
                m2(destin, source) = m2(destin, source) + rate;

                prop=p.slum_cov_screen;
                
                rate=ishneg*islum*p.slum_enhance*competing_hazard*prop/(1-prop)* p.slum_boost*p.slum_screen_dim; 
                
                m13(destin, source) = m13(destin, source) + rate;

                % Regimen default
                source = Pf; destin = Lf;
                rate =  r.tpt_default(ia);
                m(destin, source) = m(destin, source) + rate;

                source = Ps; destin = Ls;
                rate = r.tpt_default(ia);
                m(destin, source) = m(destin, source) + rate;



                % Regimen completion
                if ismdr~=1
                    source = Pf; destin = Qf; rate = r.outReg*(1-p.potency);
                    m(destin, source) = m(destin, source) + rate;

                    source = Ps; destin = Qs; rate = r.outReg*(1-p.potency);
                    m(destin, source) = m(destin, source) + rate;

                    source = Pf; destin = Qu; rate = r.outReg*(p.potency);
                    m(destin, source) = m(destin, source) + rate;

                    source = Ps; destin = Qu; rate = r.outReg*(p.potency);
                    m(destin, source) = m(destin, source) + rate;

                else
                    source = Pf; destin = Qf; rate = r.outReg*(1-p.potency)*p.rif_free;
                    m(destin, source) = m(destin, source) + rate;

                    source = Ps; destin = Qs; rate = r.outReg*(1-p.potency)*p.rif_free;
                    m(destin, source) = m(destin, source) + rate;

                    source = Pf; destin = Qu; rate = r.outReg*(p.potency)*p.rif_free;
                    m(destin, source) = m(destin, source) + rate;

                    source = Ps; destin = Qu; rate = r.outReg*(p.potency)*p.rif_free;
                    m(destin, source) = m(destin, source) + rate;

                    source = Pf; destin = Qf; rate = r.outReg*(1-p.rif_free)*2;
                    m(destin, source) = m(destin, source) + rate;

                    source = Ps; destin = Qs; rate = r.outReg*(1-p.rif_free)*2;
                    m(destin, source) = m(destin, source) + rate;


                end

                % --- Fast progression and LTBI stabilisation
                source  = Lf;
                destins = [I,                 Ls];
                rates   = [r.progression(ir,ig)*r.RRprog_age(ia), r.LTBI_stabil(ir,ig)];
                m(destins, source) = m(destins, source) + rates';

                % --- Reactivation
                source = Ls; destin = I; rate = r.reactivation(ir,ig);
                m(destin, source) = m(destin, source) + rate;

                % --- Reactivation with PT effect
               
                source = Pf; 
                destin = [I,  Ps];
                rates = [r.progression(ir,ig)*r.RRprog_age(ia)*ismdr*(1-tpt_efficacy*p.effondr) +...
                    r.progression(ir,ig)*r.RRprog_age(ia)*(1-ismdr)*(1-tpt_efficacy),...
                    r.LTBI_stabil(ir,ig)];
                m(destin, source) = m(destin, source) + rates';

                source  = Qf;
                destins = [I,                 Qs];
                m(destins, source) = m(destins, source) + rates';


                source = Ps; destin = I;
                rate = r.reactivation(ir,ig)*ismdr*(1-tpt_efficacy*p.effondr) +...
                    r.reactivation(ir,ig)*(1-ismdr)*(1-tpt_efficacy);
                m(destin, source) = m(destin, source) + rate;

                source = Qs; destin = I; 
                m(destin, source) = m(destin, source) + rate;


                % Careseeking

                source = I; destin = Dx;
                rate= r.careseeking*r.cseek_fac(ig)*r.cseek_fac_slum(ir);
                m(destin, source) = m(destin, source) + rate;


                % Primary casefinding
                p_firstline = p.Tx_init * ( ...
                    ( (1-ismdr)*p.xpert*p.xpert_sens) +...
                    ( (1-p.xpert)*p.smear_sens ) ...
                    );

                p_secondline=p.Tx_init2*ismdr*p.xpert*p.xpert_sens;

                p_ltfu  = 1-( p.Tx_init * ( ...
                    ( (1-ismdr)*p.xpert*p.xpert_sens) +...
                    ( (1-p.xpert)*p.smear_sens ) ...
                    ) + ...
                    (...
                    p.Tx_init2*ismdr*p.xpert*p.xpert_sens...
                    )...
                    );

                source = Dx;
                destin = [Tx,       Tx2,             I];
                rate   = r.Dx*[p_firstline p_secondline p_ltfu];
                m(destin, source) = m(destin, source) + rate';


                %% FOT
                rate=r.Dx*(p_firstline+p_secondline);
                sw=p.cfy>0;
                fot(1,I)= rate*sw;


                %% Screening I to Tx and Tx2
                source = I;
                destin = [Tx  Tx2];
                competing_hazard= r.careseeking*r.cseek_fac(ig)*r.cseek_fac_slum(ir);
                
                %~~~~~~ Case finding (PLHIV)
                prop=p.case_finding_I(ia,ir,ig)*p.hiv_cov_screen;
                rate= isnoslum * competing_hazard*prop/(1-prop)* p.plhiv_boost;
                rates= [rate*(1-ismdr) , rate*ismdr];

                m14(destin, source) = m14(destin, source) + rates';


                %~~~~~~ Case finding HHC Case finding
                if(ismdr==0); st='Ids'; else; st='Idr'; end
                r_i=get_HHC_rate(ismdr,p,ia,competing_hazard,st,p.I_age_hhc,'screen');

                prop= p.case_findingHHC_I(ia,ir,ig) * p.base_hhc_fac;

                rates=[...
                    r_i * prop * (1-ismdr) *...
                    r.hhc_scale(ia)*p.hhc_boost*p.hhcIfac,...
                    r_i * prop * (ismdr) *...
                    r.hhc_scale(ia)*p.hhc_boost*p.hhcIfac];
                m4(destin, source) = m4(destin, source) + rates';

                r_i=get_HHC_rate(ismdr,p,ia,competing_hazard,st,p.I_age_hhc,'screen');

                prop=1;

                rates=[...
                    r_i * prop * (1-ismdr) *...
                    r.hhc_scale(ia)*p.hhc_boost*p.hhcIfac,...
                    r_i * prop * (ismdr) *...
                    r.hhc_scale(ia)*p.hhc_boost*p.hhcIfac];

                m9(destin, source) = m9(destin, source) + rates';


                %~~~~~~ Case finding (slum)
                prop=p.case_findingslum_I(ia,ir,ig)*p.slum_cov_screen;
                rate = ishneg * competing_hazard*prop/(1-prop)* p.slum_boost ;
                rates= [rate*(1-ismdr) , rate*ismdr ];

                m2(destin, source) = m2(destin, source) + rates';


                prop= p.slum_cov_screen;
                rate =ishneg * islum * competing_hazard*prop/(1-prop)* p.slum_boost ;
                rates= [rate*(1-ismdr) , rate*ismdr];
        
                m13(destin, source) = m13(destin, source) + rates';



                % --- Treatment outcomes
                r.Tx          = [12/6 12/(18*(1-p.sl_short) + 6*p.sl_short)];
                p.Tx_complete = [p.fl_suc p.sl_short*0.89 + (1-p.sl_short)*0.6];
                r.default     = r.Tx.*(1-p.Tx_complete)./p.Tx_complete;
                if ismdr~=1

                    source  = Tx;
                    destins = [R  ,  Rhi ,i.Tx.(age).(risk).(hiv).('mdr')];
                    rates   = [r.Tx(1),  r.default(1)*(1-r.MDR_acqu) r.default(1)*r.MDR_acqu];
                    m(destins, source) = m(destins, source) + rates';
                else

                    source  = Tx;
                    destins = [Tx2             Rhi];
                    rates   = [r.Tx(1)*p.SL_trans,   r.Tx(1)*(1-p.SL_trans) + r.default(1)];
                    m(destins, source) = m(destins, source) + rates';

                end


                source  = Tx2;
                destins = [R         Rhi      ];
                rates   = [r.Tx(2),    r.default(2)];
                m(destins, source) = m(destins, source) + rates';

                % --- Relapse
                sources = [Rhi, R];
                destin  = Dx;
                rates   = r.relapse;
                m(destin, sources) = m(destin, sources) + rates;


                sources = Rhi;
                destin  = R;
                rates   = 0.5;
                m(destin, sources) = m(destin, sources) + rates;

                % --- Self cure
                rate = r.self_cure(ig);

                source = I;
                destin = Rhi;
                m(destin, source) = m(destin, source) + rate;


                source = Dx;
                destin = Rhi;
                m(destin, source) = m(destin, source) + rate;



            end
        end

    end
end



%% --- Ageing
sources = s.a0_4;
destins = s.a5_9;
rates   = r.ageing(1);
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
m(inds) = m(inds) + rates;


sources = s.a5_9;
destins = s.a10_15;
rates   = r.ageing(2);
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
m(inds) = m(inds) + rates;

sources = s.a10_15;
destins = s.a15p;
rates   = r.ageing(3);
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
m(inds) = m(inds) + rates;



%% HIV cascade
% HIV acquisition
sources = intersect(s.a0_4,intersect(s.hneg, s.no_slum));
destins = intersect(s.a0_4,intersect(s.hpos, s.no_slum));
rates   = p.hivcq_ad(1);
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
m3(inds) = m3(inds) + rates;

sources = intersect(s.a5_9,intersect(s.hneg, s.no_slum));
destins = intersect(s.a5_9,intersect(s.hpos, s.no_slum));
rates   = p.hivcq_ad(2);
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
m3(inds) = m3(inds) + rates;

sources = intersect(s.a10_15,intersect(s.hneg, s.no_slum));
destins = intersect(s.a10_15,intersect(s.hpos, s.no_slum));
rates   = p.hivcq_ad(3);
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
m3(inds) = m3(inds) + rates;

sources = intersect(s.a15p,intersect(s.hneg, s.no_slum));
destins = intersect(s.a15p,intersect(s.hpos, s.no_slum));
rates   = p.hivcq_ad(4);
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
m3(inds) = m3(inds) + rates;



%% Interventions among PLHIV

for ia = 1:length(gps.age)
    age = gps.age{ia};

    hpos=intersect(s.hpos, s.no_slum);
    hart=intersect(s.hart, s.no_slum);

    % U tp U
    sources = intersect(intersect(hpos,s.U),s.(age));
    destins = intersect(intersect(hart,s.U),s.(age));
    rates   =  r.ART_init *...
        (1 - (p.IPThiv + (1-p.IPThiv) * p.hiv_cov_tpt * p.case_findingPLHIV_U(ia)));
    inds    = sub2ind([i.nstates, i.nstates],destins,sources);
    m(inds) = m(inds) + rates;

    % U to Pu
    sources = intersect(intersect(hpos,s.U),s.(age));
    destins = intersect(intersect(hart,s.Pu),s.(age));
    rates   =  r.ART_init *...
        ( p.IPThiv + (1-p.IPThiv) * p.hiv_cov_tpt * p.case_findingPLHIV_U(ia));
    inds    = sub2ind([i.nstates, i.nstates],destins,sources);
    m(inds) = m(inds) + rates;

    % count screen
    rates   =  r.ART_init *...
        ( p.IPThiv + (1-p.IPThiv) * p.hiv_cov_screen);
    inds    = sub2ind([i.nstates, i.nstates],destins,sources);
    m8(inds) = m8(inds) + rates;

    % Lf to :Lf
    sources = intersect(intersect(hpos,s.Lf),s.(age));
    destins = intersect(intersect(hart,s.Lf),s.(age));
    rates   =  r.ART_init *...
        (1 - (p.IPThiv + (1-p.IPThiv) * p.hiv_cov_tpt * p.case_findingPLHIV_Lf(ia)));
    inds    = sub2ind([i.nstates, i.nstates],destins,sources);
    m(inds) = m(inds) + rates;

    % Lf to P
    sources = intersect(intersect(hpos,s.Lf),s.(age));
    destins = intersect(intersect(hart,s.Pf),s.(age));
    rates   =  r.ART_init *...
        ( p.IPThiv + (1-p.IPThiv) * p.hiv_cov_tpt * p.case_findingPLHIV_Lf(ia));
    inds    = sub2ind([i.nstates, i.nstates],destins,sources);
    m(inds) = m(inds) + rates;

    % count screen
    rates   =  r.ART_init *...
        ( p.IPThiv + (1-p.IPThiv) * p.hiv_cov_screen );
    inds    = sub2ind([i.nstates, i.nstates],destins,sources);
    m8(inds) = m8(inds) + rates;

    % Ls to Ls
    sources = intersect(intersect(hpos,s.Ls),s.(age));
    destins = intersect(intersect(hart,s.Ls),s.(age));
    rates   =  r.ART_init *...
        (1 - (p.IPThiv + (1-p.IPThiv) * p.hiv_cov_tpt * p.case_findingPLHIV_Ls(ia)));
    inds    = sub2ind([i.nstates, i.nstates],destins,sources);
    m(inds) = m(inds) + rates;

    % Ls to P
    sources = intersect(intersect(hpos,s.Ls),s.(age));
    destins = intersect(intersect(hart,s.Ps),s.(age));
    rates   =  r.ART_init *...
        ( p.IPThiv + (1-p.IPThiv) * p.hiv_cov_tpt * p.case_findingPLHIV_Ls(ia));
    inds    = sub2ind([i.nstates, i.nstates],destins,sources);
    m(inds) = m(inds) + rates;

    % Count screen
    rates   =  r.ART_init *...
        ( p.IPThiv + (1-p.IPThiv) * p.hiv_cov_screen);
    inds    = sub2ind([i.nstates, i.nstates],destins,sources);
    m8(inds) = m8(inds) + rates;

    % I to I
    sources = intersect(intersect(hpos,s.I),s.(age));
    destins = intersect(intersect(hart,s.I),s.(age));
    rates   = r.ART_init *...
    (1- p.case_findingPLHIV_I(ia) * p.hiv_cov_screen);
    inds    = sub2ind([i.nstates, i.nstates],destins,sources);
    m(inds) = m(inds) + rates;

    % Dx to Dx
    sources = intersect(intersect(hpos,s.Dx),s.(age));
    destins = intersect(intersect(hart,s.Dx),s.(age));
    rates   = r.ART_init *...
    (1- p.case_findingPLHIV_I(ia) * p.hiv_cov_screen);
    inds    = sub2ind([i.nstates, i.nstates],destins,sources);
    m(inds) = m(inds) + rates;

    % I to Tx
    sources = intersect(intersect(intersect(hpos,s.I),s.(age)),s.ds);
    destins = intersect(intersect(intersect(hart,s.Tx),s.(age)),s.ds);
    rates   = r.ART_init *...
    p.case_findingPLHIV_I(ia) * p.hiv_cov_screen;
    inds    = sub2ind([i.nstates, i.nstates],destins,sources);
    m(inds) = m(inds) + rates;

    % count screen
    rates   = r.ART_init *p.hiv_cov_screen;
    inds    = sub2ind([i.nstates, i.nstates],destins,sources);
    m8(inds) = m8(inds) + rates;

    % Dx to Tx
    % sources = intersect(intersect(intersect(s.hpos,s.Dx),s.(age)),s.ds);
    % destins = intersect(intersect(intersect(s.hart,s.Tx),s.(age)),s.ds);
    % rates   = p.case_findingPLHIV_I(ia) * p.hiv_cov_screen;
    % inds    = sub2ind([i.nstates, i.nstates],destins,sources);
    % m(inds) = m(inds) + rates;
    % 
    % % count screen
    % rates   = p.hiv_cov_screen;
    % inds    = sub2ind([i.nstates, i.nstates],destins,sources);
    % m8(inds) = m8(inds) + rates;


    % I to Tx2
    sources = intersect(intersect(intersect(hpos,s.I),s.(age)),s.mdr);
    destins = intersect(intersect(intersect(hart,s.Tx2),s.(age)),s.mdr);
    rates   = r.ART_init*p.case_findingPLHIV_I(ia) * p.hiv_cov_screen;
    inds    = sub2ind([i.nstates, i.nstates],destins,sources);
    m(inds) = m(inds) + rates;

    % count screen
    rates   = r.ART_init*p.hiv_cov_screen;
    inds    = sub2ind([i.nstates, i.nstates],destins,sources);
    m8(inds) = m8(inds) + rates;

    % Dx to Tx2
    % sources = intersect(intersect(intersect(s.hpos,s.Dx),s.(age)),s.mdr);
    % destins = intersect(intersect(intersect(s.hart,s.Tx2),s.(age)),s.mdr);
    % rates   = p.case_findingPLHIV_I(ia) * p.hiv_cov_screen;
    % inds    = sub2ind([i.nstates, i.nstates],destins,sources);
    % m(inds) = m(inds) + rates;
    % 
    % 
    % % count screen
    % rates   = p.hiv_cov_screen;
    % inds    = sub2ind([i.nstates, i.nstates],destins,sources);
    % m8(inds) = m8(inds) + rates;

end

%------- slum/no_slum turnover (applies only for Brazil prison)

no_slum=intersect(s.hneg, s.no_slum);
slum=intersect(s.hneg, s.slum);


sources = intersect(s.a15p,no_slum);
destins = intersect(s.a15p,slum);
rates   = p.entry_highrisk;
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
m(inds) = m(inds) + rates;

sources = intersect(s.a15p,slum);
destins = intersect(s.a15p,no_slum);
rates   = p.out_highrisk;
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
m(inds) = m(inds) + rates;





% --- Bring them together
M.lin    = sparse(m - diag(sum(m,1)));
M.slumlin  = sparse(m2 - diag(sum(m2,1)));
M.linHIV   = sparse(m14 - diag(sum(m14,1)));
M.linHHC = sparse(m4 - diag(sum(m4,1)));
M.screen_hiv = sparse(m8 - diag(sum(m8,1)));
M.screen_hhc = sparse(m9 - diag(sum(m9,1)));
M.screen_slum = sparse(m13 - diag(sum(m13,1)));
M.nlinHIV = sparse(m3 - diag(sum(m3,1)));
M.fot = fot;

% --- Get the nonlinear rates
for ia = 1:length(gps.age)
    age = gps.age{ia};


    for istr = 1:length(gps.strain)
        strain = gps.strain{istr};

        m = zeros(i.nstates);


        for ir = 1:length(gps.georisk)

            risk = gps.georisk{ir};



            for ig = 1:length(gps.hiv)

                hiv = gps.hiv{ig};


                getso = @(st) [i.(st).(age).(risk).(hiv).('ds') ,...
                    i.(st).(age).(risk).(hiv).('mdr') ];

                getdes = @(st) i.(st).(age).(risk).(hiv).(strain);

                U = i.U.(age).(risk).(hiv);
                Qu = i.U.(age).(risk).(hiv);

                Qfde  = getdes('Qf');
                Qfso  = getso('Qf');
                Lfde  = getdes('Lf');
                Lfso  = getso('Lf');
                Ls  = getso('Ls');
                Qs  = getso('Qs');
                Rhi = getso('Rhi');
                R   = getso('R');

                % m(Lfde, [U Ls Lfso Rhi R]) = 1;
                % m(Qfde, [Qu Qfso Qs]) = 1;%(1-tpt_efficacy);%1;

                m(Lfde, [U Ls Qu Qfso Qs Rhi R]) = 1*r.RR_ari_ch(ia);
                
                % Adjust for imunuty
                cols = [Lfso, Ls, Qfso, Qs, Rhi, R];
                m(:,cols) = m(:,cols)*(1-p.imm(ig));
            end

        end
        m(:,[s.hpos s.hart]) = m(:,[s.hpos s.hart])*p.HIVlam;
        m(:,s.slum) = m(:,s.slum)*(r.RRslum*r.slum_on + (1-r.slum_on));
        M.nlin.(age).(strain) = sparse(m - diag(sum(m,1)));
     end
end



% --- Getting force-of-infection for children and adults settings
%1) 0_4-low
%2) 0_4-high
%3) 5_9-low
%4) 5_9-high
%5) 10_15-low
%6) 10_15-high
%7) 15p-low
%8) 15p-high


% --- Getting force-of-infection for DS and DR
% DS


m_ds = zeros(4,i.nstates);
m_dr = zeros(4,i.nstates);
for ii = 1:length(gps.age)

    m_ds(ii,intersect(s.infectious,s.ds))  = r.beta*p.transmission(ii);
    m_dr(ii,intersect(s.infectious,s.mdr)) = r.beta_mdr*p.transmission(ii);

end
%Strain specific beta
m_ds(:,intersect(s.infectious,s.ch))  = m_ds(:,intersect(s.infectious,s.ch))*r.bRRpaed;
m_dr(:,intersect(s.infectious,s.ch))  = m_dr(:,intersect(s.infectious,s.ch))*r.bRRpaed;
m_ds(:,intersect(s.infectious,s.hivpositive))=m_ds(:,intersect(s.infectious,s.hivpositive)).*r.bRRhiv;
m_dr(:,intersect(s.infectious,s.hivpositive))=m_dr(:,intersect(s.infectious,s.hivpositive)).*r.bRRhiv;



M.lambda_ds = m_ds;
M.lambda_dr = m_dr;

% M.lambda_ds = sparse(m_ds);
% M.lambda_dr = sparse(m_dr);


%% -- - Get the mortality rates- for fit
% First column is for non-TB-related mortality

mfac=r.mfac;

m = zeros(i.nstates,3);
m(s.a0_4,1)   = r.mort(1);
m(s.a5_9,1)   = r.mort(2);
m(s.a10_15,1) = r.mort(3);
m(s.a15p,1)   = r.mort(4);

% hiv positive
m(s.hpos,1) = r.mort_h*r.HIV_mort;
% add postTB RR
m(s.postTByll,1) = m(s.postTByll,1)*1.1;

% Second column is for HIV-ve TB mortality
inds = intersect(s.hneg, s.infectious);
m(inds,2) = r.mort_TB(1);
inds = intersect(s.hneg, s.treated);
m(inds,2) = r.mort_TB(1)*mfac;

% Third column is for HIV+ve TB mortality
inds = intersect(s.hpos, s.infectious);
m(inds,3) = r.mort_TB(2);
inds = intersect(s.hpos, s.treated);
m(inds,2) = r.mort_TB(2)*mfac;

inds = intersect(s.hart, s.infectious);
m(inds,3) = r.mort_TB(3);
inds = intersect(s.hart, s.treated);
m(inds,2) = r.mort_TB(3)*mfac;


M.mortvec = m;

%% Mortlity for model output and Dalys
m = zeros(1,i.nstates);
m(s.a0_4) = r.mort(1);
m(s.a5_9) = r.mort(2);
m(s.a10_15) = r.mort(3);
m(s.a15p) = r.mort(4);

% hiv positive
m(s.hpos) = r.mort_h*r.HIV_mort;

% add postTB RR
m(s.postTByll) = m(s.postTByll)*1.1;

% Second column is for HIV-ve TB mortality
inds = intersect(s.hneg, s.infectious);
m(inds) = r.mort_TB(1);
inds = intersect(s.hneg, s.treated);
m(inds) = r.mort_TB(1)*mfac;

% Third column is for HIV+ve TB mortality
inds = intersect(s.hpos, s.infectious);
m(inds) = r.mort_TB(2);
inds = intersect(s.hpos, s.treated);
m(inds) = r.mort_TB(2)*mfac;

inds = intersect(s.hart, s.infectious);
m(inds) = r.mort_TB(3);
inds = intersect(s.hart, s.treated);
m(inds) = r.mort_TB(3)*mfac;

M.mu = m';
M.p=p;
M.r=r;
