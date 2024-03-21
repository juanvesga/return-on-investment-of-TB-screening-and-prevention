%% This Funtions constructs Model structure with trasnitions (linear and non-linear)

function M = make_model3(p, r, i, s, gps)

% --- Structures to hold matrices
% --- Get the linear rates ------------------------------------------------
m  = zeros(i.nstates);     % Matrix for linear transitions in the model (i.e. except infection)
m2 = zeros(i.nstates);    % slum Linear
m14 = zeros(i.nstates);   % Linear components of HIV tx and prevention
m3 = zeros(i.nstates);   % Linear components HIV cascade
m4 = zeros(i.nstates);   % force of transition for HHC case finding U
m5 = zeros(i.nstates);   % force of transition for HHC case finding Lf
m6 = zeros(i.nstates);   % force of transition for HHC case finding Ls
m7 = zeros(i.nstates);   % force of transition for HHC case finding I
m8 = zeros(i.nstates);   % matrix for counting screening on plhiv
m9 = zeros(i.nstates);   % matrix for counting screening on hhcU
m10 = zeros(i.nstates);   % matrix for counting screening on hhcLf
m11 = zeros(i.nstates);   % matrix for counting screening on hhcLs
m12 = zeros(i.nstates);   % matrix for counting screening on hhcI
m13 = zeros(i.nstates);   % matrix for counting screening on slum


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
            rate = p.case_finding_U(ia,ir,ig)*p.hiv_cov_tpt;
            m14(destin, source) = m14(destin, source) + rate;

            % rate = p.hiv_cov*ishivart;
            % m8(destin, source) = m8(destin, source) + rate;


            % HHC Case finding
            source = U; destin = Pu;
            rate = p.case_findingHHC_U(ia,ir,ig) * p.cfy*p.hhc_a_cov_tpt(ia);
            m4(destin, source) = m4(destin, source) + rate;

            rate = p.cfy*p.hhc_a_cov_screen(ia)*ishneg*isnoslum;
            m9(destin, source) = m9(destin, source) + rate;

            % Case finding (slums)
            source = U; destin = Pu;
            rate = p.case_findingslum_U(ia,ir,ig)*p.slum_cov_tpt*p.slum_enhance;
            m2(destin, source) = m2(destin, source) + rate;

            rate = p.slum_cov_screen*islum;
            m13(destin, source) = m13(destin, source) + rate;



            % Regimen completion
            r.outReg        =  12/(6*(1-p.tpt_short) + 3*p.tpt_short);
            p.tpt_complete(1) = (p.tpt_base_comp*(1-p.tpt_short) + p.tpt_short/2*p.tpt_target_comp + p.tpt_short/2*p.tpt_target_compa0);
            p.tpt_complete(2) = (p.tpt_base_comp*(1-p.tpt_short) + p.tpt_short*p.tpt_target_comp); 
            p.tpt_complete(3) = (p.tpt_base_comp*(1-p.tpt_short) + p.tpt_short*p.tpt_target_comp); 
            p.tpt_complete(4) = (p.tpt_base_comp*(1-p.tpt_short) + p.tpt_short*p.tpt_target_compa15p); 

            r.tpt_default  = r.outReg.*(1-p.tpt_complete)./p.tpt_complete;
            source = Pu; destin = Qu; rate = r.outReg;
            m(destin, source) = m(destin, source) + rate;

             % Regimen default
             source = Pu; destin = Qu; rate = r.tpt_default(ia);
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
                I    = gi('I');
                Dx    = gi('Dx');
                Tx    = gi('Tx');
                Tx2   = gi('Tx2');
                Rhi   = gi('Rhi');
                R     = gi('R');


                ismdr   = strcmp(strain, 'mdr');

                % Case finding (PLHIV)
                source = Lf; destin = Pf;
                rate = p.case_finding_Lf(ia,ir,ig)*p.hiv_cov_tpt;
                m14(destin, source) = m14(destin, source) + rate;

                % rate = p.hiv_cov*ishivart;
                % m8(destin, source) = m8(destin, source) + rate;

                source = Ls; destin = Ps;
                rate = p.case_finding_Ls(ia,ir,ig)*p.hiv_cov_tpt;
                m14(destin, source) = m14(destin, source) + rate;

                % rate = p.hiv_cov*ishivart;
                % m8(destin, source) = m8(destin, source) + rate;


                % HHC Case finding
                source = Lf; destin = Pf;
                rate = p.case_findingHHC_Lf(ia,ir,ig) * p.cfy*p.hhc_a_cov_tpt(ia);
                m5(destin, source) = m5(destin, source) + rate;

                rate = p.cfy*p.hhc_a_cov_screen(ia)*ishneg*isnoslum;
                m10(destin, source) = m10(destin, source) + rate;

                source = Ls; destin = Ps;
                rate = p.case_findingHHC_Ls(ia,ir,ig) * p.cfy*p.hhc_a_cov_tpt(ia) ;
                m6(destin, source) = m6(destin, source) + rate;

                rate =  p.cfy*p.hhc_a_cov_screen(ia)*ishneg*isnoslum;
                m11(destin, source) = m11(destin, source) + rate;

                % Slum Case finding
                source = Lf; destin = Pf;
                rate = p.case_findingslum_Lf(ia,ir,ig)*p.slum_cov_tpt*p.slum_enhance;
                m2(destin, source) = m2(destin, source) + rate;

                rate = p.slum_cov_screen*islum;
                m13(destin, source) = m13(destin, source) + rate;

                source = Ls; destin = Ps;
                rate = p.case_findingslum_Ls(ia,ir,ig)*p.slum_cov_tpt*p.slum_enhance;
                m2(destin, source) = m2(destin, source) + rate;

                rate = p.slum_cov_screen*islum;
                m13(destin, source) = m13(destin, source) + rate;

                 % Regimen default
                 source = Pf; destin = Qf;
                 rate =  r.tpt_default(ia);
                m(destin, source) = m(destin, source) + rate;

                source = Ps; destin = Qs;
                rate = r.tpt_default(ia);
                m(destin, source) = m(destin, source) + rate;



                % Regimen completion
                if ismdr~=1
                    source = Pf; destin = Qs; rate = r.outReg*(1-p.potency);
                    m(destin, source) = m(destin, source) + rate;

                    source = Ps; destin = Qs; rate = r.outReg*(1-p.potency);
                    m(destin, source) = m(destin, source) + rate;

                    source = Pf; destin = Qu; rate = r.outReg*(p.potency);
                    m(destin, source) = m(destin, source) + rate;

                    source = Ps; destin = Qu; rate = r.outReg*(p.potency);
                    m(destin, source) = m(destin, source) + rate;

                else
                    source = Pf; destin = Qs; rate = r.outReg*(1-p.potency)*p.rif_free;
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

                source  = Qf;
                destins = [I,                 Qs];
                rates   = [r.progression(ir,ig)*r.RRprog_age(ia), r.LTBI_stabil(ir,ig)];
                m(destins, source) = m(destins, source) + rates';

                % --- Reactivation
                source = Qs; destin = I; rate = r.reactivation(ir,ig);
                m(destin, source) = m(destin, source) + rate;


                % --- Reactivation with PT effect
                source = Pf; destin = I;
                rate = r.progression(ir,ig)*r.RRprog_age(ia)*ismdr*(1-p.pteffi*p.effondr) +...
                    r.progression(ir,ig)*r.RRprog_age(ia)*(1-ismdr)*(1-p.pteffi);
                m(destin, source) = m(destin, source) + rate;

                source = Ps; destin = I;
                rate = r.reactivation(ir,ig)*ismdr*(1-p.pteffi*p.effondr) +...
                    r.reactivation(ir,ig)*(1-ismdr)*(1-p.pteffi);
                m(destin, source) = m(destin, source) + rate;


                % Careseeking 

                source = I; destin = Dx;
                rate= r.careseeking;
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


                % FOT
                rate=r.Dx*(p_firstline+p_secondline);
                sw=p.cfy>0;
                fot(1,I)= rate*sw;



                % Interventionn Case finding (PLHIV)
                source = I;
                destin = [Tx  Tx2];
                rates = [p.case_finding_I(ia,ir,ig)*(1-ismdr)*p.hiv_cov_screen ,...
                    p.case_finding_I(ia,ir,ig)*(ismdr)*p.hiv_cov_screen];
                m14(destin, source) = m14(destin, source) + rates';

                % rates = [p.hiv_cov*(1-ismdr)*ishivart ,...
                %     p.hiv_cov*(ismdr)*ishivart];
                % m8(destin, source) = m8(destin, source) + rates';



                %Interventionn Case finding HHC Case finding
                source = I;
                destin = [Tx  Tx2];
                rates = [p.case_findingHHC_I(ia,ir,ig)*(1-ismdr)* p.cfy*p.hhc_a_cov_screen(ia) ,...
                    p.case_findingHHC_I(ia,ir,ig)*(ismdr)* p.cfy*p.hhc_a_cov_screen(ia)];
                m7(destin, source) = m7(destin, source) + rates';

                rates = [(1-ismdr)* p.cfy*p.hhc_a_cov_screen(ia)*ishneg*isnoslum ,...
                    (ismdr)* p.cfy*p.hhc_a_cov_screen(ia)*ishneg*isnoslum];
                m12(destin, source) = m12(destin, source) + rates';

                source = Dx;
                destin = [Tx  Tx2];
                rates = [p.case_findingHHC_I(ia,ir,ig)*(1-ismdr)* p.cfy*p.hhc_a_cov_screen(ia) ,...
                    p.case_findingHHC_I(ia,ir,ig)*(ismdr)* p.cfy*p.hhc_a_cov_screen(ia)];
                m7(destin, source) = m7(destin, source) + rates';

                rates = [(1-ismdr)* p.cfy*p.hhc_a_cov_screen(ia)*ishneg*isnoslum ,...
                    (ismdr)* p.cfy*p.hhc_a_cov_screen(ia)*ishneg*isnoslum];
                m12(destin, source) = m12(destin, source) + rates';

                % Interventionn Case finding (slum)
                source = I;
                destin = [Tx  Tx2];
                rates = [p.case_findingslum_I(ia,ir,ig)*(1-ismdr)*p.slum_cov_screen ,...
                    p.case_findingslum_I(ia,ir,ig)*(ismdr)*p.slum_cov_screen];
                m2(destin, source) = m2(destin, source) + rates';

                rates = [(1-ismdr)*p.slum_cov_screen*islum ,...
                    (ismdr)*p.slum_cov_screen*islum];
                m13(destin, source) = m13(destin, source) + rates';


                source = Dx;
                destin = [Tx  Tx2];
                rates = [p.case_findingslum_I(ia,ir,ig)*(1-ismdr)*p.slum_cov_screen ,...
                    p.case_findingslum_I(ia,ir,ig)*(ismdr)*p.slum_cov_screen];
                m2(destin, source) = m2(destin, source) + rates';

                rates = [(1-ismdr)*p.slum_cov_screen*islum ,...
                    (ismdr)*p.slum_cov_screen*islum];
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
                destin  = I;
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
sources = intersect(s.a0_4,s.hneg);
destins = intersect(s.a0_4,s.hpos);
rates   = p.hivcq_ad(1);
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
m3(inds) = m3(inds) + rates;

sources = intersect(s.a5_9,s.hneg);
destins = intersect(s.a5_9,s.hpos);
rates   = p.hivcq_ad(2);
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
m3(inds) = m3(inds) + rates;

sources = intersect(s.a10_15,s.hneg);
destins = intersect(s.a10_15,s.hpos);
rates   = p.hivcq_ad(3);
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
m3(inds) = m3(inds) + rates;

sources = intersect(s.a15p,s.hneg);
destins = intersect(s.a15p,s.hpos);
rates   = p.hivcq_ad(4);
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
m3(inds) = m3(inds) + rates;



%% Interventions among PLHIV

for ia = 1:length(gps.age)
    age = gps.age{ia};

    % U tp U
    sources = intersect(intersect(s.hpos,s.U),s.(age));
    destins = intersect(intersect(s.hart,s.U),s.(age));
    rates   =  r.ART_init *...
        (1 - (p.IPThiv + (1-p.IPThiv) * p.hiv_cov_tpt * p.case_findingPLHIV_U(ia)));
    inds    = sub2ind([i.nstates, i.nstates],destins,sources);
    m(inds) = m(inds) + rates;

    % U to Pu
    sources = intersect(intersect(s.hpos,s.U),s.(age));
    destins = intersect(intersect(s.hart,s.Pu),s.(age));
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
    sources = intersect(intersect(s.hpos,s.Lf),s.(age));
    destins = intersect(intersect(s.hart,s.Lf),s.(age));
    rates   =  r.ART_init *...
        (1 - (p.IPThiv + (1-p.IPThiv) * p.hiv_cov_tpt * p.case_findingPLHIV_Lf(ia)));
    inds    = sub2ind([i.nstates, i.nstates],destins,sources);
    m(inds) = m(inds) + rates;

    % Lf to P
    sources = intersect(intersect(s.hpos,s.Lf),s.(age));
    destins = intersect(intersect(s.hart,s.Pf),s.(age));
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
    sources = intersect(intersect(s.hpos,s.Ls),s.(age));
    destins = intersect(intersect(s.hart,s.Ls),s.(age));
    rates   =  r.ART_init *...
        (1 - (p.IPThiv + (1-p.IPThiv) * p.hiv_cov_tpt * p.case_findingPLHIV_Ls(ia)));
    inds    = sub2ind([i.nstates, i.nstates],destins,sources);
    m(inds) = m(inds) + rates;

    % Ls to P
    sources = intersect(intersect(s.hpos,s.Ls),s.(age));
    destins = intersect(intersect(s.hart,s.Ps),s.(age));
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
    sources = intersect(intersect(s.hpos,s.I),s.(age));
    destins = intersect(intersect(s.hart,s.I),s.(age));
    rates   = 1- p.case_findingPLHIV_I(ia) * p.hiv_cov_screen;
    inds    = sub2ind([i.nstates, i.nstates],destins,sources);
    m(inds) = m(inds) + rates;


     % I to I
    sources = intersect(intersect(s.hpos,s.I),s.(age));
    destins = intersect(intersect(s.hart,s.I),s.(age));
    rates   = 1- p.case_findingPLHIV_I(ia) * p.hiv_cov_screen;
    inds    = sub2ind([i.nstates, i.nstates],destins,sources);
    m(inds) = m(inds) + rates;


    % Dx to Dx
    sources = intersect(intersect(s.hpos,s.Dx),s.(age));
    destins = intersect(intersect(s.hart,s.Dx),s.(age));
    rates   = 1- p.case_findingPLHIV_I(ia) * p.hiv_cov_screen;
    inds    = sub2ind([i.nstates, i.nstates],destins,sources);
    m(inds) = m(inds) + rates;

    % I to Tx
    sources = intersect(intersect(intersect(s.hpos,s.I),s.(age)),s.ds);
    destins = intersect(intersect(intersect(s.hart,s.Tx),s.(age)),s.ds);
    rates   = p.case_findingPLHIV_I(ia) * p.hiv_cov_screen;
    inds    = sub2ind([i.nstates, i.nstates],destins,sources);
    m(inds) = m(inds) + rates;

    % Dx to Tx
    sources = intersect(intersect(intersect(s.hpos,s.Dx),s.(age)),s.ds);
    destins = intersect(intersect(intersect(s.hart,s.Tx),s.(age)),s.ds);
    rates   = p.case_findingPLHIV_I(ia) * p.hiv_cov_screen;
    inds    = sub2ind([i.nstates, i.nstates],destins,sources);
    m(inds) = m(inds) + rates;



    % count screen
    rates   = p.hiv_cov_screen;
    inds    = sub2ind([i.nstates, i.nstates],destins,sources);
    m8(inds) = m8(inds) + rates;

    % I to Tx2
    sources = intersect(intersect(intersect(s.hpos,s.I),s.(age)),s.mdr);
    destins = intersect(intersect(intersect(s.hart,s.Tx2),s.(age)),s.mdr);
    rates   = p.case_findingPLHIV_I(ia) * p.hiv_cov_screen;
    inds    = sub2ind([i.nstates, i.nstates],destins,sources);
    m(inds) = m(inds) + rates;

    % Dx to Tx2
    sources = intersect(intersect(intersect(s.hpos,s.Dx),s.(age)),s.mdr);
    destins = intersect(intersect(intersect(s.hart,s.Tx2),s.(age)),s.mdr);
    rates   = p.case_findingPLHIV_I(ia) * p.hiv_cov_screen;
    inds    = sub2ind([i.nstates, i.nstates],destins,sources);
    m(inds) = m(inds) + rates;


    % count screen
    rates   = p.hiv_cov_screen;
    inds    = sub2ind([i.nstates, i.nstates],destins,sources);
    m8(inds) = m8(inds) + rates;

end

%------- slum/no_slum turnover (applies only for Brazil prison)

sources = intersect(s.a15p,s.no_slum);
destins = intersect(s.a15p,s.slum);
rates   = p.entry_highrisk;
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
m(inds) = m(inds) + rates;

sources = intersect(s.a15p,s.slum);
destins = intersect(s.a15p,s.no_slum);
rates   = p.out_highrisk;
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
m(inds) = m(inds) + rates;





% --- Bring them together
M.lin    = sparse(m - diag(sum(m,1)));
M.slumlin  = sparse(m2 - diag(sum(m2,1)));
M.linHIV   = sparse(m14 - diag(sum(m14,1)));
% M.nlinHIV = sparse(m3 - diag(sum(m3,1)));
% M.fot = sparse(fot);
M.cfy_U = sparse(m4 - diag(sum(m4,1)));
M.cfy_Lf = sparse(m5 - diag(sum(m5,1)));
M.cfy_Ls = sparse(m6 - diag(sum(m6,1)));
M.cfy_I = sparse(m7 - diag(sum(m7,1)));
M.screen_hiv = sparse(m8 - diag(sum(m8,1)));
M.screen_U = sparse(m9 - diag(sum(m9,1)));
M.screen_Lf = sparse(m10 - diag(sum(m10,1)));
M.screen_Ls = sparse(m11 - diag(sum(m11,1)));
M.screen_I = sparse(m12 - diag(sum(m12,1)));
M.screen_slum = sparse(m13 - diag(sum(m13,1)));

% M.lin    = (m - diag(sum(m,1)));
% M.slumlin  = (m2 - diag(sum(m2,1)));
M.nlinHIV = sparse(m3 - diag(sum(m3,1)));
M.fot = fot;
% M.cfy_U = m4 - diag(sum(m4,1));
% M.cfy_Lf = m5 - diag(sum(m5,1));
% M.cfy_Ls = m6 - diag(sum(m6,1));
% M.cfy_I = m7 - diag(sum(m7,1));
% M.screen_hiv = (m8 - diag(sum(m8,1)));
% M.screen_U = (m9 - diag(sum(m9,1)));
% M.screen_Lf = (m10 - diag(sum(m10,1)));
% M.screen_Ls = (m11 - diag(sum(m11,1)));
% M.screen_I = (m12 - diag(sum(m12,1)));
% M.screen_slum = (m13 - diag(sum(m13,1)));

 % spr=@(M) 1e2*numel(nonzeros(M))/prod(size(M));
% 
% spr(M.lin   )% (m - diag(sum(m,1)));
% spr(M.slumlin )% (m2 - diag(sum(m2,1)));
% spr(M.nlinHIV)% (m3 - diag(sum(m3,1)));
% spr(M.fot)% fot;
% spr(M.cfy_U)% m4 - diag(sum(m4,1));
% spr(M.cfy_Lf)% m5 - diag(sum(m5,1));
% spr(M.cfy_Ls)% m6 - diag(sum(m6,1));
% spr(M.cfy_I)% m7 - diag(sum(m7,1));
% spr(M.screen_hiv)% (m8 - diag(sum(m8,1)));
% spr(M.screen_U)% (m9 - diag(sum(m9,1)));
% spr(M.screen_Lf)% (m10 - diag(sum(m10,1)));
% spr(M.screen_Ls)% (m11 - diag(sum(m11,1)));
% spr(M.screen_I)% (m12 - diag(sum(m12,1)));
% spr(M.screen_slum)% (m13 - diag(sum(m13,1)));

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

                m(Lfde, [U Ls Lfso Rhi R]) = 1;
                m(Qfde, [Qu Qfso Qs]) = 1;

                % Adjust for imunuty
                cols = [Lfso, Ls, Qfso, Qs, Rhi, R];
                m(:,cols) = m(:,cols)*(1-p.imm(ig));
            end

        end
        m(:,[s.hpos s.hart]) = m(:,[s.hpos s.hart])*p.HIVlam;
        m(:,s.slum) = m(:,s.slum)*r.RRslum;
          M.nlin.(age).(strain) = sparse(m - diag(sum(m,1)));
         % M.nlin.(age).(strain) = (m - diag(sum(m,1)));
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

% Third column is for HIV+ve TB mortality
inds = intersect(s.hpos, s.infectious);
m(inds,3) = r.mort_TB(2);
inds = intersect(s.hart, s.infectious);
m(inds,3) = r.mort_TB(3);
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

% Third column is for HIV+ve TB mortality
inds = intersect(s.hpos, s.infectious);
m(inds) = r.mort_TB(2);
inds = intersect(s.hart, s.infectious);
m(inds) = r.mort_TB(3);

M.mu = m';
