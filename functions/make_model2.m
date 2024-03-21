%% This Funtions constructs Model structure with trasnitions (linear and non-linear)

function M = make_model2(p, r, i, s, gps)

% --- Structures to hold matrices
% --- Get the linear rates ------------------------------------------------
m  = zeros(i.nstates);     % Matrix for linear transitions in the model (i.e. except infection)
m2 = zeros(i.nstates);    % Separate linear matrix for linkage to treatment - this will be multiplied by a time-dependent factor in 'goveqs_basis', to reflect COVID-related disruptions
m3 = zeros(i.nstates);   % Linear components HIV cascade


% --- Apply diagnostic transitions to NTS and Preventive therapy for all model
for ia = 1:length(gps.age)
    age = gps.age{ia};

    for ir = 1:length(gps.georisk)
        risk = gps.georisk{ir};

        for ig = 1:length(gps.hiv)
            hiv = gps.hiv{ig};
            getind1 = @(st) i.(st).(age).(risk).(hiv);

            U    = getind1('U');
            Pu   = getind1('Pu');
            Qu   = getind1('Qu');


            % Recruitment into PT
            source = U; destin = Pu;
            rate = (1-p.LTBItest)*p.IPT(ig)*p.housedist(1);
            m(destin, source) = m(destin, source) + rate;

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
            source = Pu; destin = Qu;
            rate = r.tpt_default(ia);
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

                % Recruitment into PT
                source = Lf; destin = Pf;
                rate = p.IPT(ig)*p.housedist(2);
                m(destin, source) = m(destin, source) + rate;

                source = Ls; destin = Ps;
                rate =p.IPT(ig)*p.housedist(3);
                m(destin, source) = m(destin, source) + rate;

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
                m2(destin, source) = m2(destin, source) + rate';


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
% --- Get HIV cascade linear transitions  ------ to ART
sources = intersect(s.hpos,s.U); destins = intersect(s.hart,s.U);
rates   =  r.ART_init*(1-p.IPThiv);
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
m(inds) = m(inds) + rates;

sources = intersect(s.hpos,s.U); destins = intersect(s.hart,s.Pu);
rates   =  r.ART_init*p.IPThiv;
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
m(inds) = m(inds) + rates;

sources = intersect(s.hpos,s.Lf); destins = intersect(s.hart,s.Lf);
rates   = r.ART_init*(1-p.IPThiv);
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
m(inds) = m(inds) + rates;

sources = intersect(s.hpos,s.Lf); destins = intersect(s.hart,s.Pf);
rates   = r.ART_init*p.IPThiv;
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
m(inds) = m(inds) + rates;

sources = intersect(s.hpos,s.Ls); destins = intersect(s.hart,s.Ls);
rates   = r.ART_init*(1-p.IPThiv);
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
m(inds) = m(inds) + rates;

sources = intersect(s.hpos,s.Ls); destins = intersect(s.hart,s.Ps);
rates   = r.ART_init*p.IPThiv;
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
m(inds) = m(inds) + rates;

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
M.Dxlin  = sparse(m2 - diag(sum(m2,1)));
M.nlinHIV = sparse(m3 - diag(sum(m3,1)));


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
                cols = [Lfso, Ls, Qfso, Qs,  Rhi, R];
                m(:,cols) = m(:,cols)*(1-p.imm(ig));
            end

        end
        m(:,[s.hpos s.hart]) = m(:,[s.hpos s.hart])*p.HIVlam;
        m(:,s.slum) = m(:,s.slum)*r.RRslum;
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



M.lambda_ds = (m_ds);
M.lambda_dr = (m_dr);


%% -- - Get the mortality rates
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

