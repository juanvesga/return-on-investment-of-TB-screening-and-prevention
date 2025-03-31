
%% Map population and model structure
gps.age    = population_table.Properties.VariableNames(1:cols-1);%{'a0_4','a5_9','a10_15','a15p'};
gps.strain = {'ds' ,'mdr'};
gps.georisk   = {'no_slum' ,'slum'};
gps.hiv    = {'hneg' ,'hpos','hart'};

states0 = {...
    'U',...  % TB uninfected
    'Pu',...
    'Qu'...
    };                                          % Structured by geography alone


states1 = {...
    'Lf',...% Latent fast progressors
    'Ls',...% Latent slow progressors
    'Pf',...% TPT for fast progs
    'Ps',...% TPT for fast progs
    'Qf',...% post TPT
    'Qs',...% post TPT
    'I',...  % Active TB symptomatic
    'Dx',... % Seeking TB care 
    'Tx',... % FL treatment
    'Tx2',...% SL treatment (MDR)
    'Rhi',...% Non TX completion or wrong treatment
    'R'};    % Recovered stable

% Create multidimension index structure
[i, s, d, lim] = get_addresses({states0, gps.age, gps.georisk, gps.hiv}, [], [], [], 0);
[i, s, d, lim] = get_addresses({states1, gps.age, gps.georisk, gps.hiv,gps.strain}, i, s, d, lim);
d = char(d);

%Include the auxiliaries
auxnames = {...
    'inc',...
    'notif',...
    'mort',...
    'pt',...
    'newart',...
    'hiv', ...
    'inc_ds',...
    'inc_dr',...
    'inc_plhiv',...
    'inc_slum',...
    'screen_all',...
    'screen_plhiv',...
    'screen_hhc',...
    'screen_slum',...
    'tst_all',...
    'tst_plhiv',...
    'tst_hhc',...
    'tst_slum',...
    'tpt_all',...
    'tpt_plhiv',...
    'tpt_hhc',...
    'tpt_slum',...
    'tx_fl_all',...
    'tx_fl_plhiv',...
    'tx_fl_hhc',...
    'tx_fl_slum',...
    'tx_sl_all',...
    'tx_sl_plhiv',...
    'tx_sl_hhc',...
    'tx_sl_slum',...
    'daly_all',...
    'daly_plhiv',...
    'daly_slum',...
    'yld_all',...
    'yld_plhiv',...
    'yld_slum',...
     'yld_ptb',...
    'yll_all',...
    'yll_plhiv',...
    'yll_slum',...
    'yll_ptb',...
    'mu_ds',...
    'mu_dr',...
    'mu_ptb',...
    'lam_ds',...
    'lam_dr',...
    'screen_hhc_source'...
    };

auxinds  = [...
    4 ,... inc
    2,...notif
    2, ... mort
    3,...pt
    2,... newart
    4, ...hiv
    3*numel(gps.age),... inc_ds
    3*numel(gps.age),... inc_dr
    1,... inc_plhiv
    1,... inc_slum
    numel(gps.age),... screen
    numel(gps.age),... screen
    numel(gps.age),... screen
    numel(gps.age),... screen
    numel(gps.age),... tst
    numel(gps.age),... tst
    numel(gps.age),... tst
    numel(gps.age),... tst
    numel(gps.age),... tpt
    numel(gps.age),... tpt
    numel(gps.age),... tpt
    numel(gps.age),... tpt
    numel(gps.age),... tx_fl
    numel(gps.age),... tx_fl
    numel(gps.age),... tx_fl
    numel(gps.age),... tx_fl
    numel(gps.age),... tx_sl
    numel(gps.age),... tx_sl
    numel(gps.age),... tx_sl
    numel(gps.age),... tx_sl
    numel(gps.age),... daly
    numel(gps.age),... daly
    numel(gps.age),... daly
    numel(gps.age),... yld
    numel(gps.age),... yld
    numel(gps.age),... yld
    numel(gps.age),... yld
    numel(gps.age),... yll
    numel(gps.age),... yll
    numel(gps.age),... yll
    numel(gps.age),... yll
    3*numel(gps.age),... mu_ds
    3*numel(gps.age),... mu_dr
     numel(gps.age),... mu_ptb
     numel(gps.age),...%lam_ds
     numel(gps.age),...%lam_ds
     4]; % screen HHCsource 



for ii = 1:length(auxnames)
    inds = lim + (1:auxinds(ii));
    i.aux.(auxnames{ii}) = inds;
    lim = inds(end);
end
i.nx = lim;
s.nstates=(1:i.nstates);


% General states
s.ch = [s.a0_4, s.a5_9, s.a10_15];
s.ad = s.a15p;
s.hivpositive= [s.hpos, s.hart];
if (strcmp(country,"GEO")||strcmp(country,"BRA"))
s.hhc        = intersect(s.nstates, s.no_slum);%
else
s.hhc        = s.nstates;
end

s.infectious = [s.I,s.Dx,intersect(s.Tx,s.mdr)];
s.prevalent  = [s.I, s.Dx];
s.ltbiprevalent  = [s.Lf s.Ls s.Pf s.Ps s.Qf s.Qs];
s.TBmort     = [s.I s.Dx s.Tx s.Tx2];
s.ipt        = [s.Pu s.Pf s.Ps];
s.treated    = [s.Tx s.Tx2];
s.notif      = [s.Dx];
s.tbcare     = [s.Tx s.Tx2];
s.notbcare   = [s.Pu, s.Qu, s.Pf, s.Qf, s.Qs, s.Ps ,s.I, s.Rhi  s.R];
s.postTB     = [s.Rhi s.R];
s.postTByll  = [s.Rhi];
%% Matrix selectors and aggregators for model transitions

% block unallowed transitions
tmp = zeros(i.nstates);
tmp(:,:) = 1;
for id = 1:numel(gps.age)
    for is = 1:numel(gps.age)
        dest= s.(gps.age{id});
        sour= s.(gps.age{is});
        if id==is; continue;
        else
            tmp(dest,sour)=0;
        end
    end
end
for id = 1:numel(gps.hiv)
    for is = 1:numel(gps.hiv)
        dest= s.(gps.hiv{id});
        sour= s.(gps.hiv{is});
        if id==is; continue;
        else
            tmp(dest,sour)=0;
        end
    end
end
for id = 1:numel(gps.georisk)
    for is = 1:numel(gps.georisk)
        dest= s.(gps.georisk{id});
        sour= s.(gps.georisk{is});
        if id==is; continue;
        else
            tmp(dest,sour)=0;
        end
    end
end
tmp(s.treated,[s.U,s.Lf,s.Ls]) = 0;
tmp(s.ipt,s.I) = 0;
check=tmp;

%% Selectors for fitting

% inc
tmp = zeros(4,i.nstates);
tmp(1,s.I) = 1;                           %All
tmp(2,intersect(s.I,s.mdr)) = 1;          %MDR
tmp(3,intersect(s.I,s.hivpositive)) = 1;%HIV/TB
tmp(4,intersect(s.I,s.slum)) = 1;%HIV/TB
agg.inc = sparse(tmp);

tmp = zeros(i.nstates);
tmp(s.I,:) = 1;
tmp=tmp.*check;
sel.inc = sparse(tmp - diag(diag(tmp)));

tmp = zeros(i.nstates);
tmp(intersect(s.Tx,s.mdr),intersect(s.Tx,s.ds)) = 1;
tmp=tmp.*check;
sel.acqu = sparse(tmp - diag(diag(tmp)));

% TB notification
tmp = zeros(2,i.nstates);
tmp(1,s.notif)  = 1;
tmp(2,intersect(s.notif,s.mdr))  = 1;
agg.notif = sparse(tmp);

tmp = zeros(i.nstates);
tmp(s.notif,:) = 1;
tmp=tmp.*check;
sel.notif = sparse(tmp - diag(diag(tmp)));

%PT
tmp = zeros(3,i.nstates);
tmp(1,s.ipt)  = 1;
tmp(2,intersect(s.ipt,s.hart))  = 1;% PT
tmp(3,intersect(s.ipt,[s.hneg s.hpos]))  = 1;% PT a
agg.ipt = sparse(tmp);

tmp = zeros(i.nstates);
tmp(s.Pu,s.U) = 1;
tmp(s.Pf,s.Lf) = 1;
tmp(s.Ps,s.Ls) = 1;
tmp=tmp.*check;
for ii = 1:(numel(gps.age))
    age=gps.age{ii};
    tmp(intersect(intersect(s.Pu,s.hart),s.(age)),intersect(intersect(s.hpos,s.U ),s.(age))) = 1;
    tmp(intersect(intersect(s.Pf,s.hart),s.(age)),intersect(intersect(s.hpos,s.Lf ),s.(age))) = 1;
    tmp(intersect(intersect(s.Ps,s.hart),s.(age)),intersect(intersect(s.hpos,s.Ls ),s.(age))) = 1;
end
tmp(s.ds,s.mdr)=0;
tmp(s.mdr,s.ds)=0;
tmp(s.slum,s.no_slum)=0;
tmp(s.no_slum,s.slum)=0;
sel.ipt = sparse(tmp - diag(diag(tmp)));


tmp = zeros(i.nstates);

tmp=tmp.*check;

 tmp(intersect(s.Pu,s.hart),intersect(s.hpos,s.U )) = 1;
 tmp(intersect(s.Pf,s.hart),intersect(s.hpos,s.Lf )) = 1;
  tmp(intersect(s.Ps,s.hart),intersect(s.hpos,s.Ls )) = 1;
tmp(s.ds,s.mdr)=0;
tmp(s.mdr,s.ds)=0;
tmp(s.slum,s.no_slum)=0;
tmp(s.no_slum,s.slum)=0;
sel.newartipt = sparse(tmp - diag(diag(tmp)));


%% IPT HHC
tmp = zeros(i.nstates);
tmp(intersect(s.Pu,s.hhc),intersect(s.U,s.hhc)) = 1;
tmp(intersect(s.Pf,s.hhc),intersect(s.Lf,s.hhc)) = 1;
tmp(intersect(s.Ps,s.hhc),intersect(s.Ls,s.hhc)) = 1;
tmp=tmp.*check;

tmp(s.ds,s.mdr)=0;
tmp(s.mdr,s.ds)=0;
tmp(s.slum,s.no_slum)=0;
tmp(s.no_slum,s.slum)=0;
sel.ipt_hhc = sparse(tmp - diag(diag(tmp)));

%% IPT Slum
tmp = zeros(i.nstates);
tmp(intersect(s.Pu,s.slum),intersect(s.U,s.slum)) = 1;
tmp(intersect(s.Pf,s.slum),intersect(s.Lf,s.slum)) = 1;
tmp(intersect(s.Ps,s.slum),intersect(s.Ls,s.slum)) = 1;
tmp=tmp.*check;

tmp(s.ds,s.mdr)=0;
tmp(s.mdr,s.ds)=0;
tmp(s.slum,s.no_slum)=0;
tmp(s.no_slum,s.slum)=0;
sel.ipt_slum = sparse(tmp - diag(diag(tmp)));


% Selectors for New on ART
tmp = zeros(2,i.nstates);
tmp(1,[intersect(s.hart,[s.U s.Ls s.Lf]) intersect(s.hart,s.ipt)]  ) = 1;
tmp(2,s.hart)  = 1;
agg.art = sparse(tmp);

tmp = zeros(i.nstates);
sources = s.hpos;
destins = s.hart;
rates   = 1;
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
tmp(inds) = tmp(inds) + rates;
for ii = 1:(numel(gps.age))
    age=gps.age{ii};
    tmp(intersect(intersect(s.Pu,s.hart),s.(age)),intersect(intersect(s.hpos,s.U ),s.(age))) = 1;
    tmp(intersect(intersect(s.Pf,s.hart),s.(age)),intersect(intersect(s.hpos,s.Lf ),s.(age))) = 1;
    tmp(intersect(intersect(s.Ps,s.hart),s.(age)),intersect(intersect(s.hpos,s.Ls ),s.(age))) = 1;
end
tmp(s.ds,s.mdr)=0;
tmp(s.mdr,s.ds)=0;
tmp(s.slum,s.no_slum)=0;
tmp(s.no_slum,s.slum)=0;
sel.art = sparse(tmp - diag(diag(tmp)));


% Selectors for  HIV Incidence
tmp = zeros(4,i.nstates);
tmp(1,s.hpos)  = 1;
tmp(2,intersect(s.hpos,s.a0_4 ) )  = 1;
tmp(3,intersect(s.hpos,s.a5_9 ) )  = 1;
tmp(4,intersect(s.hpos,[s.a10_15 s.a15p] ) )  = 1;
agg.hiv = sparse(tmp);

tmp = zeros(i.nstates);
sources = s.hneg;
destins = s.hpos;
rates   = 1;
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
tmp(inds) = tmp(inds) + rates;
sel.hiv = sparse(tmp - diag(diag(tmp)));


%% Selectors for model output interventions
groups={'nstates','hivpositive','slum'};
groups2={'nstates','hivpositive','slum','hhc'};

% inc_ds
tmp = zeros(3*numel(gps.age),i.nstates);
for jj = 0:2
    group=s.(groups{jj+1});
    for ii = 0:(numel(gps.age)-1)
        age=gps.age{ii+1};
        index= 1+(jj * numel(gps.age) + ii);

        tmp(index,intersect(...
            intersect(s.I,s.ds),intersect(s.(age),group)...
            )) = 1;
    end
end
agg.inc_ds = sparse(tmp);

tmp = zeros(i.nstates);
tmp(intersect(s.I,s.ds),s.ds) = 1;
tmp=tmp.*check;
sel.inc_ds = sparse(tmp - diag(diag(tmp)));

% inc_dr
tmp = zeros(3*numel(gps.age),i.nstates);
for jj = 0:2
    group=s.(groups{jj+1});
    for ii = 0:(numel(gps.age)-1)
        age=gps.age{ii+1};
        index= 1+(jj * numel(gps.age) + ii);

        tmp(index,intersect(...
            intersect(s.I,s.mdr),intersect(s.(age),group)...
            )) = 1;
    end
end
agg.inc_dr = sparse(tmp);

tmp = zeros(i.nstates);
tmp(intersect(s.I,s.mdr),s.mdr) = 1;
tmp=tmp.*check;
sel.inc_dr = sparse(tmp - diag(diag(tmp)));



% Inc PLHIV
tmp = zeros(1,i.nstates);
tmp(1,intersect(s.I,s.hivpositive))  = 1;
agg.inc_plhiv = sparse(tmp);

tmp = zeros(i.nstates);
tmp(intersect(s.I,s.hivpositive),s.hivpositive) = 1;
tmp=tmp.*check;
sel.inc_plhiv = sparse(tmp - diag(diag(tmp)));


% Inc Slums
tmp = zeros(1,i.nstates);
tmp(1,intersect(s.I,s.slum))  = 1;
agg.inc_slum = sparse(tmp);

tmp = zeros(i.nstates);
tmp(intersect(s.I,s.slum),s.slum) = 1;
tmp=tmp.*check;
sel.inc_slum = sparse(tmp - diag(diag(tmp)));




%% Aggregators for screen, tst, tpt, treatment counts
[tmp0, tmp1, tmp2, tmp3,...
    tmp4, tmp5, tmp6, tmp7,...
    tmp8 ,tmp9 ,tmp10 ,tmp11,...
    tmp12, tmp13, tmp14, tmp15,...
    tmp16, tmp17, tmp18, tmp19] = deal(zeros(numel(gps.age),i.nstates));

for ii = 1:(numel(gps.age))
    age=gps.age{ii};
    %screen
    tmp0(ii,intersect([s.treated,s.ipt],intersect(s.(age),s.nstates))) = 1;
    tmp1(ii,intersect([s.treated,s.ipt],intersect(s.(age),s.hivpositive))) = 1;
    tmp2(ii,intersect([s.treated,s.ipt],intersect(s.(age),s.hhc))) = 1;
    tmp3(ii,intersect([s.treated,s.ipt],intersect(s.(age),s.slum))) = 1;
    % TST
    val=1;
    if strcmp(age,'a0_4')
        val=0;
    end
    tmp4(ii,intersect([s.ipt],intersect(s.(age),s.nstates))) = val;
    tmp5(ii,intersect([s.ipt],intersect(s.(age),s.hivpositive))) = 0;
    tmp6(ii,intersect([s.ipt],intersect(s.(age),s.hhc))) = val;
    tmp7(ii,intersect([s.ipt],intersect(s.(age),s.slum))) = val;
    % TPT
    tmp8(ii,intersect([s.ipt],intersect(s.(age),s.nstates))) = 1;
    tmp9(ii,intersect([s.ipt],intersect(s.(age),s.hivpositive))) = 1;
    tmp10(ii,intersect([s.ipt],intersect(s.(age),s.hhc))) = 1;
    tmp11(ii,intersect([s.ipt],intersect(s.(age),s.slum))) = 1;

    % TXfl
    tmp12(ii,intersect([s.Tx],intersect(s.(age),s.nstates))) = 1;
    tmp13(ii,intersect([s.Tx],intersect(s.(age),s.hivpositive))) = 1;
    tmp14(ii,intersect([s.Tx],intersect(s.(age),s.hhc))) = 1;
    tmp15(ii,intersect([s.Tx],intersect(s.(age),s.slum))) = 1;

    % TXsl
    tmp16(ii,intersect([s.Tx2],intersect(s.(age),s.nstates))) = 1;
    tmp17(ii,intersect([s.Tx2],intersect(s.(age),s.hivpositive))) = 1;
    tmp18(ii,intersect([s.Tx2],intersect(s.(age),s.hhc))) = 1;
    tmp19(ii,intersect([s.Tx2],intersect(s.(age),s.slum))) = 1;

end

agg.screen_all = sparse(tmp0);
agg.screen_plhiv = sparse(tmp1);
agg.screen_hhc = sparse(tmp2);
agg.screen_slum = sparse(tmp3);

agg.tst_all = sparse(tmp4);
agg.tst_plhiv = sparse(tmp5);
agg.tst_hhc = sparse(tmp6);
agg.tst_slum = sparse(tmp7);

agg.tpt_all = sparse(tmp8);
agg.tpt_plhiv = sparse(tmp9);
agg.tpt_hhc = sparse(tmp10);
agg.tpt_slum = sparse(tmp11);

agg.tx_fl_all = sparse(tmp12);
agg.tx_fl_plhiv = sparse(tmp13);
agg.tx_fl_hhc = sparse(tmp14);
agg.tx_fl_slum = sparse(tmp15);

agg.tx_sl_all = sparse(tmp16);
agg.tx_sl_plhiv = sparse(tmp17);
agg.tx_sl_hhc = sparse(tmp18);
agg.tx_sl_slum = sparse(tmp19);

tmp = zeros(4,i.nstates);
tmp(1,intersect(s.Pu,s.hhc)) = 1;
tmp(2,intersect(s.Pf,s.hhc)) = 1;
tmp(3,intersect(s.Ps,s.hhc)) = 1;
tmp(4,intersect([s.Tx s.Tx2],s.hhc)) = 1;
agg.screen_hhc_source = sparse(tmp);

%% Selectors for tsts, tpt screen and treatment counts
tmp = zeros(i.nstates);
tmp(s.Pu,s.U) = 1;
tmp(s.Pf,s.Lf) = 1;
tmp(s.Ps,s.Ls) = 1;
tmp(s.treated,s.I) = 1;
tmp=tmp.*check;
for ii = 1:(numel(gps.age))
    age=gps.age{ii};
    tmp(intersect(intersect(s.Pu,s.hart),s.(age)),intersect(intersect(s.hpos,s.U ),s.(age))) = 1;
    tmp(intersect(intersect(s.Pf,s.hart),s.(age)),intersect(intersect(s.hpos,s.Lf ),s.(age))) = 1;
    tmp(intersect(intersect(s.Ps,s.hart),s.(age)),intersect(intersect(s.hpos,s.Ls ),s.(age))) = 1;
    tmp(intersect(intersect(s.treated,s.hart),s.(age)),intersect(intersect(s.hpos,[s.I] ),s.(age))) = 1;
end
tmp(s.ds,s.mdr)=0;
tmp(s.mdr,s.ds)=0;
tmp(s.slum,s.no_slum)=0;
tmp(s.no_slum,s.slum)=0;
sel.screen = sparse(tmp - diag(diag(tmp)));

%% Screen HHC
tmp = zeros(i.nstates);
tmp(intersect(s.Pu,s.hhc),intersect(s.U,s.hhc)) = 1;
tmp(intersect(s.Pf,s.hhc),intersect(s.Lf,s.hhc)) = 1;
tmp(intersect(s.Ps,s.hhc),intersect(s.Ls,s.hhc)) = 1;
tmp(intersect(s.treated,s.hhc),intersect(s.I,s.hhc)) = 1;

tmp=tmp.*check;
tmp(s.ds,s.mdr)=0;
tmp(s.mdr,s.ds)=0;
tmp(s.slum,s.no_slum)=0;
tmp(s.no_slum,s.slum)=0;
sel.screen_hhc = sparse(tmp - diag(diag(tmp)));


%% Screen slum
tmp = zeros(i.nstates);
tmp(intersect(s.Pu,s.slum),intersect(s.U,s.slum)) = 1;
tmp(intersect(s.Pf,s.slum),intersect(s.Lf,s.slum)) = 1;
tmp(intersect(s.Ps,s.slum),intersect(s.Ls,s.slum)) = 1;
tmp(intersect(s.treated,s.slum),intersect(s.I,s.slum)) = 1;

tmp=tmp.*check;
tmp(s.ds,s.mdr)=0;
tmp(s.mdr,s.ds)=0;
tmp(s.slum,s.no_slum)=0;
tmp(s.no_slum,s.slum)=0;
sel.screen_slum = sparse(tmp - diag(diag(tmp)));


tmp = zeros(i.nstates);
tmp(s.Tx,:) = 1;
tmp=tmp.*check;
for ii = 1:(numel(gps.age))
    age=gps.age{ii};
    tmp(intersect(intersect(s.Tx,s.hart),s.(age)),intersect(intersect(s.hpos,[s.I s.Dx] ),s.(age))) = 1;
end
tmp(s.ds,s.mdr)=0;
tmp(s.mdr,s.ds)=0;
tmp(s.slum,s.no_slum)=0;
tmp(s.no_slum,s.slum)=0;
sel.tx_fl = sparse(tmp - diag(diag(tmp)));

tmp = zeros(i.nstates);
tmp(s.Tx2,:) = 1;
tmp=tmp.*check;
for ii = 1:(numel(gps.age))
    age=gps.age{ii};
    tmp(intersect(intersect(s.Tx2,s.hart),s.(age)),intersect(intersect(s.hpos,[s.I s.Dx] ),s.(age))) = 1;
end
tmp(s.ds,s.mdr)=0;
tmp(s.mdr,s.ds)=0;
tmp(s.slum,s.no_slum)=0;
tmp(s.no_slum,s.slum)=0;
sel.tx_sl = sparse(tmp - diag(diag(tmp)));


%% Selectors for DALYs
tmp0 = zeros(numel(gps.age) ,i.nstates);
yldptb=tmp0;
%YLD
for ii = 1:(numel(gps.age))
    age=gps.age{ii};
tmp0(ii,intersect(s.infectious,s.(age))) = 0.333; % Weight TB disease
tmp0(ii,intersect(intersect(s.hivpositive,s.infectious),s.(age)))=0.408; % Weight TB/HIV
% ind=setdiff(s.hpos,s.infectious);
% tmp0(ii,intersect(ind,s.(age))) = 0.274; % Weight HIV w/o ART
% ind=setdiff(s.hart,s.infectious);
% tmp0(ii,intersect(ind,s.(age))) = 0.078; % Weight HIV w ART
ind=s.postTB;
tmp0(ii,intersect(ind,s.(age))) = 0.036; % Weight postTB
yldptb(ii,intersect(ind,s.(age))) = 0.036; % Weight postTB
end

tmp1=tmp0;
ind=setdiff(s.nstates,s.hivpositive);
tmp1(:,ind) = 0; % remove non-hiv
tmp1(:,s.slum) = 0; % remove slum

tmp2=tmp0;
tmp2(:,s.no_slum) = 0; % remove hiv
tmp2(:,s.hivpositive) = 0; % remove hiv

agg.yld_all = sparse(tmp0);
agg.yld_plhiv = sparse(tmp1);
agg.yld_slum = sparse(tmp2);
agg.yld_ptb = sparse(yldptb);

%YLL
%p.weights_age=[87.007 , 82.03, 76.05, 37.7];
pars  =load_data('data/params.xlsx');
weights_age=[
pars{country,'yll0_4'},...
pars{country,'yll5_9'},...
pars{country,'yll10_14'},...
pars{country,'yll15p'}];
tmp0 = zeros(numel(gps.age) ,i.nstates);
ptbyll = zeros(numel(gps.age) ,i.nstates);
ptbmu = zeros(numel(gps.age) ,i.nstates);

for ii = 1:(numel(gps.age))
age=gps.age{ii};
tmp0(ii,intersect(s.infectious,s.(age))) = weights_age(ii); % TB by age
% tmp0(ii,intersect(s.hpos,s.(age)))= weights_age(ii); % HIV 
tmp0(ii,intersect(s.postTByll,s.(age)))= weights_age(ii); %  
ptbyll(ii,intersect(s.postTByll,s.(age)))= weights_age(ii); % PTBD 
ptbmu(ii,intersect(s.postTByll,s.(age)))= 1; % PTBD 
end

% tmp0(1,intersect(s.infectious,s.a0_4))=weights_age(1);
% tmp0(2,intersect(s.infectious,s.a5_9))=weights_age(2);
% tmp0(3,intersect(s.infectious,s.a10_15))=weights_age(3);
% tmp0(4,intersect(s.infectious,s.a15p))=weights_age(4);
% 
% tmp0(1,intersect(s.hpos,s.a0_4))=weights_age(1);
% tmp0(2,intersect(s.hpos,s.a5_9))=weights_age(2);
% tmp0(3,intersect(s.hpos,s.a10_15))=weights_age(3);
% tmp0(4,intersect(s.hpos,s.a15p))=weights_age(4);
% 
% tmp0(1,intersect(s.postTByll,s.a0_4))=weights_age(1);
% tmp0(2,intersect(s.postTByll,s.a5_9))=weights_age(2);
% tmp0(3,intersect(s.postTByll,s.a10_15))=weights_age(3);
% tmp0(4,intersect(s.postTByll,s.a15p))=weights_age(4);

tmp1=tmp0;
ind=setdiff(s.nstates,s.hivpositive);
tmp1(:,ind) = 0; % remove non-hiv
tmp1(:,s.slum) = 0; % remove slum

tmp2=tmp0;
tmp2(:,s.no_slum) = 0; % remove hiv
tmp2(:,s.hivpositive) = 0; % remove hiv

agg.yll_all = sparse(tmp0);
agg.yll_plhiv = sparse(tmp1);
agg.yll_slum = sparse(tmp2);
agg.yll_ptb = sparse(ptbyll);
agg.mu_ptb = sparse(ptbmu);
% Sel for Mortality DS
tmp = zeros(3*numel(gps.age) ,i.nstates);

for jj = 0:2
    group=s.(groups{jj+1});
    for ii = 0:(numel(gps.age)-1)
        age=gps.age{ii+1};
        index= 1+(jj * numel(gps.age) + ii);

        tmp(index,intersect(...
            intersect(s.infectious,s.ds),intersect(s.(age),group)...
            )) = 1;
    end
end
agg.mu_ds = sparse(tmp);


% Sel for Mortality DR
tmp = zeros(3*numel(gps.age) ,i.nstates);

for jj = 0:2
    group=s.(groups{jj+1});
    for ii = 0:(numel(gps.age)-1)
        age=gps.age{ii+1};
        index= 1+(jj * numel(gps.age) + ii);

        tmp(index,intersect(...
            intersect(s.infectious,s.mdr),intersect(s.(age),group)...
            )) = 1;
    end
end
agg.mu_dr = sparse(tmp);


% Reference structure
ref.i = i; ref.s = s; ref.d = d;

clear auxinds auxnames check inds lim tmp i s d lim dest sour tmp0 tmp1 tmp2 tmp3 tmp4  tmp5 tmp6 tmp7 tmp8 tmp9 tmp10 tmp11 tmp12 tmp13 tmp14 tmp15 tmp16 tmp17 tmp18 tmp19