function results=return_output_single(fx,x,ref,prm,location,file,type,varargin)

i = ref.i; s = ref.s;
r = prm.r; p = prm.p;
xvals=x;
nsims=1;
runtype=type;
if ~isempty(varargin)
    if (strcmp(type,'mcmc'))
        tmp=varargin{1}; burn=tmp(1);thin=tmp(2);
        xvals=x(burn:thin:end,:);
    end
    nsims=size(xvals,1);
end

ii=nsims;
x0=xvals(ii,:);
%
% profile on
aux = fx(x0);
% profile viewer
% profile off

soln=aux.solnL;
model=aux.model;

t=soln(:,end);
tL=aux.solnL(:,end);


tmp = interp1(t,sum(soln(:,1:i.nstates),2),t(1):t(end));
pop = (tmp(1:end-1)+tmp(2:end))/2;

tmp = interp1(t,sum(soln(:,s.a15p),2),t(1):t(end));
pop_ad = (tmp(1:end-1)+tmp(2:end))/2;

tmp = interp1(t,sum(soln(:,s.slum),2),t(1):t(end));
pop_slum = (tmp(1:end-1)+tmp(2:end))/2;

tmp = interp1(tL,sum(soln(:,1:i.nstates),2),tL(1):tL(end));
popL = (tmp(1:end-1)+tmp(2:end))/2;

tmp = interp1(tL,sum(soln(:,s.slum),2),tL(1):tL(end));
popslumL = (tmp(1:end-1)+tmp(2:end))/2;


% Prevalence
%Notif (All)
tmp  = diff(interp1(t,soln(:,i.aux.notif),t(1):t(end)),1);
notif=1e5*(tmp(:,1)./pop')';

%Prevalence
% tmp = interp1(t,sum(soln(:,intersect(s.prevalent,s.a15p)),2),t(1):t(end));
tmp = interp1(t,sum(soln(:,s.prevalent),2),t(1):t(end));
num = (tmp(1:end-1)+tmp(2:end))/2;
prev= 1e5*(num./pop_ad);

%Prevalence high risk
tmp = interp1(t,sum(soln(:,intersect(s.prevalent,s.slum)),2),t(1):t(end));
num = (tmp(1:end-1)+tmp(2:end))/2;
prev_hi= 1e5*(num./pop_slum);


%Incidence
%Inc (All , MDR, TB/HIV)
tmp = diff(interp1(t,soln(:,i.aux.inc),t(1):t(end)),1);
inc=1e5*(tmp./[pop' pop' pop' pop_slum'])';
inc_all= inc(1,:);
inc_mdr= inc(2,:);
inc_tbhiv= inc(3,:);
inc_slum= inc(4,:);


%  Mort TB in HIV negative    Mort TB in HIV +
tmp  = diff(interp1(t,soln(:,i.aux.mort),t(1):t(end)),1);
mort= 1e5*[tmp(:,1)./pop' ,...
    tmp(:,2)./pop']';

mort_tbhn= mort(1,:);
mort_tbhiv= mort(2,:);


%Incidence
tmp = diff(interp1(tL,soln(:,i.aux.inc),tL(1):tL(end)),1);
incL = 1e5.*([tmp(:,1) , tmp(:,2), tmp(:,3) ]./[popL; popL; popL]');
incL_all=incL(:,1);
incL_mdr=incL(:,2);
incL_tbhiv=incL(:,3);

%cases
tmp = diff(interp1(t,soln(:,i.aux.inc),t(1):t(end)),1).*p.pop1970;
cases_all=tmp(:,1);
cases_mdr=tmp(:,2);
cases_tbhiv=tmp(:,3);

%HIV Incidence per 1000
tmp = diff(interp1(t,soln(:,i.aux.hiv),t(1):t(end)),1);
hiv_all=1e3*(tmp(:,1)'./pop);

%HIV
tpt_hhc_u5=getfield(model,'tpt_hhc_u5','est');
tpt_hhc=getfield(model,'tpt_hhc','est');

%HHC
pr_onart=getfield(model,'pr_onart','est');
pr_onipt=getfield(model,'pr_onipt','est');

txcov=getfield(model,'txcov','est');
%Arrays
if  (isfield(aux,'llk'))

    results(ii).llk=aux.llk;
end
results(ii).popu= pop;
results(ii).popu_slum= pop_slum;
results(ii).popslumL=popslumL;
results(ii).prev= prev;
results(ii).prev_hi= prev_hi;
results(ii).inc_all= inc_all;
results(ii).inc_mdr= inc_mdr;
results(ii).inc_tbhiv= inc_tbhiv;
results(ii).inc_slum= inc_slum;
results(ii).mort_tbhn= mort_tbhn;
results(ii).mort_tbhiv= mort_tbhiv;
results(ii).notif= notif;
results(ii).cases_all= cases_all;
results(ii).cases_mdr= cases_mdr;
results(ii).cases_tbhiv= cases_tbhiv;
results(ii).incL_all= incL_all;
results(ii).incL_mdr= incL_mdr;
results(ii).incL_tbhiv= incL_tbhiv;
results(ii).pr_onart=pr_onart;
results(ii).pr_onipt=pr_onipt;
results(ii).tpt_hhc_u5=tpt_hhc_u5;
results(ii).tpt_hhc=tpt_hhc;
results(ii).hiv_all=hiv_all;
results(ii).tL=tL(1):tL(end);
results(ii).txcov=txcov;

results(ii).sfin= soln(end-2,1:i.nx);
results(ii).x=x0;
% results(ii).hhc_fac= aux.hhc_fac;
%results(ii).M=aux.M;

%% Run Sims

for ii=1:nsims-1


    ppm.increment();
    x0=xvals(ii,:);

    aux = fx(x0);
    model=aux.model;
    soln=aux.solnL;
    t=soln(:,end);
    tL=aux.solnL(:,end);


    tmp = interp1(t,sum(soln(:,1:i.nstates),2),t(1):t(end));
    pop = (tmp(1:end-1)+tmp(2:end))/2;

    tmp = interp1(t,sum(soln(:,s.slum),2),t(1):t(end));
    pop_slum = (tmp(1:end-1)+tmp(2:end))/2;

    tmp = interp1(tL,sum(soln(:,1:i.nstates),2),tL(1):tL(end));
    popL = (tmp(1:end-1)+tmp(2:end))/2;

    % Prevalence
    %Notif (All)
    tmp  = diff(interp1(t,soln(:,i.aux.notif),t(1):t(end)),1);
    notif=1e5*(tmp(:,1)./pop')';

    %Prevalence
    % tmp = interp1(t,sum(soln(:,intersect(s.prevalent,s.a15p)),2),t(1):t(end));
    tmp = interp1(t,sum(soln(:,s.prevalent),2),t(1):t(end));
    num = (tmp(1:end-1)+tmp(2:end))/2;
    prev= 1e5*(num./pop_ad);

    %Prevalence high risk
    tmp = interp1(t,sum(soln(:,intersect(s.prevalent,s.slum)),2),t(1):t(end));
    num = (tmp(1:end-1)+tmp(2:end))/2;
    prev_hi= 1e5*(num./pop_slum);


    %Incidence
    %Inc (All , MDR, TB/HIV)
    tmp = diff(interp1(t,soln(:,i.aux.inc),t(1):t(end)),1);
    inc=1e5*(tmp./[pop' pop' pop' pop_slum'])';
    inc_all= inc(1,:);
    inc_mdr= inc(2,:);
    inc_tbhiv= inc(3,:);
    inc_slum= inc(4,:);


    %  Mort TB in HIV negative    Mort TB in HIV +
    tmp  = diff(interp1(t,soln(:,i.aux.mort),t(1):t(end)),1);
    mort= 1e5*[tmp(:,1)./pop' ,...
        tmp(:,2)./pop']';

    mort_tbhn= mort(1,:);
    mort_tbhiv= mort(2,:);


    %Incidence
    tmp = diff(interp1(tL,soln(:,i.aux.inc),tL(1):tL(end)),1);
    incL = 1e5.*([tmp(:,1) , tmp(:,2), tmp(:,3) ]./[popL; popL; popL]');
    incL_all=incL(:,1);
    incL_mdr=incL(:,2);
    incL_tbhiv=incL(:,3);

    %cases
    tmp = diff(interp1(t,soln(:,i.aux.inc),t(1):t(end)),1).*p.pop1970;
    cases_all=tmp(:,1);
    cases_mdr=tmp(:,2);
    cases_tbhiv=tmp(:,3);

    %HIV Incidence per 1000
    tmp = diff(interp1(t,soln(:,i.aux.hiv),t(1):t(end)),1);
    hiv_all=1e3*(tmp(:,1)'./pop);

    %HIV
    tpt_hhc_u5=getfield(model,'tpt_hhc_u5','est');
    tpt_hhc=getfield(model,'tpt_hhc','est');

    %HHC
    pr_onart=getfield(model,'pr_onart','est');
    pr_onipt=getfield(model,'pr_onipt','est');

    txcov=getfield(model,'txcov','est');
    %Arrays
    if  (isfield(aux,'llk'))

        results(ii).llk=aux.llk;
    end
    results(ii).popu= pop;
    results(ii).popu_slum= pop_slum;
    results(ii).popslumL=popslumL;
    results(ii).prev= prev;
    results(ii).prev_hi= prev_hi;
    results(ii).inc_all= inc_all;
    results(ii).inc_mdr= inc_mdr;
    results(ii).inc_tbhiv= inc_tbhiv;
    results(ii).inc_slum= inc_slum;
    results(ii).mort_tbhn= mort_tbhn;
    results(ii).mort_tbhiv= mort_tbhiv;
    results(ii).notif= notif;
    results(ii).cases_all= cases_all;
    results(ii).cases_mdr= cases_mdr;
    results(ii).cases_tbhiv= cases_tbhiv;
    results(ii).incL_all= incL_all;
    results(ii).incL_mdr= incL_mdr;
    results(ii).incL_tbhiv= incL_tbhiv;
    results(ii).pr_onart=pr_onart;
    results(ii).pr_onipt=pr_onipt;
    results(ii).tpt_hhc_u5=tpt_hhc_u5;
    results(ii).tpt_hhc=tpt_hhc;
    results(ii).hiv_all=hiv_all;
    results(ii).tL=tL(1):tL(end);
    results(ii).txcov=txcov;


    results(ii).sfin= soln(end-2,1:i.nx);
    results(ii).x=x0;
    % results(ii).hhc_fac= aux.hhc_fac;
    %results(ii).M=aux.M;

end
%%Save files

if (file~=0)
    save_results(results,file,location,runtype);
end

