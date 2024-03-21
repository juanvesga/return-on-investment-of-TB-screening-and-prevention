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


model=aux.model;

t=aux.soln(:,end);
tL=aux.solnL(:,end);


tmp = interp1(t,sum(aux.soln(:,1:i.nstates),2),t(1):t(end));
pop = (tmp(1:end-1)+tmp(2:end))/2;

tmp = interp1(t,sum(aux.soln(:,s.a0_4),2),t(1):t(end));
pop_0_4 = (tmp(1:end-1)+tmp(2:end))/2;

tmp = interp1(t,sum(aux.soln(:,s.a5_9),2),t(1):t(end));
pop_5_9 = (tmp(1:end-1)+tmp(2:end))/2;

tmp = interp1(t,sum(aux.soln(:,[s.a10_15 s.a15p]),2),t(1):t(end));
pop_10p = (tmp(1:end-1)+tmp(2:end))/2;

tmp = interp1(t,sum(aux.soln(:,s.a15p),2),t(1):t(end));
pop_ad = (tmp(1:end-1)+tmp(2:end))/2;

tmp = interp1(t,sum(aux.soln(:,s.slum),2),t(1):t(end));
pop_slum = (tmp(1:end-1)+tmp(2:end))/2;

tmp = interp1(tL,sum(aux.solnL(:,1:i.nstates),2),tL(1):tL(end));
popL = (tmp(1:end-1)+tmp(2:end))/2;

tmp = interp1(tL,sum(aux.solnL(:,s.slum),2),tL(1):tL(end));
popslumL = (tmp(1:end-1)+tmp(2:end))/2;


% Prevalence
% tmp = interp1(t,sum(aux.soln(:,intersect(s.prevalent,s.a15p)),2),t(1):t(end));
% pre = (tmp(1:end-1)+tmp(2:end))/2;
% prev= 1e5*(pre./pop_ad);

prev=getfield(model,'prev','est');

prev_hi=getfield(model,'prev_hi','est');

%Incidence
% tmp = diff(interp1(t,aux.soln(:,i.aux.inc),t(1):t(end)),1);
% inc = 1e5.*([tmp(:,1) , tmp(:,2), tmp(:,3) ]./[pop; pop; pop]');
inc_all=getfield(model,'inc_all','est');
inc_mdr=getfield(model,'inc_mdr','est');
inc_tbhiv=getfield(model,'inc_tbhiv','est');
inc_slum=getfield(model,'inc_slum','est');
 
%Incidence
tmp = diff(interp1(tL,aux.solnL(:,i.aux.inc),tL(1):tL(end)),1);
incL = 1e5.*([tmp(:,1) , tmp(:,2), tmp(:,3) ]./[popL; popL; popL]');
incL_all=incL(:,1);
incL_mdr=incL(:,2);
incL_tbhiv=incL(:,3);

%cases
tmp = diff(interp1(t,aux.soln(:,i.aux.inc),t(1):t(end)),1).*p.pop1970;
cases_all=tmp(:,1);
cases_mdr=tmp(:,2);
cases_tbhiv=tmp(:,3);

%HIV Incidence per 1000
tmp = diff(interp1(t,aux.soln(:,i.aux.hiv),t(1):t(end)),1);
hiv_all=1e3*[tmp(:,1)'./pop];
hiv_a0_4=1e3*[tmp(:,2)'./pop_0_4];
hiv_a5_9=1e3*[tmp(:,3)'./pop_5_9];
hiv_a10p=1e3*[tmp(:,4)'./pop_10p];

%Notified
% tmp = diff(interp1(t,aux.soln(:,i.aux.notif),t(1):t(end)),1);
% notif = 1e5.*([tmp(:,1) , tmp(:,2)]./[pop; pop]');
% notif_all=notif(:,1);
% notif_mdr=notif(:,2);
notif=getfield(model,'notif','est');


%Mortality
% tmp  = diff(interp1(t,aux.soln(:,i.aux.mort),t(1):t(end)),1);
% mort_tb=1e5*(tmp(:,2)./pop');
% mort_tbhiv=1e5.*(tmp(:,3)./pop');
% mort_tbmdr=1e5.*(tmp(:,4)./pop');
mort_tbhn=getfield(model,'mort_tbhn','est');
mort_tbhiv=getfield(model,'mort_tbhiv','est');

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

results(ii).sfin= aux.soln(end-2,1:i.nx);
results(ii).x=x0;


%% Run Sims

for ii=1:nsims-1


    ppm.increment();
    x0=xvals(ii,:);

    aux = fx(x0);
model=aux.model;

t=aux.soln(:,end);
tL=aux.solnL(:,end);
tL = tL(1):tL(end);

tmp = interp1(t,sum(aux.soln(:,1:i.nstates),2),t(1):t(end));
pop = (tmp(1:end-1)+tmp(2:end))/2;

tmp = interp1(t,sum(aux.soln(:,s.a15p),2),t(1):t(end));
pop_ad = (tmp(1:end-1)+tmp(2:end))/2;

tmp = interp1(t,sum(aux.soln(:,s.slum),2),t(1):t(end));
pop_slum = (tmp(1:end-1)+tmp(2:end))/2;

tmp = interp1(tL,sum(aux.solnL(:,1:i.nstates),2),tL(1):tL(end));
popL = (tmp(1:end-1)+tmp(2:end))/2;

% Prevalence
% tmp = interp1(t,sum(aux.soln(:,intersect(s.prevalent,s.a15p)),2),t(1):t(end));
% pre = (tmp(1:end-1)+tmp(2:end))/2;
% prev= 1e5*(pre./pop_ad);

prev=getfield(model,'prev','est');
prev_hi=getfield(model,'prev_hi','est');

%Incidence
% tmp = diff(interp1(t,aux.soln(:,i.aux.inc),t(1):t(end)),1);
% inc = 1e5.*([tmp(:,1) , tmp(:,2), tmp(:,3) ]./[pop; pop; pop]');
inc_all=getfield(model,'inc_all','est');
inc_mdr=getfield(model,'inc_mdr','est');
inc_tbhiv=getfield(model,'inc_tbhiv','est');
inc_slum=getfield(model,'inc_slum','est');

%Incidence
tmp = diff(interp1(tL,aux.solnL(:,i.aux.inc),tL(1):tL(end)),1);
incL = 1e5.*([tmp(:,1) , tmp(:,2), tmp(:,3) ]./[popL; popL; popL]');
incL_all=incL(:,1);
incL_mdr=incL(:,2);
incL_tbhiv=incL(:,3);

%cases
tmp = diff(interp1(t,aux.soln(:,i.aux.inc),t(1):t(end)),1).*p.pop1970;
cases_all=tmp(:,1);
cases_mdr=tmp(:,2);
cases_tbhiv=tmp(:,3);

%Notified
% tmp = diff(interp1(t,aux.soln(:,i.aux.notif),t(1):t(end)),1);
% notif = 1e5.*([tmp(:,1) , tmp(:,2)]./[pop; pop]');
% notif_all=notif(:,1);
% notif_mdr=notif(:,2);
notif=getfield(model,'notif','est');


%Mortality
% tmp  = diff(interp1(t,aux.soln(:,i.aux.mort),t(1):t(end)),1);
% mort_tb=1e5*(tmp(:,2)./pop');
% mort_tbhiv=1e5.*(tmp(:,3)./pop');
% mort_tbmdr=1e5.*(tmp(:,4)./pop');
mort_tbhn=getfield(model,'mort_tbhn','est');
mort_tbhiv=getfield(model,'mort_tbhiv','est');

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


results(ii).sfin= aux.soln(end-2,1:i.nx);
results(ii).x=x0;


end
%%Save files

if (file~=0)
    save_results(results,file,location,runtype);
end

