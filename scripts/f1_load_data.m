% Load Contact matrices (synth matricces in Prem et al 2021 PlosComBiol
file=sprintf('%s','/data/contact_',country,'.csv');
contact_matrix = readtable(file);
C=contact_matrix{:,:};

%Load population
population_table=readtable('/data/population.csv');
cols=size(population_table,2);
id=find(strcmp(population_table{:,cols},country));
pop_vec= population_table{id,1:(cols-1)};

%% Proocess contact matrix
%Symmetric matrix
c_sym = (C + C')/2; 

% Per capita per year symmetric contact matrix
c_mat = c_sym/max(c_sym(:));%./repmat(pop_vec',1,cols-1);

%% Load HIV points for HIV acqusition rate
hivrts=load_data('data/hivrates2.xlsx');
tmp=hivrts{country,:};
id=find(tmp==-99);
tmp(id)=[];

period=5;
hivpoints=interp1([1980 (1990+numel(id)):2022],[0 tmp],1980:1:2022,'pchip');

if strcmp(country,'ZAF')

hivdecline= 0.015;%((hivpoints(end-period)-hivpoints(end))/hivpoints(end-period))/period;

else

hivdecline= ((hivpoints(end-period)-hivpoints(end))/hivpoints(end-period))/period;

end


tperiod=2022:1:(endyear+5);
hivext=tperiod*0;

% Use compound interest to find HIV rate in coming years
for ii=1:numel(tperiod)

    hivext(ii)= hivpoints(end)  - ((hivpoints(end)*(1 + hivdecline)^ii)-hivpoints(end));

end

hivpoints =[hivpoints , hivext(2:end)];

% figure;plot(1980:1:(endyear+5), hivpoints*1000)



% ylim([-0.1,0.2])

%% Load Fiting data

target_data=readtable('/data/who_data.csv');
id=find(strcmp(target_data{:,2},country));
target_data_full=target_data(id,:);
target_data=target_data(id,[2:4,end]);
next_rs=[1:7]+size(target_data,1);
target_data(next_rs,:)=target_data(next_rs-7,:);
target_data_full(next_rs,:)=target_data_full(next_rs-7,:);

% Prevalence surveys
prev_data  =load_data(strcat(pwd,'/data/prev_surveys.xlsx'));
prev_data{country,"year"};
target_data(next_rs(1),'year')=table(prev_data{country,"year"});
target_data(next_rs(1),'estimate')=table(prev_data{country,"prevalence"});
target_data(next_rs(1),'var')={'prev'};

target_data_full(next_rs(1),'year')=table(prev_data{country,"year"});
target_data_full(next_rs(1),'estimate')=table(prev_data{country,"prevalence"});
target_data_full(next_rs(1),'low')=table(prev_data{country,"low"});
target_data_full(next_rs(1),'high')=table(prev_data{country,"high"});
target_data_full(next_rs(1),'var')={'prev'};

% Prevalence high risk surveys
prev_data  =load_data(strcat(pwd,'/data/prevhi_surveys.xlsx'));
prev_data{country,"year"};
target_data(next_rs(2),'year')=table(prev_data{country,"year"});
target_data(next_rs(2),'estimate')=table(prev_data{country,"prevalence"});
target_data(next_rs(2),'var')={'prev_hi'};

target_data_full(next_rs(2),'year')=table(prev_data{country,"year"});
target_data_full(next_rs(2),'estimate')=table(prev_data{country,"prevalence"});

if prev_data{country,"low"}==-99
target_data_full(next_rs(2),'low')= table(max(-99,prev_data{country,"prevalence"} * 0.8));
target_data_full(next_rs(2),'high')=table(max(-99,prev_data{country,"prevalence"} * 1.2));
else
target_data_full(next_rs(2),'low')=table(prev_data{country,"low"});
target_data_full(next_rs(2),'high')=table(prev_data{country,"high"});
end
target_data_full(next_rs(2),'var')={'prev_hi'};



%% HIV

tbhiv_data  =readtable('data/tbhiv_data.csv');
id=find(strcmp(tbhiv_data{:,1},country));
tbhiv_data=tbhiv_data(id,:);
tbhiv_data(:,"estimate")=table(tbhiv_data{:,"estimate"});
target_data(next_rs(3),'year')=tbhiv_data(1,'year');
target_data(next_rs(3),'estimate')=tbhiv_data(1,'estimate');
target_data(next_rs(3),'var')=tbhiv_data(1,'indicator');

target_data_full(next_rs(3),'year')=tbhiv_data(1,'year');
target_data_full(next_rs(3),'estimate')=tbhiv_data(1,'estimate');
target_data_full(next_rs(3),'low')=table(-99);
target_data_full(next_rs(3),'high')=table(-99);
target_data_full(next_rs(3),'var')=tbhiv_data(1,'indicator');


target_data(next_rs(4),'year')=tbhiv_data(2,'year');
target_data(next_rs(4),'estimate')=tbhiv_data(2,'estimate');
target_data(next_rs(4),'var')=tbhiv_data(2,'indicator');

target_data_full(next_rs(4),'year')=tbhiv_data(2,'year');
target_data_full(next_rs(4),'estimate')=tbhiv_data(2,'estimate');
target_data_full(next_rs(4),'low')=table(-99);
target_data_full(next_rs(4),'high')=table(-99);
target_data_full(next_rs(4),'var')=tbhiv_data(2,'indicator');

% Incidence High Risk
hrdata  =load_data(strcat(pwd,'/data/high_risk_inc_data.xlsx'));
hrdata{country,"year"};
target_data(next_rs(5),'year')=table(hrdata{country,"year"});
target_data(next_rs(5),'estimate')=table(hrdata{country,"incidence"});
target_data(next_rs(5),'var')={'inc_slum'};

target_data_full(next_rs(5),'year')=table(hrdata{country,"year"});
target_data_full(next_rs(5),'estimate')=table(hrdata{country,"incidence"});

if hrdata{country,"low"}==-99
target_data_full(next_rs(5),'low')= table(max(-99,hrdata{country,"incidence"} * 0.8));
target_data_full(next_rs(5),'high')=table(max(-99,hrdata{country,"incidence"} * 1.2));
else
target_data_full(next_rs(5),'low')=table(hrdata{country,"low"});
target_data_full(next_rs(5),'high')=table(hrdata{country,"high"});
end
target_data_full(next_rs(5),'var')={'inc_slum'};

% HHC TPT
hhc_data  =readtable('data/hhc_data.csv');
id=find(strcmp(hhc_data{:,1},country));
hhc_data=hhc_data(id,:);
hhc_data(:,"estimate")=table(hhc_data{:,"estimate"});
target_data(next_rs(6),'year')=hhc_data(1,'year');
target_data(next_rs(6),'estimate')=hhc_data(1,'estimate');
target_data(next_rs(6),'var')=hhc_data(1,'indicator');

target_data_full(next_rs(6),'year')=hhc_data(1,'year');
target_data_full(next_rs(6),'estimate')=hhc_data(1,'estimate');
target_data_full(next_rs(6),'low')=table(-99);
target_data_full(next_rs(6),'high')=table(-99);
target_data_full(next_rs(6),'var')=hhc_data(1,'indicator');


target_data(next_rs(7),'year')=hhc_data(2,'year');
target_data(next_rs(7),'estimate')=hhc_data(2,'estimate');
target_data(next_rs(7),'var')=hhc_data(2,'indicator');

target_data_full(next_rs(7),'year')=hhc_data(2,'year');
target_data_full(next_rs(7),'estimate')=hhc_data(2,'estimate');
target_data_full(next_rs(7),'low')=table(-99);
target_data_full(next_rs(7),'high')=table(-99);
target_data_full(next_rs(7),'var')=hhc_data(2,'indicator');


%Remove Notification data and then edit plot_targets replace with
%notificationcoverage distribution
notifiedrate=target_data{strcmp(target_data.var,'notif'),"estimate"}';
toDelete = strcmp(target_data.var, 'notif');
target_data(toDelete,:) = [];
