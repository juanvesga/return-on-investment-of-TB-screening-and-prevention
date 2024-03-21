
function LLK = get_llk(data,model)

%Number of data categories to fit
cats=unique(data{:,"var"});
n_cats=numel(cats);
%n_points=size(data,1);
noise=exprnd(1/1e6,n_cats,numel(getfield(model,'inc_all','year')));

llks=zeros(n_cats,1);


for c = 1:n_cats

    field=cats{c};
    dtyr= data{strcmp(data.var,field),"year"}';
    ids=ismember( getfield(model,field,'year'),dtyr) ;
    mod=getfield(model,field,'est'); mod=mod(ids);
    dat=data{strcmp(data.var,field),"estimate"}';
    weight=numel(dat);

    if(dat~=-99)
        if (strcmp(field,"txcov") ||...
                strcmp(field,"pr_onart") ||...
                strcmp(field,"pr_onipt") || ...
                strcmp(field,"prev") || strcmp(field,"prev_hi"))

            llks(c)=sum(log(binopdf(...
                round(dat*1e2),... % Observed events
                1e2,...                     % Observed trials
                min(mod/1e5,1))+noise(1:numel(mod))))/weight;                         % Modelled probability

        else

            if (strcmp(field,"inc_all"))

                llks(c)=sum( ...
                    (log(nbinpdf(...
                    round(dat),... % Observed events
                    1e5,...                     % Observed trials
                    min(mod/1e5,1))+noise(1:numel(mod)))) ...
                    )/weight;                         % Modelled probability
            else

                llks(c)=(sum( log(poisspdf(round(dat),mod) + ...
                    noise(c,numel(mod))) ))/weight;
            end


        end

    else
        llks(c)=0;
    end

end

LLK=llks;

end


