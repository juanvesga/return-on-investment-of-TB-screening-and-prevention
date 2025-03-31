
function LLK = get_llk(model,data,fns,alpha)

%Number of data categories to fit
cats=(fields(fns));
n_cats=numel(cats);
%n_points=size(data,1);
noise=exprnd(1/1e6);%,n_cats,numel(getfield(model,'inc_all','year')));

llks=zeros(n_cats,1);


for c = 1:n_cats

    field=cats{c};
    dtyr= data{strcmp(data.var,field),"year"}';
    ids=ismember( getfield(model,field,'year'),dtyr) ;
    mod=getfield(model,field,'est');


    if (strcmp(field,"txcov") ||...
            strcmp(field,"pr_onart") ||...
            strcmp(field,"pr_onipt") || ...
            strcmp(field,"prev") || strcmp(field,"prev_hi"))

        mod=mod(ids)/1e5;

    else

        mod=mod(ids);

    end

    dat=data{strcmp(data.var,field),"estimate"}';

    %Weight likelihood if more than 1 point
    weight= (numel(dat) * 1* (numel(dat)>1) ) + ( 1-(numel(dat)>1) );

  

    


    llks_id=zeros(numel(dat),1);

    if(dat~=-99)

        for jj = 1:weight
            llks_id(jj)=fns.(field){jj}(mod(jj)+noise);
        end

        llks(c)=sum(llks_id)/weight;





    else

        llks(c)=0;

    end



end

LLK=llks;

end


