
function MLS = get_mls(data,model)

%Number of data categories to fit
cats=unique(data{:,"var"});
n_cats=numel(cats);
%n_points=size(data,1);
noise=exprnd(1/1e6,n_cats,numel(getfield(model,'inc_all','year')));

mlss=zeros(n_cats,1);


for c = 1:n_cats

    field=cats{c};
    dtyr= data{strcmp(data.var,field),"year"}';
    ids=ismember( getfield(model,field,'year'),dtyr) ;
    mod=getfield(model,field,'est'); mod=mod(ids);
    dat=data{strcmp(data.var,field),"estimate"}';


    if(dat~=-99)
        if (strcmp(field,"pr_onart") || strcmp(field,"pr_onipt")...
                ||strcmp(field,"prev") || strcmp(field,"prev_hi"))

            mlss(c)= -sum( ((dat*1e5+ noise(c,numel(mod)))-mod).^2)/numel(dtyr);

        else

            mlss(c)= -sum( ((dat+ noise(c,numel(mod)))-mod).^2)/numel(dtyr);

        end

    else

        mlss(c)=0;
    end



end

MLS=mlss;
