
function lhd =log_likelihoods(data,data_prev,location,quant)

allpoints=data{location,:};
prevpoints=data_prev{location,:};

if prevpoints(1)~=-99

    allpoints=[allpoints prevpoints(1)];
    
present=(find(~isnan(allpoints)));

noise=exprnd(1/1e6,numel(present),1); % Exponentail random noise to avoid iNF

f={};


if (strcmp(quant,'llk'))
    for g=1:numel(present)
        if (allpoints(g)>1)
            
            f{g}=@(x)log(poisspdf(allpoints(g),x)+noise(g));
            
        else
            
            f{g}=@(x)log(binopdf(...
                round(allpoints(g)*1000),... % Observed events
                1000,...                     % Observed trials
                x)+noise(g));                         % Modelled probability
            
        end
        
    end
else
    for g=1:numel(present)
        
        f{g}= @(x) -sum( ((allpoints(g)+noise(g))-x).^2);
        
    end
    
    
end

lhd = @(model)[...
    f{1}(model(1)) ,...
    f{2}(model(2))  ,...
    f{3}(model(3)) ,...
    f{4}(model(4))  ,...
    f{5}(model(5))  ,...
    f{6}(model(6))  ,...
    f{7}(model(7)) ,...
    f{8}(model(8))  ,...
    f{9}(model(9)) ,...
    f{10}(model(10))];
else
present=(find(~isnan(allpoints)));

noise=exprnd(1/1e6,numel(present),1); % Exponentail random noise to avoid iNF

f={};


if (strcmp(quant,'llk'))
    for g=1:numel(present)
        if (allpoints(g)>1)
            
            f{g}=@(x)log(poisspdf(allpoints(g),x)+noise(g));
            
        else
            
            f{g}=@(x)log(binopdf(...
                round(allpoints(g)*1000),... % Observed events
                1000,...                     % Observed trials
                x)+noise(g));                         % Modelled probability
            
        end
        
    end
else
    for g=1:numel(present)
        
        f{g}= @(x) -sum( ((allpoints(g)+noise(g))-x).^2);
        
    end
    
    
end

lhd = @(model)[...
    f{1}(model(1)) ,...
    f{2}(model(2))  ,...
    f{3}(model(3)) ,...
    f{4}(model(4))  ,...
    f{5}(model(5))  ,...
    f{6}(model(6))  ,...
    f{7}(model(7)) ,...
    f{8}(model(8))  ,...
    f{9}(model(9))];
end





