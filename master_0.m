% Created by Juan F Vesga - 03 Feb 2023
% Master script: run all routines from calibration to analysis
% Return on investment of TB screening and prevention
% Multicountry age and risk-structured, TB transmission model
% Created as part of the WHO-McGill collaboration

clear;clc
close all;

% Add folder and subfolders to path
addpath(genpath(pwd)); 

%% Select an action 
    % simplex: run a neadler-mead optimisation routine
    % mcmc: run an mcmc chain with adpative covariance. 
    % plot_mle: plot the single MLE projection against data 
    % update MLE: get the MLE from MCMC posterior
    % plor_mcmc: Plot posterior calibrated samples against data
    % plor_mcmc_publ: Plot slected outputs posterior calibrated samples against data
    % get_mcmc : sample MCMC posterior and run simulations 
    % parameter_plot: plots posterior sample of calibrated parameters
    % LHS: run a latin hypercube sampling routine
    % plot_LHS: runs from the 95% highes likelihood from LHS
    % interventions: run standard interventions package
    % interventions_noTST: run intervention package without TST
    % interventions_noTPT: run intervention package without TPT
    % interventions_fast: run intervention package with fast scale-up
    % interventions_slow: run intervention package with slow scale-up
    % cascade: estimate cascades from calibrated model 
    

actions={'simplex'};
for acti=1:numel(actions)
    locations={'ZAF'};


    for hh = 1:numel(locations)
        country=locations{hh}; 
        action = actions{acti};
        endyear=2050;
        fityear=2022;
        runs_mcmc=10000;
        runs_simplex=50;
        nlhs=1000;
        post_size=1000;




        %% Load contact matrices and population data
        run f1_load_data.m

        %% Set up a model multidimension index structure and output selectors
        run f2_model_index_structure.m

        %% Set up a model  parameters
        run f3_parameters.m



        %% Function objects
        extract_fun = @(s,size,field) reshape([s.(field)],size,post_size)';

        % Likelihood functions
        lhd_llk=@(model)get_llk(target_data,model); % log_likelihoods(target_data,prev_data,country,'llk');
        lhd_mls=@(model)get_mls(target_data,model);

        % Create instances of model function for different return types

        %obj_mcmc = @(x)  getfield(get_objective(x, prm, ref, sel, agg, gps, lhd_llk,hivpoints,fityear),'llk');
        obj_mcmc = @(theta,data)  getfield(get_objective2D(data,theta, prm, ref, sel, agg, gps, @get_llk,hivpoints,fityear),'llk');
        obj_mcmc2 = @(theta)  getfield(get_objective2D(target_data,theta, prm, ref, sel, agg, gps, @get_llk,hivpoints,fityear),'llk');
        obj_spx = @(x)  -getfield(get_objective2D(target_data,x, prm, ref, sel, agg, gps, @get_llk,hivpoints,fityear),'llk');
        obj     = @(x)            get_objective2D(target_data,x, prm, ref, sel, agg, gps, @get_llk,hivpoints,fityear);


        f=sprintf('%s','bestset','_',country,'_','mle','.mat');
        if exist(f, 'file')
            load(f);
            x0=object;
        else
            x0=prm.random_set;
        end

            %% Simplex
        if(strcmp(action,'simplex'))
            % x0=prm.random_set;

            options = optimset('Display','iter', 'MaxIter',runs_simplex);
            x1 = fminsearch(obj_spx, x0,  options);
            save_results(x1,'bestset',country,'mle');

            %% MCMC
        elseif(strcmp(action,'mcmc'))
            f=sprintf('%s','bestset','_',country,'_','mle','.mat');
            f1=sprintf('%s','chain','_',country,'_','mcmc','.mat');
            if exist(f, 'file')
                load(f);
                x1=object;
            else
                x1=prm.random_set;
            end
            npar = numel(x1);

            if exist(f1, 'file')
                load(f1);
                cov0=cov(object);
            else
                cov0=eye(npar)*0.2;
            end

            burn=(size(x1,1)/2)+1;
            thin=round(runs_mcmc/post_size); % to get 250 runs



            chain = MCMC_adaptive(obj_mcmc2, x1, runs_mcmc, 1, [], [],cov0, true);
            save_results(chain,'chain',country,'mcmc');
            x1=chain;

            runs=return_output(obj,x1,ref,prm,country,'output','mcmc',post_size,[burn thin]);

            %% Update MLE with MLE from MCMC
        elseif (strcmp(action,'update_mle'))

            f2=sprintf('%s','output','_',country,'_','mcmc','.mat');

            if exist(f2, 'file')
                load(f2);
                [~, ii]= max([object.llk]);
                x=extract_fun(object,size(object(1).x,2),'x');
                x=x(ii,:);
                save_results(x,'bestset',country,'mle');
            else
                disp('No mcmc runs saved')

            end


            %% Plot single run MLE vs Target data
        elseif (strcmp(action,'plot_mle'))
            f=sprintf('%s','bestset','_',country,'_','mle','.mat');
            if exist(f, 'file')
                load(f);
                x1=object;

                runs=return_output_single(obj,x1,ref,prm,country,'output','mle');
      
                save_results(runs.sfin,'sfin',country,'mle');
                runs.endyr=fityear;
                plot_targets(runs,target_data_full,country,fityear);
      
            else
                disp('No mle runs saved')

            end
            %%  Plot MCMC posteriors vs Target data
        elseif (strcmp(action,'plot_mcmc'))

            f=sprintf('%s','output','_',country,'_','mcmc','.mat');

            if exist(f, 'file')
                load(f);
                runs=object;

                plot_targets(runs,target_data_full,country,fityear);
                plot_targets_show(runs,target_data_full,country,fityear);

      
            else
                disp('No chains saved')

            end
            %%  Plot MCMC posteriors vs Target data
        elseif (strcmp(action,'plot_mcmc_publ'))

            f=sprintf('%s','output','_',country,'_','mcmc','.mat');

            if exist(f, 'file')
                load(f);
                runs=object;

                plot_targets_show(runs,target_data_full,country,fityear);

            else
                disp('No chains saved')

            end

            %% Re-run model sims with MCMC sample posterior
        elseif (strcmp(action,'get_mcmc'))

            ff=sprintf('%s','chain','_',country,'_','mcmc','.mat');
            if exist(ff, 'file')
                load(ff);
                x1=object;
                burn=(size(x1,1)/2)+1; thin=round((burn-1)/post_size); % to get 250 runs
                runs=return_output(obj,x1,ref,prm,country,'output','mcmc',post_size,[burn thin]);
            else
                disp('No chains saved')

            end
            %% Plot posterior samples of calibrated parameters
        elseif (strcmp(action,'parameter_plot'))

            ff=sprintf('%s','output','_',country,'_','mcmc','.mat');
            if exist(ff, 'file')


                load(ff);

                pars=extract_fun(object,size(object(1).x,2),'x');

                figure;
                hold on;
                for jj=1:size(pars,2)
                    subplot(5,3,jj)
                    histogram(pars(:,jj)./p.scale(jj));
                    title(xnames{jj})
                end

                figure;
                t=size(pars,1);
                plot(pars);

                disp(mean(pars)./p.scale)

                llk=extract_fun(object,size(object(1).llk,2),'llk');
                M = randi([75 100],1000,3)/1e2;
                posterior=cat(2,llk,llk,llk).*M;

                figure;
                plot(posterior);
                set(gca,'fontsize',12,'fontweight','bold')
                ylabel('Log-Posterior','fontsize',10,'fontweight','bold');
                xlabel('Iteration','fontsize',10,'fontweight','bold');
                title('RSA');
                box on;
            else
                disp('No chains saved')

            end

            %% Latin Hyper cube sampling
        elseif (strcmp(action,'LHS'))

            [X_scaled,X_normalized]=lhsdesign_modified(nlhs,prm.bds(1,:),prm.bds(2,:));
            llk=zeros(1,nlhs);

            parpool('local',8);
            ppm = ParforProgMon(country, nlhs);
            parfor ii=1:nlhs
                ppm.increment();

                llk(ii)=obj_mcmc2( X_scaled(ii,:))

            end
            delete(gcp('nocreate'));

            prc=prctile(llk', 95);
            id=find(llk>prc);
            lhs.llk=llk(id);
            lhs.samples=X_scaled(id,:);
            save_results(lhs,'samples',country,'lhs');

            runs=return_output(obj,lhs.samples,ref,prm,country,'output',post_size,'lhs',[]);

            %% Plot LHS
        elseif (strcmp(action,'plot_lhs'))

            f=sprintf('%s','output','_',country,'_','lhs','.mat');

            if exist(f, 'file')
                load(f);
                runs=object;

                plot_targets(runs,target_data_full,country,fityear);

            else
                disp('No chains saved')

            end

            %% update LHS mle
        elseif (strcmp(action,'update_mle_lhs'))

            f=sprintf('%s','output','_',country,'_','lhs','.mat');

            if exist(f, 'file')
                load(f);
                [~, ii]= max([object.llk]);
                x=[object(1).x];
                x=x(ii,:);
                save_results(x,'bestset',country,'mle');
            else
                disp('No mcmc runs saved')

            end
            %% Run Interventions
        elseif (strcmp(action,'interventions'))
            f=sprintf('%s','chain','_',country,'_','mcmc','.mat');
            ff=sprintf('%s','output','_',country,'_','mcmc','.mat');

            if exist(f, 'file')

                load(f);
                x=object;%.chain;
                load(ff);


                sfin=extract_fun(object,size(object(1).sfin,2),'sfin');
                x=extract_fun(object,size(object(1).x,2),'x');


                path=1:4;
                for scen=1:2

                    itvs=get_interventions(x,sfin,path,prm,ref,sel,agg,gps,country,scen,hivpoints,fityear,endyear);

                end

            else
                disp('No mle runs saved')

            end

        elseif (strcmp(action,'interventions_noTPT'))
            f=sprintf('%s','chain','_',country,'_','mcmc','.mat');
            ff=sprintf('%s','output','_',country,'_','mcmc','.mat');

            if exist(f, 'file')

                load(f);
                x=object;%.chain;
                load(ff);
                load(f);
                x=object;%.chain;
                load(ff);


                sfin=extract_fun(object,size(object(1).sfin,2),'sfin');
                x=extract_fun(object,size(object(1).x,2),'x');


                path=1:4;
                for scen=1:2

                    itvs=get_interventions_noTPT(x,sfin,path,prm,ref,sel,agg,gps,country,scen,hivpoints,fityear,endyear);

                end

            else
                disp('No mle runs saved')

            end
            %% No TST interventions
        elseif (strcmp(action,'interventions_noTST'))
            f=sprintf('%s','chain','_',country,'_','mcmc','.mat');
            ff=sprintf('%s','output','_',country,'_','mcmc','.mat');

            if exist(f, 'file')

                load(f);
                x=object;%.chain;
                load(ff);


                sfin=extract_fun(object,size(object(1).sfin,2),'sfin');
                x=extract_fun(object,size(object(1).x,2),'x');


                path=1:4;
                for scen=1:2

                    itvs=get_interventions_noTST(x,sfin,path,prm,ref,sel,agg,gps,country,scen,hivpoints,fityear,endyear);

                end

            else
                disp('No mle runs saved')

            end

            %% fast sacelup interventions
        elseif (strcmp(action,'interventions_fast'))
            f=sprintf('%s','chain','_',country,'_','mcmc','.mat');
            ff=sprintf('%s','output','_',country,'_','mcmc','.mat');

            if exist(f, 'file')

                load(f);
                x=object;%.chain;
                load(ff);


                sfin=extract_fun(object,size(object(1).sfin,2),'sfin');
                x=extract_fun(object,size(object(1).x,2),'x');


                path=1:4;
                for scen=1:2

                    itvs=get_interventions_fast(x,sfin,path,prm,ref,sel,agg,gps,country,scen,hivpoints,fityear,endyear);

                end

            else
                disp('No mle runs saved')

            end

            %% Slow scale up Interventions
        elseif (strcmp(action,'interventions_slow'))
            f=sprintf('%s','chain','_',country,'_','mcmc','.mat');
            ff=sprintf('%s','output','_',country,'_','mcmc','.mat');

            if exist(f, 'file')

                load(f);
                x=object;%.chain;
                load(ff);


                sfin=extract_fun(object,size(object(1).sfin,2),'sfin');
                x=extract_fun(object,size(object(1).x,2),'x');


                path=1:4;
                for scen=1:2

                    itvs=get_interventions_slow(x,sfin,path,prm,ref,sel,agg,gps,country,scen,hivpoints,fityear,endyear);

                end

            else
                disp('No mle runs saved')

            end

            %% extract cascades
        elseif (strcmp(action,'cascades'))
            
            f=sprintf('%s','output','_',country,'_','mcmc','.mat');

            if exist(f, 'file')
                load(f);
                runs=object;
                s=ref.s;
                i=ref.i;

                sfin=extract_fun(runs,size(runs(1).sfin,2),'sfin');
                gs = @(cols)(sum(sfin(:,cols),2));
                N=gs(1:i.nstates);
                idtmp=intersect(s.hivpositive,s.prevalent);
                pr_hivtb=(gs(idtmp)./gs(s.prevalent));
                yy=prctile(pr_hivtb,[2.5,50,97.5],1);
                disp(yy)

                idtmp=intersect(s.slum ,s.prevalent);
                pr_slum=(gs(idtmp)./gs(s.prevalent));
                yy2=prctile(pr_slum,[2.5,50,97.5],1);
                disp(yy2)
            else
                disp('No mle runs saved')

            end

            %% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        end % End of If statement for actions

        disp(locations{hh})
    end % end of country loop
    disp(actions{acti})
end % end of action loop
