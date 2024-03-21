function out=prcc(pars,outcome)


dataset=cat(2,pars,outcome);

corrs=partialcorr(dataset);

out=corrs(end,1:end-1);
end