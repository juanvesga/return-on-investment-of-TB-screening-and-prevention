function save_results(object,type,location,runtype)

filename = sprintf('%s',type,'_',location,'_',runtype,'.mat');
f = fullfile('output',filename);
save(f,'object');

