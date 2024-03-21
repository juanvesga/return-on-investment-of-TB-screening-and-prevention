function save_fig(object,type,location,runtype,folder)


fig = get(groot,'CurrentFigure');
if ~isempty(fig)
filename = sprintf('%s',type,'_',location,'_',runtype,'.fig');
f = fullfile(folder,filename);
savefig(object,f,'compact'); 
end
