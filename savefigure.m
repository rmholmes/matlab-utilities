function savefigure(fname)
% This script saves the current figure in .fig and .png files.

saveas(gcf,[fname '.fig']);
export_fig([fname '.png'],'-painters');

end
