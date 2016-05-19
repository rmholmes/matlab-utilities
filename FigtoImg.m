function FigtoImg(directory,type)
%This function turns every .fig file in directory into an image
%file with extension type using export_fig

tmp = dir(directory);
for f=3:length(tmp)
    fname = [directory tmp(f).name];
    if (strcmp(fname((end-3):end),'.fig'))
        outname = [fname(1:(end-3)) type];
        if (exist(outname)~=2)
        openfig(fname,'new','invisible');
        export_fig(outname,'-painters');
        close all;
        end
    end
end

end
