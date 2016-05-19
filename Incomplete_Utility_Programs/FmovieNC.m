function A = FmovieNC(fname_VARS,fname_DIAS,VarOp,slice,clim)
%This function makes and plays a movie of the Variable given by VARID
%taken from a netcdf file. The inputs are the same as pcolPlotNC
%(except the t slice entry should not be included). 
%clim is a two
%element vector containing the caxis.

    fig1 = figure;%('visible','off');
    set(fig1, 'Position', get(0,'Screensize')); % Maximize figure. 
    winsize = get(fig1,'Position');
    ncid = netcdf.open(fname_VARS,'NC_NOWRITE');
    ID = netcdf.inqDimID(ncid,'timeD');
    [temp numframes] = netcdf.inqDim(ncid,ID);
    netcdf.close(ncid);
    A = moviein(numframes,fig1,winsize);
    set(fig1,'NextPlot','replacechildren');
for t=1:numframes
    [Var,iVar,jVar] = pcolPlotNC(fname_VARS,fname_DIAS,VarOp,[slice t],1);
    caxis(clim);
    A(:,t)=getframe(fig1);
end
end
    