function pcolMovNC(hname,dname,VarOp,slice,varargin)
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%
% This function makes and plays a movie of a variable taking the
% same arguments as pcolPlotNC with the time slice element being
% the sequence of times wanted.
%
% INPUTS: standard pcolPlotNC
%
% vargin = clim,axs,oname
%
% clim = caxis limits, = 0 -> automatic, = 1 -> max/min = vec -> set
%
% axs = axis limits (4 element vector, or 0 for automatic);
% 
% oname = file name base for exporting through export_fig
%---------------------------------------------------------------------
%
% Dependencies; none. 
%
% Ryans ROMS Matlab and netcdf Utilities 15/8/13
%
%---------------------------------------------------------------------
%---------------------------------------------------------------------

    nargin
    if (nargin == 4)
        csw = 0;
        exp = 0;
    else
        if (length(varargin{1}) == 2)
            csw = 2;
        else
            csw = varargin{1};
            
        end
    end
    if (nargin >= 6)
        if (length(varargin{2}) == 4)
            saxs = 1;
            axs = varargin{2};
        else
            saxs = 0;
        end
    end
    if (nargin == 7)
        exp = 1;
        oname = varargin{3};
    else
        exp = 0;
    end
    
    tlim = slice{4};
    if (round(tlim(1)) == tlim(1))
        tvec = tlim(1):1:tlim(2);
    else
        ncid = netcdf.open(hname,'NC_NOWRITE');
        tmp = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'ocean_time'), ...
                             [0],[2],'double')/86400;
        dt = tmp(2)-tmp(1);
        tvec = tlim(1):dt:tlim(2);
    end
    
    set(gcf,'NextPlot','replacechildren');
for ti = 1:length(tvec)
    slice{4} = tvec(ti);
    [Var,iVar,jVar] = pcolPlotNC(hname,dname,VarOp,slice,1);
    if (csw == 0)
        colorbar;
    elseif (csw == 1)
        caxis([min(min(Var)) max(max(Var))]);
    else
        caxis(varargin{1});
    end
    drawnow;
    if (exp == 1)
        if (ti<10)
            name = [oname '_00' num2str(ti) '.png'];
        elseif (ti<100)
            name = [oname '_0' num2str(ti) '.png'];
        else
            name = [oname '_' num2str(ti) '.png'];
        end
        export_fig(name,'-painters');
    end
end

end
    