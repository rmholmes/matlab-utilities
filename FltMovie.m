
function FltMovie(hname,dname,VarOp,slice,fname,flts,expname,caxs);
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%
% FltMovie(hname,dname,VarOp,slice,fname,flts);
%
% This function generates a series of plots or a single plot of a
% variable as for pcolPlotNC, except that the float positions (flts)
% from the floats file fname are plotted over the top
%
% INPUTS:
%
% hname,dname,VarOp,slice -> as for pcolPlotNC.
%
% fname = filename of floats file.
%
% flts = float index numbers (vector)
%
% expname = 0 (don't output png's), = 'name' output pngs with
% name_tindex.png...
%
%---------------------------------------------------------------------
%
% Dependencies; pcolPlotNC, FltPos
%
% Ryans ROMS Matlab and netcdf Utilities 19/9/13
%
%---------------------------------------------------------------------
%---------------------------------------------------------------------
    
    
%Derive time limits (for now this is just year days):
   tlims = slice{4};
   
%Get float pos vectors:
flt_pos = FltPos(fname,flts,tlims,0);
tL = length(flt_pos(1,:,1));

%Start time loop:
for t = 1:tL
    cla;
    slice{4} = flt_pos(1,t,1)+1e-5;
    %Plot field:
    pcolPlotNC(hname,dname,VarOp,slice,1);
    
    if (length(caxs) ~= 1)
        caxis(caxs);
    end

    %Plot floats (for now assumes lat/lon):
    hold on;
    plot(flt_pos(:,t,6),flt_pos(:,t,7),'*k');
    
    drawnow;
    if (length(expname)>1)
    %Export if needed:
    if (t<10)
    name = [expname '_00' num2str(t) '.png'];
    elseif (t<100)
    name = [expname '_0' num2str(t) '.png'];
    else
    name = [expname '_' num2str(t) '.png'];
    end
    export_fig(name,'-painters');
    end
end

end


    