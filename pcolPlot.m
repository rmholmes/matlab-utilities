function h = pcolPlot(lon,lat,F,paxes)
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%
% This function plots the pcolor plot for field F with x = lon and y =
% lat
%
% INPUTS:
% 
% lon = longitude field (2D)
%
% lat = latitude field (2D)
%
% F = field (2D)
%
%---------------------------------------------------------------------
%
% Dependencies; none. 
%
% Ryans ROMS Matlab and netcdf Utilities 17/7/13
%
%---------------------------------------------------------------------
%---------------------------------------------------------------------
    if (nargin==4)
        h = pcolor(lon,lat,F,'Parent',paxes);
        shading(paxes,'flat');
    else
        h = pcolor(lon,lat,F);
% $$$         colormap(dark);
        colormap(flipud(lbmap(51,'RedBlue')));
        colorbar;
        shading flat;
    end
end