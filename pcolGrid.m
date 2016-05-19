function pcolGrid(lonsp,latsp,varargin)
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
hold on;
for ln = -260:lonsp:-80
plot([ln ln],[-30 30]);
end
for lt = -30:latsp:30
plot([-260 -80],[lt lt]);
end
end
