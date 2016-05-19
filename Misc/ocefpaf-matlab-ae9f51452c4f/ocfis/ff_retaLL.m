function [loni,lati,ui,vi] = ff_retaLL(lon,lat,u,v,np);
% ff_retaLL.m -> creates a line for a given lat/lon vector
%
%   call:  [loni,lati,ui,vi] = ff_retaLL(lon,lat,u,v,np);
%
%  input:  lon, lat -> latitude/longitude points
%              u, v -> property (velocity)
%                np -> number of points for interpolated line
%
% output:  loni, lati -> interpolated latitude/longitude points
%              ui, vi -> property (velocity)
%
% example: [output] = ff_retaLL(input);
%
% subfunctions: none
%   MAT-files:  none
%

%
% ff_retaLL.m
% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     09-Jan-2010
% modified: 09-Jan-2010
%
% obs: change to any property instead of [u,v]
%

% linear fit
[P,S] = polyfit(lon,lat,1);

%  longitude increments
dlon = abs( min(lon) - max(lon) );
loni = ( min(lon):dlon/np:max(lon) );

% evaluating the fit
lati = polyval(P,loni);

% computing minimun distance to the fit
uqi = []; vqi = [];
for k = 1:np+1
  % minimun distance point
  dist    = sqrt( (lon-loni(k) ).^2 + (lat-lati(k) ).^2 );
  [dm,km] = min(dist);
  ui(k) = u(km);
  vi(k) = v(km);
end