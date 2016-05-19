function [un,vn] = ff_nrmcmp(lat, lon, u, v)
% ff_nrmcmp.m -> compute normal component to a velocity section
%
%   call:  [un,vn] = ff_nrmcmp(lat, lon, u, v)
%
%  input:  lat, lon -> latitude/longitude points
%              u, v -> velocity components
%
% output:  un, vn -> normal velocity components
%
% example: [un,vn] = ff_nrmcmp(lat, lon, u, v);
%
% subfunctions: sw_dist.m from seawater
%

%
% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     09-Jan-2010
% modified: 09-Jan-2010
%
% obs:
%

[dist,ang] = sw_dist(lat,lon,'km');

lam = ( lat(1:end-1) + lat(2:end) )/2;
lom = ( lon(1:end-1) + lon(2:end) )/2;
uu  = ( u(:,1:end-1) + u(:,2:end) )/2;
vv  = ( v(:,1:end-1) + v(:,2:end) )/2;

ang = 180-ang;
%  ang = ang;

un = uu.*cos( (ang*pi)/180 ) - vv.*sin( (ang*pi)/180 );
vn = vv.*cos( (ang*pi)/180 ) + uu.*sin( (ang*pi)/180 );

%  un = uu.*cos( (ang*pi)/180 ) + vv.*sin( (ang*pi)/180 );
%  vn = vv.*cos( (ang*pi)/180 ) - uu.*sin( (ang*pi)/180 );