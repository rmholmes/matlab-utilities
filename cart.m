
function [x,y] = cart(lon,lat)
% this function roughly calculates the cartesian coordinates given
% longitude latitude coordinates in the equatorial Pacific
Re = 6378000;CF = pi/180*Re;

y = CF*lat;
x = CF*(lon+132).*cos(pi/180*lat);
end
