
function [lon,lat] = lonlat(x,y)
% this function roughly calculates the cartesian coordinates given
% longitude latitude coordinates in the equatorial Pacific
Re = 6378000;CF = pi/180*Re;

lat = y/CF;
lon = x./cos(pi/180*lat)/CF-132;
end
