function [lon,lat] = Midpoint(lon1,lat1,lon2,lat2);
%This function calculates the midpoint between two points
%given the longitude and latitudes of the two points.

%Radius of Earth:
Re=6371000.0; %m

%calculate D:
lon1 = lon1*pi/180;
lon2 = lon2*pi/180;
lat1 = lat1*pi/180;
lat2 = lat2*pi/180;

[lat,lon] = meanm([lat1 lat2],[lon1 lon2]);

lon = lon/pi*180;
lat = lat/pi*180;
end

