function D = Haversine(lon1,lat1,lon2,lat2);
%This function calculates the distance between two points in meters
%given the longitude and latitudes of the two points.

%Radius of Earth:
Re=6371000.0; %m

%calculate D:
lon1 = lon1*pi/180;
lon2 = lon2*pi/180;
lat1 = lat1*pi/180;
lat2 = lat2*pi/180;

D = 2*Re*asin(sqrt(sin((lat2-lat1)/2)^2+...
                   cos(lat1)*cos(lat2)* ...
                   sin((lon2-lon1)/2)^2));

end

