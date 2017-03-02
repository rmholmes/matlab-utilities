function theta = Heading(lon1,lat1,lon2,lat2);
%This function calculates the rhumb line heading between two
%points

%calculate D:
    'ERRORRRRR!!! This doesnt work'
lon1 = lon1*pi/180;
lon2 = lon2*pi/180;
lat1 = lat1*pi/180;
lat2 = lat2*pi/180;

dpsi = log(tan(pi/4+lon2/2)/tan(pi/4+lon1/2));
theta = atan2(lat2-lat1,dpsi);
end

