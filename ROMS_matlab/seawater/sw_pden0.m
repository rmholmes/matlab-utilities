function [pden] = sw_pden0(salt,temp,depth,latitude)

%This function is simply shorthand for
%sw_pden(salt,temp,sw_pres(depth,lat),sea level pressure)-1000
%
% It gives the resulting potential density anomaly.
pden = sw_pden(salt,temp,sw_pres(depth,latitude),1/9.87*ones(size(depth)))-1000;

end
