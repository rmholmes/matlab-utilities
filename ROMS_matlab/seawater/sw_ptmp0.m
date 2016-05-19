function [ptmp] = sw_ptmp0(salt,temp,depth,latitude)

%This function is simply shorthand for
%sw_ptmp(salt,temp,sw_pres(depth,lat),sea level pressure)
%
% It gives the resulting potential temperature.
%
% Ryan Holmes 11-10-13
ptmp = sw_ptmp(salt,temp,sw_pres(depth,latitude),1/9.87*ones(size(depth)));

end
