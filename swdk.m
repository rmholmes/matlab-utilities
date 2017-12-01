function [frac] = swdk(Z)
% This function returns the fraction of solar radiation penetrating
% below the specified depth Z (any shape) for ROMS.

mu1 = 0.35;
mu2 = 23.0;
r1  = 0.58;
Zscale = -1;

fac1 = Zscale/mu1;
fac2 = Zscale/mu2;
fac3 = r1;
frac = exp(Z*fac1)*fac3 + exp(Z*fac2)*(1-fac3);

end



