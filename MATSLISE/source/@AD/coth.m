function x = coth(x)
% Coth for AD objects.

c=coth(x.tc(1));
x.tc=[c,(1-c^2)*x.tc(2),-2*c*(1-c^2)*x.tc(2)^2+(1-c^2)*x.tc(3)];