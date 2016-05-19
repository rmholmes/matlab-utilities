function x = sinh(x)
% Sinh for AD objects.

s=sinh(x.tc(1));
x.tc=[s,x.tc(2)*cosh(x.tc(1)),x.tc(3)*cosh(x.tc(1))+x.tc(2)^2*sinh(x.tc(1))];