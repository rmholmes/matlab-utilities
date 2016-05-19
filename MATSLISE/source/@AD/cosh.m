function x = cosh(x)
% Cosh for AD objects.
c=cosh(x.tc(1));
s=sinh(x.tc(1));
 x.tc = [c,s*x.tc(2),c*x.tc(2)^2+s*x.tc(3)];
