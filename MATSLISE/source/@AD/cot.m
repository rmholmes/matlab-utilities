function x =cot(x)
% Cot for AD objects.

c=cot(x.tc(1));
x.tc=[c,(-1-c^2)*x.tc(2),-2*c*(-1-c^2)*x.tc(2)^2+(-1-c^2)*x.tc(2)];