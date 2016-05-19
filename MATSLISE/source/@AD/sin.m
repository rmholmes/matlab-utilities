function x = sin(x)
% SIN for AD objects.
s=sin(x.tc(1));
x.tc = [s,cos(x.tc(1))*x.tc(2),-s*(x.tc(2)^2)+cos(x.tc(1))*x.tc(3)];