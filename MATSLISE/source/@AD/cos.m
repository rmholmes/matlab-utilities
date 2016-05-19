function x = cos(x)
% COS for AD objects.
c=cos(x.tc(1));
x.tc = [c,-sin(x.tc(1))*x.tc(2),-c*(x.tc(2)^2)-sin(x.tc(1))*x.tc(3)];
end
