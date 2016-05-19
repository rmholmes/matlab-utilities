function y = times(x,y)
% Multiply AD objects.
% Error checking : it assumes that if X or Y is not a AD
% then it is a (scalar) numeric object.
if ~isa(x,'AD') %assume X is a number
  y.tc = x * y.tc;
elseif ~isa(y,'AD') %assume Y is a number
  ydup = y; y = x;
  y.tc = ydup * x.tc;
else
  y.tc= [x.tc(1)*y.tc(1),x.tc(2)*y.tc(1)+x.tc(1)*y.tc(2),x.tc(3)*y.tc(1)+2*x.tc(2)*y.tc(2)+x.tc(1)*y.tc(3)];
end
