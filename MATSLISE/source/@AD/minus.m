function y = minus(x,y)
% Subtract AD objects.

if ~isa(x,'AD') %assume X is a number.
  y.tc = -y.tc;
  y.tc(1) = x + y.tc(1);
elseif ~isa(y,'AD') %assume Y is a number.
  x.tc(1) = x.tc(1) - y;
  y = x;
else
  y.tc = x.tc - y.tc;
end
