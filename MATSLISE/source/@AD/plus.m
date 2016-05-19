function y = plus(x,y)

% Add AD objects.
% Error checking is poor: it assumes that if X or Y is not a AD
% then it is a (scalar) numeric object.
if ~isa(x,'AD') %assume X is a number, add it to constant term.
  y.tc(1) = x + y.tc(1);
elseif ~isa(y,'AD') %assume Y is a number, swap, do same.
  ydup = y; y = x;
  y.tc(1) = ydup + x.tc(1);
else
  y.tc = x.tc + y.tc;
end
