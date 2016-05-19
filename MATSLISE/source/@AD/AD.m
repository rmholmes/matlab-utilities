function x =  AD(a)
% Constructor for AD class.
switch nargin
case 0
  x.tc = [0, 1, 0];
  x = class(x, 'AD');
case 1
if ~(isnumeric(a) && isscalar(a))
    error('AD: Invalid input A')
else
    x.tc = [a, 1, 0];
end
x=class(x,'AD');
end
