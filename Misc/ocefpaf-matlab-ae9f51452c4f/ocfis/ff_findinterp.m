function [y_x] = ff_findinterp(x,y,x_y)
% ff_findinterp - find a value using linear interpolation
%
%   use:  [y_x] = ff_findinterp(x,y,y_x);
% input:
%            x: vector 1 (temperature)
%            y: vector 2 (depth)
%          x_y: x value to interpolate
%
% output:
%          y_x: y(x)
%
% example:
%    [d10] = findinterp(temp,depth,10);
%     finds depth where temperature is equal to 10
%
% other m-files required: near.m
%
% findinterp.m
% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     09-Apr-2009
% modified: 09-Apr-2009
%
% obs: improve "if's"
%

fx = find(isnan(x));   fy   = find(isnan(y));
f  = [fx fy]; y(f)=[]; x(f) = [];

idx = near(x,x_y,1); % find nearest value to x_y

% define where it starts where it ends
if idx == length(x)
  disp('extrapolating!!!')
      ys = y(idx-1); ye = y(idx);
      xs = x(idx-1); xe = x(idx);
else
  pc = isposrealscalar(diff(x(1:2)));
  if pc == 1
    if x(idx) >= x_y
	ys = y(idx-1); ye = y(idx);
	xs = x(idx-1); xe = x(idx);
	else
	ys = y(idx); ye = y(idx+1);
	xs = x(idx); xe = x(idx+1);
      end
  else
    if x(idx) >= x_y
      ys = y(idx+1); ye = y(idx);
      xs = x(idx+1); xe = x(idx);
      else
      ys = y(idx); ye = y(idx-1);
      xs = x(idx); xe = x(idx-1);
    end
  end
end

y_x = ys+(((x_y-xs)./(xe-xs))*(ye-ys)); % linear interpolation/extrapolation