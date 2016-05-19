function [c] = ff_cov(x,y)
% ff_cov.m -> Compute covariance
%
%   call:  [x] = ff_cov(x,y);
%
%  input:  x,y -> data sets x and y
%
% output:  c -> covariance of x and y
%
% example: [x] = ff_cov(x,y);
%

% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     27-Mar-2009
% modified: 21-Nov-2009
%
% obs: unbiased only 1/(n-1)
%

x = x(:); % make x and y one column
y = y(:);
n = length(x); % N -> length of the series
x = x-mean(x); % remove mean
y = y-mean(y);
c = (x'*y)/(n-1); % compute covariance