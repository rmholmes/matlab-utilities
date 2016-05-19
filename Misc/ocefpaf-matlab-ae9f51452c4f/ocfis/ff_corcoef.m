function function [r,cu,cd] = ff_corcoef(x,y,alpha)
% ff_corcoef.m -> compute correlation coefficient
%
%   call:  [r,cu,cd] = ff_corcoef(x,y,alpha);
%
%  input:  x,y   -> data sets x and y
%          alpha -> significance level
%
% output:  r -> Correlation Coefficient rho
%          cu,cd -> upper and lower confidence interval
%
% example:
%    x = [1.0 2.3 3.1 4.8 5.6 6.3];
%    y = [2.6 2.8 3.1 4.7 4.1 5.3];
%    [r,cu,cd] = corrcoef (x,y,0.01)
%
% subfunctions: ff_cov.m

% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     27-Mar-2009
% modified: 21-Nov-2009
%
% obs: corrcoef(x,y) = cov(x,y)/(std(x)*std(y))
%

c = ff_cov(x,y);   % uses ff_cov to compute covariance
s = std(x)*std(y); % standart deviations
r = c./s;          % Pearson's r or coefficient of correlation

%% confidence interval
n  = length(x);                               % n -> length of the series
%  dz = (-erfinv(alpha-1)).*sqrt(2)./sqrt(n-3);  % compute delta z (matlab)
q = 1-alpha;
dz = 1./sqrt(n-3)*icdf('norm',((q+1)/2),0,1); % compute delta z (test)
z  = 0.5*log((1+r)./(1-r));                   % fisher's z-transform of r
cu = tanh(z+dz); cd = tanh(z-dz);             % upper and lower limits for r

%  PHI^{-1}((q+1)/2) is the half-width of the confidence interval
%  measured in standard deviations, where q is the confidence expressed
%  as a decimal (e.g., q = 0.95 for 95% confidence)