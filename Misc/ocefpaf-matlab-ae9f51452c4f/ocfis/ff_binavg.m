function [ybin,xbin,n,s] = ff_binavg(x,y,db)
% ff_binavg.m - bin average y(x) onto ybin(xbin)
%
% use:   [ybin,xbin,n,s] = ff_binavg(x,y,db)
% input:
%    y(x) - data (vel/prof)
%    x    - data axis (time/depth)
%    db   - bin interval in x units
%
% output:
%    ybin - binned data
%    xbin - new binned axis
%    n    - number of points in each bin
%    s    - std of data in each bin
%
% example:
%    [ybin, n, s] = ff_binavg(time, uvel, 1)
%
% other m-files required: matlab's histc.m
%
% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     28-May-2009
% modified: 28-May-2009
%
% obs: modifyed from bindata.m
%

[yr,yc] = size(y);

x    = x(:);
y    = y(:);
idx  = find(~isnan(y));
xbin = ceil(x(1)):db:floor(x(end));

if (isempty(idx))
    b = nan * xbin;
    n = nan * xbin;
    s = nan * xbin;
end

x = x(idx);
y = y(idx);

binwidth = diff(xbin);
xx = [xbin(1)-binwidth(1)/2 xbin(1:end-1) + binwidth/2 xbin(end)+binwidth(end)/2];

% shift bins so the interval is "( ]" instead of "[ )".
bins = xx + max(eps,eps*abs(xx));

[n,bin] = histc(x,bins,1);
n = n(1:end-1); n(n==0)=NaN;

idx  = find(bin>0);
sum  = full(sparse(bin(idx),idx*0+1,y(idx)));
sum  = [sum;zeros(length(xbin)-length(sum),1)*nan]; % add zeroes to the end
ybin = sum./n;
sum  = full(sparse(bin(idx),idx*0+1,y(idx).^2));
sum  = [sum;zeros(length(xbin)-length(sum),1)*nan]; % add zeroes to the end
s    = sqrt(sum./(n-1) - ybin.^2.*n./(n-1) );