function [h] = ff_stickplt(t,u,v)
% ff_stickplt.m - plot a "nice" stickplot by scaling time with data and adding some NaN to "stretch"
%
% use:  ff_stickplt(t,u,v);
%
% input:
%    t   -> time vector
%    u,v -> velocity components
%
% output:
%   h -> figure handle
%
% example:
%  date = 1:100;
%  u    = [filter(hamming(20),1,rand(1,100)-.5)];
%  v    = [filter(hamming(20),1,rand(1,100)-.5)];
%
%  ff_stickplt(date,u,v)
%

% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     27-Oct-2009
% modified: 27-Oct-2009
%
% obs: it is almost the same as stickplot.m from oceans (whoi) toolbox,
%      except that advanced graphic handling (multiple plots, yoffset, data apsect ratio ...) 
%      is left for matlab. This avoid breaking with each new releases and give the user more options
%

%% column vectors
t = t(:); u = u(:); v = v(:);

maxmag = max( max( sqrt(u.^2 + v.^2) ) ); % maximun magnitude
xax    = t(end)-t(1); yax = xax;          % x-length
scale  = xax/maxmag/8;                    % scale

% scaled u,v by time
wu = u*scale;

xp = [ t  t+wu  NaN*ones(size(t)) ]'; xp=xp(:);
%  yp = [ yy wv+yy NaN*ones(size(t)) ]'; yp=yp(:);
yp = [ zeros(size(t)) v NaN*ones(size(t)) ]'; yp=yp(:);

% Now plot
h=plot(xp,yp,'-');