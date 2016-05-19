function [jhour] = oasp_jhour(datevec)
%
% oasp_jhour.m - description of what the function performs
%
% use:  [jhour] = oasp_jhour([yyyy mm dd hh mi se]);
% input:
%    input - yyyy  -> year
%          - mm    -> month
%          - dd    -> day
%          - hh    -> hour
%          - mi    -> minutes
%          - se    -> seconds
%
% output:
%    output - Julian hours for OASP at 0000 hours, 1 Jan, 1900
%
% example:
%    859296 = oasp_jhour([1998 01 11 00 00 00])

% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     04-Oct-2009
% modified: 04-Oct-2009
%
% obs: based on RPStuff julian.m
%

% Convert Gregorian calendar dates to corresponding
% Julian hours. Although the formal definition
% holds that Julian days start and end at noon, here
% Julian days start and end at midnight.
%
% In this convention, Julian day 0 began at 0000 hours, 1 Jan, 1900.
% 859296.00000000    jan 11 1998   0: 0:  0.000
% 2147483647
y = datevec(:,1);
m = datevec(:,2);
d = datevec(:,3);
h = datevec(:,4);
min = datevec(:,5);
sec = datevec(:,6);

h = h+(min+sec/60)/60;
mo = m+9;
yr = y-1;
i = find(m>2);
mo(i) = m(i)-3;
yr(i) = y(i);
c = floor(yr/100);
yr = yr - c*100;
j = floor((146097*c)/4) + floor((1461*yr)/4) + floor((153*mo +2)/5) + ...
              d -693902; %+1721119; here is the key date difference
jhour = j*24+h;
