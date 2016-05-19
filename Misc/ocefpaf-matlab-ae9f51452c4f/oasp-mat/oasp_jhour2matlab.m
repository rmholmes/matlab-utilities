function [matdate] = oasp_jhour2matlab(jhours)
%
% oasp_jhour2matlab.m - convert OASP jhours to matlab datenum
%
% use:  [matdate] = oasp_jhour2matlab(jhours);
% input:
%    input - jhours from OASP
%
% output:
%    output - matlab date serial date number
%
% example:
%    729766 = oasp_jhour2matlab(859296)
%
% other m-files required: oasp_greg.m

% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     05-Oct-2009
% modified: 05-Oct-2009
%
% obs: 
%

[gtime]   = oasp_greg(jhours);
matdate = datenum(gtime);