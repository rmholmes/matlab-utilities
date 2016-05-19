% 
% l.m
% purpose:  list with colors
% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     15-Oct-2004
% modified: 15-Oct-2004
%
% obs: 
%

function theResult = l(flags)

if exist('flags')
    eval(['!ls -ls -A --group-directories-first -B --color ',flags])
else
    !ls -ls -A --group-directories-first -B --color
end
