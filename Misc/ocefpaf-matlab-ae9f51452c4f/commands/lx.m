% 
% lx.m
% purpose:  list by extension
% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     15-Oct-2004
% modified: 15-Oct-2004
%
% obs: 
%

function theResult = lx(flags)

if exist('flags')
    eval(['!ls -lhBpX --color ',flags])
else
    !ls -lhBpX --color
end