% 
% lk.m
% purpose:  list by size
% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     07-Oct-2010
% modified: 07-Oct-2010
%
% obs: 
%

function theResult = lx(flags)

if exist('flags')
    eval(['!ls -lSrh --color ',flags])
else
    !ls -lSrh --color
end
