% 
% vi.m
% purpose:  opens vi inside matlab
% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     24-Oct-2009
% modified: 24-Oct-2009
%
% obs: 
%

function theResult = vi(flags)

if exist('flags')
    eval(['!vi ',flags])
else
    !vi
end