function map = cgray(m);
% ctopo.m - colormap for positive/negative data with grayscale only
%
% use:  map = cgray(m);
% input:
%    input - colormap length (64)
%
% output:
%    output - colormap
%
% example:
%    cm1 = cgray(255);
%
% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     22-May-2009
% modified: 22-May-2009
%
% obs: original from cushman-roisin book cd-rom
%

values  = [0.:(1-0.)*1/(m-1):1].'; 
map    = [values values values];