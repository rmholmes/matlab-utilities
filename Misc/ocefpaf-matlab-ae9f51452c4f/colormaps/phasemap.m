function map = phasemap(m);
% ctopo.m - colormap periodic/circula data (phase)
%
% use:  map = phasemap(m);
% input:
%    input - colormap length (64)
%
% output:
%    output - colormap
%
% example:
%    cm1 = phasemap(255);
%
% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     22-May-2009
% modified: 22-May-2009
%
% obs: original from (?)
%

if nargin < 1, m = size(get(gcf,'colormap'),1); end

theta = 2*pi*(0:m-1)/m;
circ = exp(i*theta');

% vertices of colour triangle
vred = -2; vgreen = 1-sqrt(3)*i; vblue = 1+sqrt(3)*i;
vredc = vred-circ; vgreenc = vgreen-circ; vbluec = vblue-circ;

red   = abs(imag(vgreenc.*conj(vbluec)));
green = abs(imag(vbluec.*conj(vredc)));
blue  = abs(imag(vredc.*conj(vgreenc)));

map = 1.5*[red green blue]./abs(imag((vred-vgreen)*conj(vred-vblue)));