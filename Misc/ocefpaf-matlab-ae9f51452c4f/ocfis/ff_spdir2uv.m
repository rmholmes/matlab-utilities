function [u,v] = ff_spdir2uv(dir,spd)
% ff_spdir2uv.m -> computes u, v components from speed and direction
%
% use:  [u,v] = ff_spdir2uv(dir,spd);
% input:
%            speed
%            direction [degrees]
%
% output:
%          [u,v] -> in speed units
%
% example:
%          [u,v] = ff_spdir2uv(45,14.1421);
%
% other m-files required: cart2pol.m

%
% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     21-Oct-2009
% modified: 21-Oct-2009
%
% obs: check also ff_uv2spdir.m
%

dir   = dir.*pi./180;

[u,v] = pol2cart(dir,spd);
