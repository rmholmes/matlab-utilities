function [spd,dir] = ff_uv2spdir(u,v,mag,rot)
% ff_uv2spdir.m -> computes speed and direction from u, v components
%                  allows for rotation and magnetic declination correction
%
% use:  [spd,dir] = ff_uv2spdir(u,v,mag,rot)
% input:
%            u,v -> velocity components
%            mag -> magnetic declination [deg], can be omitted
%            rot -> rotation [deg], can be omitted
%
% output:
%             spd -> speed
%             dir -> direction [deg]
%
% example:
%    [spd,dir] = ff_uv2spdir(10,10)
%
% other m-files required: cart2pol.m

%
% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     15-Oct-2009
% modified: 15-Oct-2009
%
% obs: need some tests for rotation and magnetic declination
%      check also spddir2uv.m

if nargin < 4,
  rot = 0;
  if nargin == 2,
    mag = 0;
  end 
end

% convert from cartesian to polar
[th,spd] = cart2pol(u,v);

% 90 is east and 270 is west
fz     = find(th < 0);
th(fz) = th(fz) + 2*pi;

% convert to degrees
dir = th.*(180./pi);

% apply rotation and magnetic correction
dir = dir - mag + rot;

% old version
%  vec = u+i*v;
%  spd = abs(vec);
%  dir = angle(vec);
%  dir = dir*180/pi;
%  dir = dir-mag + rot;
%  dir = mod(90-dir,360);