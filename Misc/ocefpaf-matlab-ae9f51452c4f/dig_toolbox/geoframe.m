function [lon,lat] = geoframe
% geoframe.m -> visual selection of lon,lat box from a world coastline map
%
% output:  latitute/longitude box in degrees
%
% example: [lon,lat] = geoframe;
%

% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     18-Oct-2008
% modified: 28-Dec-2009
%
% obs: 
%

load coast
plot(long,lat,'k'), axis tight
questdlg('Select the rectangle than press enter','','Continue','Continue');

h = imrect(gca,[]);

pausa

api = iptgetapi(h);
box = api.getPosition(); api.delete();
lon = [box(1) box(1)+box(3)];
lat = [box(2) box(2)+box(4)];

close(gcf)