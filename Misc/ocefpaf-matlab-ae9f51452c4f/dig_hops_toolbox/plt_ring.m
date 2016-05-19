function [h,ringname,lon,lat,radius]=plt_ring(ringfile,color,marker)
%
% [h,ringname,lon,lat,radius]=plt_ring(ringfile,color,marker)
%
% use:
% [h,ringname,lon,lat,radius]=plt_rings('./EXPS_2008/245/ring.2008245.dig','k','x')
% lon,lat-> ring center
% radius -> radius of the ring
% h -> handle for color lengend

%
% plt_ring.m
% purpose: plot digitized ring file on SST image
% author: Filipe Fernandes
% date: 12-Oct-08
%
% obs: function version
%

warning off

fid = fopen(ringfile, 'r');
tline = fgetl(fid); x=str2num(tline);

for i=1:x
	tline = fgetl(fid);
	if ~ischar(tline), break, end
	ringname(i,:)=tline(1:2);
	ring(i,:)=sscanf(tline,'%*s %d %d %g %g %d %d %d %d',[1,inf]);
end

[a,b]=size(ring);
for i=1:a; h(i)=ff_range_ring(-ring(i,4),ring(i,3),ring(i,6),[-85 -45],[20 55]); set(h(i),'Color',color); end
plot(-ring(:,4),ring(:,3),'Color',color,'Marker',marker,'LineStyle','none');

lon=-ring(:,4); lat=ring(:,3); radius=ring(:,6); h=h(1);