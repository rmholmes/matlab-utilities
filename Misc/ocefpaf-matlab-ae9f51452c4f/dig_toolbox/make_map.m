function [h] = make_map(file)
% make_map.m -> plota mapa usando seleção por mouse para área de plotagem
%               e arquivo de texto com pontos lon,lat
%
%   call:  [h] = make_map(file);
%
%  input:  file -> arquivo com pontos de longitude e latitude
%
% output:  h -> figure handle para modificações
%
% example: [h] = make_map('lonlat.txt');
%
% m-files required: geoframe.m
%        MAT-files: coast.mat
% 

% 
% make_map.m
% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     27-Mar-2010
% modified: 27-Mar-2010
%
% obs: 
%

[lonbb,latbb] = geoframe;

%  disp('Logitude limits:')
%  disp(lonbb)
%  
%  disp('Latitude limits:')
%  disp(lonbb)

data = str2num( cell2mat( readline(file,[2 -1]') ) );
lon = data(:,1);
lat = data(:,2);

% load gshhs  fine resolution coastline
coast = '~/svn-tools/mymatlab/datasets/gshhs/gshhs_f.b';
coast = gshhs(coast,latbb,lonbb);

lonc = getfield(coast,'Lon');
latc = getfield(coast,'Lat');

figure; set(gca,'color','none'); hold on
  plot(lonc, latc, 'k');
  h = plot(lon,lat,'ro','MarkerEdgeColor','b','MarkerFaceColor','r','MarkerSize',7);
  axis([lonbb latbb])


%  url   = 'http://blackburn.whoi.edu:8081/thredds/dodsC/bathy/etopo2_v2c.nc';
%  [z,g] = nj_subsetGrid(url,'topo',[lonlim(1) -68.5 latlim]);
%  [c,h] = contour(g.lon,g.lat,double(z),[-50 -50],'k');
%  [c,h] = contour(g.lon,g.lat,double(z),[-200 -200],'k--');