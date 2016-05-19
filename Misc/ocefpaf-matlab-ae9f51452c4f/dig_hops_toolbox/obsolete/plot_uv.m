%
% plot_uv.m
%
% purpose: load sat absolute ssh and plot over fermi sst image
% author: Filipe Fernandes
% date: 14-Jul-08
%
% obs:  Changed to SNCTools with Java
%       Plot over both Gulf Stream and Storm Watch regions
%       The altimeter products were produced by Ssalto/Duacs and distributed by Aviso with support from Cnes.
%

cls

%% read and plot image using uigetfile
 [arqi,pafi] = uigetfile('gs*.png','PNG Image');
 [xysat,mapsat,aplha] = imread([pafi,arqi],'png'); xysatinv=flipud(xysat);
 xysatcrop=xysatinv(130:910,80:869); xysatcrop=flipud(xysatcrop);
 [xysatSW,mapsatSW,aplhaSW] = imread([pafi,'SW',arqi(3:end)],'png');
 xysatinvSW=flipud(xysatSW); xysatcropSW=xysatinvSW(130:1292,80:1279); xysatcropSW=flipud(xysatcropSW);

%% geostrophic velocity from AVISO:
%% NRT-near real time, DT - Delayed, UPD - updated, MSLA - mean sea level anomaly, MADT - mean absolute sea level
%   url='http://opendap.aviso.oceanobs.com/thredds/dodsC/duacs_global_nrt_msla_merged_uv'; % NRT-MSLA velocity
%   url='http://opendap.aviso.oceanobs.com/thredds/dodsC/duacs_global_nrt_madt_merged_uv'; % NRT-MADT velocity
%   url='http://opendap.aviso.oceanobs.com/thredds/dodsC/dt_upd_global_merged_msla_uv';    % UPD-MSLA velocity
 url='http://opendap.aviso.oceanobs.com/thredds/dodsC/dt_upd_global_merged_madt_uv';    % UPD-MADT velocity

%% download data
 t=datenum([arqi(9:10),'-',arqi(6:8),'-',['20',arqi(4:5)]]); % date from image filename
 bbox = [-85 -45 20 55]; bbox(1:2) = 360+bbox(1:2); % [lonmin lonmax latmin latmax] positive longitude for server
 time=nc_varget(url,'time'); time = datenum(1950, 1, time, 0, 0, 0); % days since 1950-01-01
 NbLatitudes=nc_varget(url,'NbLatitudes'); NbLongitudes=nc_varget(url,'NbLongitudes'); % 0-360

%% subset
 [d,itime]   = min(abs(time-t));
 [d,ilonmin] = min(abs(NbLongitudes-bbox(1)));
 [d,ilonmax] = min(abs(NbLongitudes-bbox(2)));
 [d,ilatmin] = min(abs(NbLatitudes-bbox(3)));
 [d,ilatmax] = min(abs(NbLatitudes-bbox(4)));
 Nblon=eval(['length(' sprintf('[%d:%d]', [ilonmin ilonmax]),')']);
 Nblat=eval(['length(' sprintf('[%d:%d]', [ilatmin ilatmax]),')']);
 %% assign values
 u=nc_varget(url,'Grid_0001',[itime ilonmin ilatmin],[1 Nblon Nblat]); u=u';
 v=nc_varget(url,'Grid_0002',[itime ilonmin ilatmin],[1 Nblon Nblat]); v=v';
 lon=nc_varget(url,'NbLongitudes',[ilonmin],[Nblon]);
 lat=nc_varget(url,'NbLatitudes',[ilatmin],[Nblat]);
 lon = [lon(find(lon<=180)); lon(find(lon>180))-360];

%% mask image blank limits
 [LON,LAT]=meshgrid(lon,lat); f = find(LON > -55 & LAT < 30); u(f)=NaN; v(f)=NaN;
%   LON(f)=NaN; LAT(f)=NaN;

%% plot vectors over image (attetion to the inversion for the latitude in imagesc)
 figure(1);set(1,'color','w'); hold on; maximize(1);
  imagesc([-80 -45],[55 30],ind2rgb(xysatcropSW,mapsatSW));
  imagesc([-85 -55],[45 20],ind2rgb(xysatcrop,mapsat));
  axis xy; axis([-85 -45 20 55])
  s=0.006; quiver(LON,LAT,u*s,v*s,0,'k');