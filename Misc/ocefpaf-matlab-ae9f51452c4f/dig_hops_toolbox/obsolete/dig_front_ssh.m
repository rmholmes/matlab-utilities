%
% dig_front_ssh.m
% purpose: SST image thermal front [lon lat] from mouse input
% author: Filipe Fernandes
% date: 16-Sep-08
%
% obs: same as dig_front but with AVISO SSH data
%

cls

disp(' '), disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
           disp('% Choose the image from GS - region %')
           disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' '),

%% read and plot image using uigetfile
 [arqi,pafi] = uigetfile('../gs*.png','Choose the image from GS - region');
 [xysat,mapsat,aplha] = imread([pafi,arqi],'png');
 lat0=45; latf=20; lon0=-85; lonf=-55; % [lon,lat] limits
 xysatinv=flipud(xysat); xysatcrop=xysatinv(130:910,80:869); xysatcrop=flipud(xysatcrop); % plot croped image

disp(' '), disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
           disp('% loading the same for SW - region  %')
           disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' '),

 [xysatSW,mapsatSW,aplhaSW] = imread([pafi,'SW',arqi(3:end)],'png');
 lat0SW=55; latfSW=30; lon0SW=-80; lonfSW=-45;  % [lon,lat] limits
 xysatinvSW=flipud(xysatSW); xysatcropSW=xysatinvSW(130:1292,80:1279); xysatcropSW=flipud(xysatcropSW); % plot croped image

figure(1);set(1,'color','w'); hold on;
 imagesc([lon0SW lonfSW],[lat0SW latfSW],ind2rgb(xysatcropSW,mapsatSW));
 imagesc([lon0 lonf],[lat0 latf],ind2rgb(xysatcrop,mapsat));
 axis xy; axis([-85 -45 20 55])
 title('Digitize GS path with the LEFT button then click the RIGHT button to end','FontSize',14,'FontWeight','Bold')

%% feature model domain
 plot([-82.5695,-59.2661],[40.9746,49.6855],'k--'),  plot([-59.2661,-50.8305],[49.6855,35.8254],'k--')
 plot([-50.8305,-74.1339],[35.8254,27.1145],'k--'),  plot([-74.1339,-82.5695],[27.1145,40.9746],'k--')
%% numerical model domain
 tgrid2= nc_varget('../HOPS_files/grids_NWA_15km.nc','tgrid2');
 plot(tgrid2(1,:,1),tgrid2(1,:,2),'k'), plot(tgrid2(:,end,1),tgrid2(:,end,2),'k')
 plot(tgrid2(:,1,1),tgrid2(:,1,2),'k'), plot(tgrid2(end,:,1),tgrid2(end,:,2),'k')
%% plot bathymetry
 load ../HOPS_mat/model_bath.mat; contour(lonb,latb,depth,[2000 1000 500],'k');

%% digitize data points
clc, disp(' '), disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
                disp('% LEFT button digitize        %')
                disp('% RIGHT button zoom in        %')
                disp('% DOUBLE LEFT button zoom out %')
                disp('% -- Press ENTER when done -- %')
                disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
     disp(' ')

[xdata,ydata] = ginput2('k.');

%  %% show the clicked front
%  close(1); figure(1);set(1,'color','w'); hold on
%   imagesc([lon0SW lonfSW],[lat0SW latfSW],ind2rgb(xysatcropSW,mapsatSW));
%   imagesc([lon0 lonf],[lat0 latf],ind2rgb(xysatcrop,mapsat));
%   axis xy; axis([-85 -45 20 55])
 title('Correct the points with the mouse - When done press enter','FontSize',14,'FontWeight','Bold')

%% correct digitized data points Absolute SSH
 t=datenum([arqi(9:10),'-',arqi(6:8),'-',['20',arqi(4:5)]]);  % date from image filename
 bbox = [-85 -45 20 55]; bbox(1:2) = 360+bbox(1:2); % [lonmin lonmax latmin latmax] positive longitude for server
 url='http://opendap.aviso.oceanobs.com/thredds/dodsC/dt_upd_global_merged_madt_h'; % UPD ABS SSH
 time=nc_varget(url,'time'); time = datenum(1950, 1, time, 0, 0, 0);
 NbLatitudes=nc_varget(url,'NbLatitudes'); NbLongitudes=nc_varget(url,'NbLongitudes');
 [d,itime]   = min(abs(time-t)); % subset
 [d,ilonmin] = min(abs(NbLongitudes-bbox(1))); [d,ilonmax] = min(abs(NbLongitudes-bbox(2)));
 [d,ilatmin] = min(abs(NbLatitudes-bbox(3)));  [d,ilatmax] = min(abs(NbLatitudes-bbox(4)));
 Nblon=eval(['length(' sprintf('[%d:%d]', [ilonmin ilonmax]),')']);
 Nblat=eval(['length(' sprintf('[%d:%d]', [ilatmin ilatmax]),')']);

%% assign values
 ssh=nc_varget(url,'Grid_0001',[itime ilonmin ilatmin],[1 Nblon Nblat]); ssh = ssh';
 lon=nc_varget(url,'NbLongitudes',[ilonmin],[Nblon]);
 lat=nc_varget(url,'NbLatitudes',[ilatmin],[Nblat]);
 lon = [lon(find(lon<=180)); lon(find(lon>180))-360];
 [c,h]=contour(lon,lat,ssh,'Color',[.25 .25 .25]); % plot ssh
 load ../HOPS_mat/model_bath.mat; contour(lonb,latb,depth,[2000 1000 500],'r'); % plot bathymetry
 h=plot(xdata,ydata,'o','MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','k'); moveplot(h,'axy'); % plot front

clc, disp(' '), disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
                disp('%    Correct the points with the mouse     %')
                disp('%          When done press enter           %')
                disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
     disp(' ')

pause; moveplot(h,'off')
xdata=get(h,'XData'); ydata=get(h,'YData');
day=zeros(length(xdata),1)+datenum([arqi(9:10),'-',arqi(6:8),'-20',arqi(4:5)])-datenum(['00-00-20',arqi(4:5)]);

% save gsfm.yyyyddd.dig file, day is year day
 aux=[day ydata' -xdata'];
 [m,n]=size(aux);
 str2='%8.4f '; str1='%6.0f ';
 fmt=[repmat(str1,1,1),repmat(str2,1,2)];
 fid1=fopen(['gsfm.20',arqi(4:5),num2str(day(1)),'.dig'],'w');
 fprintf(fid1,[fmt, '\n'],aux');
 fclose(fid1); clear fid1