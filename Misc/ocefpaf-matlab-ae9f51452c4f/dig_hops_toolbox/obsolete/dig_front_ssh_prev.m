%
% dig_front_ssh_prev.m
% purpose: digitize the surface thermal front from a fermi SST image
%          capturing the lon,lat coordinates using the mouse
% author: Filipe Fernandes
% date: 16-Sep-08
%
% obs: same as front_dig_prev but with AVISO SSH data
%

 cls

 disp(' '), disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            disp('% Choose the image from GS - region %')
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
 disp(' ')

%% read and plot image using uigetfile
 [arqi,pafi] = uigetfile('gs*.png','Choose the image from GS - region');
 [xysat,mapsat,aplha] = imread([pafi,arqi],'png');
 lat0=45; lon0=-85; latf=20; lonf=-55; % [lon,lat] limits for GS region
 xysatinv=flipud(xysat); xysatcrop=xysatinv(130:910,80:869); xysatcrop=flipud(xysatcrop); % plot croped image

disp(' '), disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
           disp('% loading the same for SW - region  %')
           disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' ')

 [xysatSW,mapsatSW,aplhaSW] = imread([pafi,'SW',arqi(3:end)],'png');
 lat0SW=55; lon0SW=-80; latfSW=30; lonfSW=-45; % [lon,lat] limits for SW region
 xysatinvSW=flipud(xysatSW); xysatcropSW=xysatinvSW(130:1292,80:1279); xysatcropSW=flipud(xysatcropSW); % plot croped image

%% load previous file
 [arq,paf] = uigetfile('gsfm.*','Choose the previous file'); load([paf,arq]);

%% show the previous front
figure(1);set(1,'color','w'); hold on
 imagesc([lon0SW lonfSW],[lat0SW latfSW],ind2rgb(xysatcropSW,mapsatSW));
 imagesc([lon0 lonf],[lat0 latf],ind2rgb(xysatcrop,mapsat));
 axis xy; axis([-85 -45 20 55])
 imday=datenum([arqi(9:10),'-',arqi(6:8),'-20',arqi(4:5)])-datenum(['00-00-20',arqi(4:5)]); frday=arq(10:end-4);
 tit=['Correct the points with the mouse - When done press enter - image day: ', num2str(imday),' and front day: ', num2str(frday)];
 title(tit,'FontSize',14,'FontWeight','Bold')

%% feature model domain for NWA
 load ../HOPS_mat/model_bath.mat; contour(lonb,latb,depth,[2000 1000 500],'r'); % plot bathymetry
 plot([-82.5695,-59.2661],[40.9746,49.6855],'k'), plot([-59.2661,-50.8305],[49.6855,35.8254],'k')
 plot([-50.8305,-74.1339],[35.8254,27.1145],'k'), plot([-74.1339,-82.5695],[27.1145,40.9746],'k')

%% correct digitized data points Absolute SSH
 t=datenum([arqi(9:10),'-',arqi(6:8),'-',['20',arqi(4:5)]]);  % date from image filename
 bbox = [-85 -45 20 55]; bbox(1:2) = 360+bbox(1:2); % [lonmin lonmax latmin latmax] positive longitude for aviso server
 url='http://opendap.aviso.oceanobs.com/thredds/dodsC/dt_upd_global_merged_madt_h'; % UPD ABS SSH
 time=nc_varget(url,'time'); time = datenum(1950, 1, time, 0, 0, 0); % days since 1950-01-01
 NbLatitudes=nc_varget(url,'NbLatitudes'); NbLongitudes=nc_varget(url,'NbLongitudes');

%% subset
 [d,itime]   = min(abs(time-t));
 [d,ilonmin] = min(abs(NbLongitudes-bbox(1))); [d,ilonmax] = min(abs(NbLongitudes-bbox(2)));
 [d,ilatmin] = min(abs(NbLatitudes-bbox(3)));  [d,ilatmax] = min(abs(NbLatitudes-bbox(4)));
 Nblon=eval(['length(' sprintf('[%d:%d]', [ilonmin ilonmax]),')']);
 Nblat=eval(['length(' sprintf('[%d:%d]', [ilatmin ilatmax]),')']);

%% assign values
 ssh=nc_varget(url,'Grid_0001',[itime ilonmin ilatmin],[1 Nblon Nblat]); ssh = ssh';
 lon=nc_varget(url,'NbLongitudes',[ilonmin],[Nblon]);
 lat=nc_varget(url,'NbLatitudes',[ilatmin],[Nblat]);
 lon = [lon(find(lon<=180)); lon(find(lon>180))-360];
 [LON,LAT]=meshgrid(lon,lat); f = find(LON > -55 & LAT < 30); ssh(f)=NaN; % mask image blank limits
 [c,h]=contour(LON,LAT,ssh,'Color',[.25 .25 .25]); % plot ssh
 h=plot(-gsfm(:,3),gsfm(:,2),'o','MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','k'); moveplot(h,'axy'); % plot front

 clc, disp(' '), disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
                 disp('%    Correct the points with the mouse     %')
                 disp('%          When done press enter           %')
                 disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
      disp(' '),

 pause; moveplot(h,'off')
 xdata=get(h,'XData'); ydata=get(h,'YData');
 day = zeros(length(xdata),1)+datenum([arqi(9:10),'-',arqi(6:8),'-20',arqi(4:5)])-datenum(['00-00-20',arqi(4:5)]);

%% save gsfm.yyyyddd.dig file, day is year day
 aux=[day ydata' -xdata'];
 [m,n]=size(aux);
 str2='%8.4f '; str1='%6.0f ';
 fmt=[repmat(str1,1,1),repmat(str2,1,2)];
 fid1=fopen(['gsfm.20',arqi(4:5),num2str(day(1)),'.dig'],'w');
 fprintf(fid1,[fmt, '\n'],aux');
 fclose(fid1); clear fid1

%% display image used
 txt=[arqi(9:10),'-',arqi(6:8),'-',['20',arqi(4:5)]];
 clc, disp(' '), disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
                 disp('%    Last image used for digitalization    %')
                 disp( ['%               ',txt,'                %'] )
                 disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
      disp(' ')