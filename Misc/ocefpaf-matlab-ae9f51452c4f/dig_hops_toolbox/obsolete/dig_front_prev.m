%
% dig_front_prev.m
% purpose: digitize the surface thermal front from a fermi SST image
%          capturing the lon,lat coordinates using the mouse
% author: Filipe Fernandes
% date: 16-Sep-08
%
% obs: only front update, for rings just re-run dig_rings.m
%

 cls

 disp(' '), disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            disp('% Choose the image from GS - region %')
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
 disp(' ')

%% read and plot image using uigetfile
 [arqi,pafi] = uigetfile('../gs*.png','Choose the image from GS - region');
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
 [arq,paf] = uigetfile('../gsfm*.dig','Choose the previous file'); load([paf,arq]);

%% show the previous front
figure(1);set(1,'color','w'); hold on
 imagesc([lon0SW lonfSW],[lat0SW latfSW],ind2rgb(xysatcropSW,mapsatSW));
 imagesc([lon0 lonf],[lat0 latf],ind2rgb(xysatcrop,mapsat));
 axis xy; axis([-85 -45 20 55])
 imday=datenum([arqi(9:10),'-',arqi(6:8),'-20',arqi(4:5)])-datenum(['00-00-20',arqi(4:5)]); frday=arq(10:end-4);
 tit=['Correct the points with the mouse - When done press enter - image day: ', num2str(imday),' and front day: ', num2str(frday)];
 title(tit,'FontSize',14,'FontWeight','Bold')

%% feature model domain
 plot([-82.5695,-59.2661],[40.9746,49.6855],'k--'),  plot([-59.2661,-50.8305],[49.6855,35.8254],'k--')
 plot([-50.8305,-74.1339],[35.8254,27.1145],'k--'),  plot([-74.1339,-82.5695],[27.1145,40.9746],'k--')
%% numerical model domain
 tgrid2= nc_varget('../HOPS_files/grids_NWA_15km.nc','tgrid2');
 plot(tgrid2(1,:,1),tgrid2(1,:,2),'k'), plot(tgrid2(:,end,1),tgrid2(:,end,2),'k')
 plot(tgrid2(:,1,1),tgrid2(:,1,2),'k'), plot(tgrid2(end,:,1),tgrid2(end,:,2),'k')
%% plot bathymetry
 load ../HOPS_mat/model_bath.mat; contour(lonb,latb,depth,[2000 1000 500],'k');

 h=plot(-gsfm(:,3),gsfm(:,2),'o','MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','k'); moveplot(h,'axy'); % plot front

 clc, disp(' '), disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
                 disp('%    Correct the points with the mouse     %')
                 disp('%          When done press enter           %')
                 disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
      disp(' ')

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

clc,  disp(' '), disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
                 disp('%    Last image used for digitalization    %')
                 disp( ['%               ',txt,'                %'] )
                 disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
      disp(' ')