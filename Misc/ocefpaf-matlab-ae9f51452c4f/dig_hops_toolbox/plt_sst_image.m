%
% plt_sst_image.m
% purpose: Plot the fermi SST image with its respective front and rings
%
% author: Filipe Fernandes
% date: 12-Oct-08
%
% obs: Plot final figure ? axis ? ...
%

 cls, disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
      disp('% Choose the image from GS - region %')
      disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%% read and plot image using uigetfile
 [arqi,pafi] = uigetfile('../gs*.png','Choose the image from GS - region');
 [xysat,mapsat,aplha] = imread([pafi,arqi],'png');
 lat0=45; lon0=-85; latf=20; lonf=-55; % [lon,lat] limits for GS region
 xysatinv=flipud(xysat); xysatcrop=xysatinv(130:910,80:869); xysatcrop=flipud(xysatcrop); % plot croped image

clc, disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
     disp('% loading the same for SW - region  %')
     disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')


 [xysatSW,mapsatSW,aplhaSW] = imread([pafi,'SW',arqi(3:end)],'png');
 lat0SW=55; lon0SW=-80; latfSW=30; lonfSW=-45; % [lon,lat] limits for SW region
 xysatinvSW=flipud(xysatSW); xysatcropSW=xysatinvSW(130:1292,80:1279); xysatcropSW=flipud(xysatcropSW); % plot croped image

%% name helpers
 year=['20',arqi(4:5)]; day=datenum([arqi(9:10),'-',arqi(6:8),'-',year])-datenum(['00-00-',year]);

%% plot everthing
figure, set(gcf,'color','w'); hold on
 imagesc([lon0SW lonfSW],[lat0SW latfSW],ind2rgb(xysatcropSW,mapsatSW));
 imagesc([lon0 lonf],[lat0 latf],ind2rgb(xysatcrop,mapsat));
 axis xy; axis([-85 -45 20 55])
 tit=['image day: ', num2str(day)]; title(tit,'FontSize',20,'FontWeight','Bold')
 set(gca,'FontSize',20)

%% numerical model domain
 tgrid2= nc_varget('../HOPS_files/grids_NWA_15km.nc','tgrid2');
 plot(tgrid2(1,:,1),tgrid2(1,:,2),'k--'), plot(tgrid2(:,end,1),tgrid2(:,end,2),'k--')
 plot(tgrid2(:,1,1),tgrid2(:,1,2),'k--'), plot(tgrid2(end,:,1),tgrid2(end,:,2),'k--')
%% plot bathymetry
 load ../HOPS_mat/model_bath.mat; contour(lonb,latb,depth,[2000 1000 500],'k');

%% plot front
 frontfile=['../EXPS_',num2str(year),'/',num2str(day),'/gsfm.',num2str(year),num2str(day),'.dig'];
 plt_front(frontfile,'k','.',5);

%% plot rings
 ringfile=['../EXPS_',num2str(year),'/',num2str(day),'/ring.',num2str(year),num2str(day),'.dig'];
 plt_ring(ringfile,'k','x');

%% save image %%
 orient(gcf,'landscape')
 fillPage(gcf, 'papersize', 'usletter', 'margins', [1 1 1 1]);
 axis([-80 -55 30 50])
%   axis([-75 -55 30 45])
%   eval(['print -dpng ',num2str(year),num2str(day),'.png'])