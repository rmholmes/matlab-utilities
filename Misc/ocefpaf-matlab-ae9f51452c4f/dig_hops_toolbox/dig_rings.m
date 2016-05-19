%
% dig_rings.m
% purpose: SST image thermal front [lon lat] from mouse input
% author: Filipe Fernandes
% date: 16-Sep-08
%
% obs:  uses plot_front.m improve this part...
%

cls

disp(' '), disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
           disp('% Choose the image from GS - region %')
           disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' ')

%% read and plot image using uigetfile
 [arqi,pafi] = uigetfile('../gs*.png','Choose the image from GS - region');
 [xysat,mapsat,aplha] = imread([pafi,arqi],'png');
 lat0=45; latf=20; lon0=-85; lonf=-55; % [lon,lat] limits
 xysatinv=flipud(xysat); xysatcrop=xysatinv(130:910,80:869); xysatcrop=flipud(xysatcrop); % crop image

disp(' '), disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
           disp('% loading the same for SW - region  %')
           disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' ')

 [xysatSW,mapsatSW,aplhaSW] = imread([pafi,'SW',arqi(3:end)],'png');
 lat0SW=55; latfSW=30; lon0SW=-80; lonfSW=-45;  % [lon,lat] limits
 xysatinvSW=flipud(xysatSW); xysatcropSW=xysatinvSW(130:1292,80:1279); xysatcropSW=flipud(xysatcropSW); % crop image

figure, set(gcf,'color','w'); hold on;
 imagesc([lon0SW lonfSW],[lat0SW latfSW],ind2rgb(xysatcropSW,mapsatSW));
 imagesc([lon0 lonf],[lat0 latf],ind2rgb(xysatcrop,mapsat));
 axis xy; axis([-85 -45 20 55])
 title('click CENTER then EDGE for each RING from LEFT to RIGHT','FontSize',20,'FontWeight','Bold')

%% feature model domain for NWA
 plot([-82.5695,-59.2661],[40.9746,49.6855],'k--'),  plot([-59.2661,-50.8305],[49.6855,35.8254],'k--')
 plot([-50.8305,-74.1339],[35.8254,27.1145],'k--'),  plot([-74.1339,-82.5695],[27.1145,40.9746],'k--')
%% plot front
 year=['20',arqi(4:5)]; day=datenum([arqi(9:10),'-',arqi(6:8),'-',year])-datenum(['00-00-',year]);
 frontfile=['../EXPS_',num2str(year),'/',num2str(day),'/gsfm.',num2str(year),num2str(day),'.dig'];
 plt_front(frontfile,'k','.',5);

%% numerical model domain
 tgrid2= nc_varget('../HOPS_files/grids_NWA_15km.nc','tgrid2');
 plot(tgrid2(1,:,1),tgrid2(1,:,2),'k'), plot(tgrid2(:,end,1),tgrid2(:,end,2),'k')
 plot(tgrid2(:,1,1),tgrid2(:,1,2),'k'), plot(tgrid2(end,:,1),tgrid2(end,:,2),'k')
%% plot bathymetry
 load ../HOPS_mat/model_bath.mat; contour(lonb,latb,depth,[2000 1000 500],'k');

hlp = menu('Plot previous ring conf. to help ?','yes','no'); pause(1);

if hlp == 1;
	[arq,paf] = uigetfile('../ring*.dig','Choose the previous file');
	plt_ring([paf,arq],'r','x')
end

clc, disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
     disp('%   Input the number of Warm Core Rings     %')
     disp('% and then click on their center and edge   %')
     disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

nw=input('Number of WCR = ');

if nw > 0;
 [xw,yw]=ginput2(nw*2,'k.');
 rw=sw_dist(yw,xw,'km'); rw=rw(1:2:end);rw80= rw-rw*.20; xw=-xw(1:2:end); yw=yw(1:2:end);
 for i=1:length(xw); h(i)=ff_range_ring(-xw(i),yw(i),rw(i),[-85 -45],[20 55]); set(h(i),'Color','k'); end; plot(-xw,yw,'kx')
end

clc, disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
     disp('%   Input the number of Cold Core Rings     %')
     disp('% and then click on their centre and border %')
     disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

nc=input('Number of CCR = ');

if nc > 0;
 [xc,yc]=ginput2(nc*2,'k.');
 rc=sw_dist(yc,xc,'km'); rc=rc(1:2:end);rc80= rc-rc*.20; xc=-xc(1:2:end); yc=yc(1:2:end);
 for i=1:length(xc); h(i)=ff_range_ring(-xc(i),yc(i),rc(i),[-85 -45],[20 55]); set(h(i),'Color','k'); end; plot(-xc,yc,'kx');
end

%% save ring.yyyyddd.dig file, day is year day
 fid2=fopen(['ring.20',arqi(4:5),num2str(day(1)),'.dig'],'w');
 fprintf(fid2,' %1.0f\n',nw+nc);

if nw > 0;
   nring = length(rw);
   for i = 1:nring;
	eval(['fprintf(fid2,''w%1.0f %3.0f  %1.0f  %2.2f %2.2f  %2.0f  %2.0f  %2.0f %4.0f\n'',[i day(1) -1 yw(i)  xw(i) 80   rw(i)    rw80(i)   1090]);'])
   end
end

if nc > 0;
   nring = length(rc);
   for i = 1:nring;
	eval(['fprintf(fid2,''c%1.0f %3.0f %1.0f  %2.2f  %2.2f %2.0f  %2.0f  %2.0f %4.0f\n'',[i day(1) 1 yc(i)  xc(i) 80   rc(i)    rc80(i)   1090]);'])
   end
end

fprintf(fid2,'\n');
fprintf(fid2,'****\n');
fprintf(fid2,' column  1: name of ring\n');
fprintf(fid2,' column  2: year day\n');
fprintf(fid2,' column  3: rotation direction\n');
fprintf(fid2,' column  4: location (latitude)\n');
fprintf(fid2,' column  5: location (west longitude)\n');
fprintf(fid2,' column  6: maximum velocity (best guess), no change ? \n');
fprintf(fid2,' column  7: radius\n');
fprintf(fid2,' column  8: 80%% of outside radius (where the maximum velocity is?)\n');
fprintf(fid2,' column  9: depth in meters , never change\n');
fclose(fid2); clear fid2