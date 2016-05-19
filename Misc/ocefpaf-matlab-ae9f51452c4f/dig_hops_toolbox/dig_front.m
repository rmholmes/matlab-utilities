%
% dig_front.m
% purpose: SST image thermal front [lon lat] from mouse input
% author: Filipe Fernandes
% date: 12-Oct-08
%
% obs: using GS and Storm watch together
%      lon: 80-40 W lat: 30-50 N
%      matlab scaling issue (use ginout2.m)
%      ssh option with automatic check added

cls, disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
     disp('% Choose the image from GS - region %')
     disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%% read and plot image using uigetfile
 [arqi,pafi] = uigetfile('../gs*.png','Choose the image from GS - region');
 [xysat,mapsat,aplha] = imread([pafi,arqi],'png');
 lat0=45; latf=20; lon0=-85; lonf=-55; % [lon,lat] limits
 xysatinv=flipud(xysat); xysatcrop=xysatinv(130:910,80:869); xysatcrop=flipud(xysatcrop); % plot croped image

clc, disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
     disp('% loading the same for SW - region  %')
     disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

 [xysatSW,mapsatSW,aplhaSW] = imread([pafi,'SW',arqi(3:end)],'png');
 lat0SW=55; latfSW=30; lon0SW=-80; lonfSW=-45;  % [lon,lat] limits
 xysatinvSW=flipud(xysatSW); xysatcropSW=xysatinvSW(130:1292,80:1279); xysatcropSW=flipud(xysatcropSW); % plot croped image

figure, set(gcf,'color','w'); hold on;
 imagesc([lon0SW lonfSW],[lat0SW latfSW],ind2rgb(xysatcropSW,mapsatSW));
 imagesc([lon0 lonf],[lat0 latf],ind2rgb(xysatcrop,mapsat));
 axis xy; axis([-85 -45 20 55])
 title('Digitize GS path with the LEFT button then click the RIGHT button to end','FontSize',20,'FontWeight','Bold')

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
 hlp = menu('Method for digitalization:','No previous front','From previous front','Previous front as a guide');
 pause(1) % bug, bug ...
 if hlp==2 | hlp==3; [arq,paf] = uigetfile('../gsfm*.dig','Choose the previous file'); end
 if hlp==2; load([paf,arq]); h=plot(-gsfm(:,3),gsfm(:,2),'k.'); end
 if hlp==3; load([paf,arq]); plot(-gsfm(:,3),gsfm(:,2),'r'); end

if hlp==1 | hlp==3;
	clc, disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
	     disp('% LEFT button digitize        %')
	     disp('% RIGHT button zoom in        %')
	     disp('% DOUBLE LEFT button zoom out %')
	     disp('% -- Press ENTER when done -- %')
	     disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
	%% calls ginput2
	[xdata,ydata] = ginput2('k.');
	%% plot the digitized front for corrections
	title('Correct the points with the mouse - When done press enter','FontSize',20,'FontWeight','Bold')
	h=plot(xdata,ydata,'o','MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','k'); % plot front
end

if hlp==2;
	imday=datenum([arqi(9:10),'-',arqi(6:8),'-20',arqi(4:5)])-datenum(['00-00-20',arqi(4:5)]); frday=arq(10:end-4);
	tit=['Correct the points with the mouse image day: ', num2str(imday),' and front day: ', num2str(frday)];
	title(tit,'FontSize',20,'FontWeight','Bold')
end

 moveplot(h,'axy'); % start move plot

clc, disp(' '), disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
                disp('%    Correct the points with the mouse     %')
                disp('%          When done press enter           %')
                disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
     disp(' '),

pause; moveplot(h,'off')
xdata=get(h,'XData'); ydata=get(h,'YData');

%% SSH date check
 url='http://opendap.aviso.oceanobs.com/thredds/dodsC/duacs_global_nrt_madt_merged_h'; % NRT MADT
 nrt=nc_varget(url,'time'); nrt = datenum(1950, 1, nrt(end), 0, 0, 0);
 url='http://opendap.aviso.oceanobs.com/thredds/dodsC/dt_upd_global_merged_madt_h'; % DT UPD MADT
 upd=nc_varget(url,'time'); upd = datenum(1950, 1, upd(end), 0, 0, 0);
 t=datenum([arqi(9:10),'-',arqi(6:8),'-',['20',arqi(4:5)]]);  % date from image filename

if t > nrt
	ssh_opt=0;
else
	ssh_opt=1;
	if t < upd
		clc, disp('Using Updated Absolute SSH')
	else
		url='http://opendap.aviso.oceanobs.com/thredds/dodsC/duacs_global_nrt_madt_merged_h'; % NRT MADT
		clc, disp('Using Near Real Time Absolute SSH')
	end
end

if ssh_opt == 1
	%% correct digitized data points with SSH
	bbox = [-85 -45 20 55]; bbox(1:2) = 360+bbox(1:2); % [lonmin lonmax latmin latmax] positive longitude for server
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
	h=plot(xdata,ydata,'r.'); moveplot(h,'axy'); pause; moveplot(h,'off')
	xdata=get(h,'XData'); ydata=get(h,'YData');
end

day=zeros(length(xdata),1)+datenum([arqi(9:10),'-',arqi(6:8),'-20',arqi(4:5)])-datenum(['00-00-20',arqi(4:5)]);

%% save gsfm.yyyyddd file, day is year day
 aux=[day ydata' -xdata'];
 [m,n]=size(aux);
 str2='%8.4f '; str1='%6.0f ';
 fmt=[repmat(str1,1,1),repmat(str2,1,2)];
 fid1=fopen(['gsfm.20',arqi(4:5),num2str(day(1)),'.dig'],'w');
 fprintf(fid1,[fmt, '\n'],aux');
 fclose(fid1); clear fid1