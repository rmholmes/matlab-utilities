
%This function plots the PV and diabatic/frictional PV terms along
%float tracks in eq40. Flt is the float number to display or float
%numbers to average over.

tic
load('Floats.mat');
fn = Floats3;
cflt = fn(50);
bad = [16 17 51 52 53 66]; %bad with 5-point horizontal
fn(bad) = [];
bad2 = [5 15 25 48]; %bad with 3-point vert.
fn(bad2) = [];

%Alternatively use choose floats:
%fn = Choose_floats(eq40r5F,cflt,[282 285],100000,15);

%Filenames:
%hname = eq40r5V;
hname = '/mnt/Data1/ryan/UW_ROMS_EQ40/run6k1/eq40run5k1VARS.nc';
%dname = eq40r5D;
dname = '/mnt/Data1/ryan/UW_ROMS_EQ40/run6k1/eq40run5k1DIAS.nc';
pname = eq40r5PV;
%fname = eq40r5F;
%fname = '/mnt/Data1/ryan/UW_ROMS_EQ40/run6k1/ocean_flt.nc';

custom_tlims = [282.0417 286.0417];


%History file properties:
ncid = netcdf.open(hname,'NC_NOWRITE');
Htime = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time'),'double');
[temp, HtimeL] = netcdf.inqDim(ncid,netcdf.inqDimID(ncid,'timeD'));
zvec = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'z'),'double');
netcdf.close(ncid);

%%%%%%%%%%%Flt file IDs:
ncidFLT = netcdf.open(fname,'NC_NOWRITE');
ftimeID = netcdf.inqVarID(ncidFLT,'ocean_time');
ftime = netcdf.getVar(ncidFLT,ftimeID,'double');
XgridFID = netcdf.inqVarID(ncidFLT,'Xgrid');
YgridFID = netcdf.inqVarID(ncidFLT,'Ygrid');
depthFID = netcdf.inqVarID(ncidFLT,'depth');
lonFID = netcdf.inqVarID(ncidFLT,'lon');
latFID = netcdf.inqVarID(ncidFLT,'lat');
rhoFID = netcdf.inqVarID(ncidFLT,'rho');
% $$$ ZgridFID = netcdf.inqVarID(ncidFLT,'Zgrid');
% $$$ tempFID = netcdf.inqVarID(ncidFLT,'temp');
% $$$ saltFID = netcdf.inqVarID(ncidFLT,'salt');
toc
tic

[temp HtimeindI] = min(abs(Htime/86400-custom_tlims(1)));
[temp HtimeindF] = min(abs(Htime/86400-custom_tlims(2)));
[temp mint] = min(abs(ftime-Htime(HtimeindI)));
[temp maxt] = min(abs(ftime-Htime(HtimeindF)));

if (mint > maxt)
    temp = mint;
    mint = maxt;
    maxt = temp;
end
flt_pos = zeros(length(fn),maxt-mint+1,7);
for f = 1:length(fn)
    fnn = fn(f);
    flt_pos(f,:,1) = ftime(mint:maxt)'/60/60/24;
    flt_pos(f,:,2) = netcdf.getVar(ncidFLT,XgridFID,[fnn-1 mint-1],[1 ...
                        maxt-mint+1],'double')'+1;
    flt_pos(f,:,3) = netcdf.getVar(ncidFLT,YgridFID,[fnn-1 mint-1],[1 ...
                        maxt-mint+1],'double')'+1;
    flt_pos(f,:,4) = netcdf.getVar(ncidFLT,depthFID,[fnn-1 mint-1],[1 ...
                        maxt-mint+1],'double')';
    flt_pos(f,:,5) = netcdf.getVar(ncidFLT,lonFID,[fnn-1 mint-1],[1 ...
                        maxt-mint+1],'double')';
    flt_pos(f,:,6) = netcdf.getVar(ncidFLT,latFID,[fnn-1 mint-1],[1 ...
                        maxt-mint+1],'double')';
     flt_pos(f,:,7) = netcdf.getVar(ncidFLT,rhoFID,[fnn-1 mint-1],[1 ...
                        maxt-mint+1],'double')';
% $$$     flt_pos(f,:,5) = abs(netcdf.getVar(ncidFLT,tempFID,[fnn-1 mint- ...
% $$$                         1],[1 maxt-mint+1],'double'))';
% $$$     flt_pos(f,:,6) = abs(netcdf.getVar(ncidFLT,saltFID,[fnn-1 mint- ...
% $$$                         1],[1 maxt-mint+1],'double'))';
end
netcdf.close(ncidFLT);
flt_pos(flt_pos>1e20) = NaN;
%Make bad values NaN's:

%Restrict time period if some of these floats have invalid values:
swi = 0;
mintaft = mint;
maxtaft = maxt;
for t = 1:length(flt_pos(1,:,1))
    if (swi == 0)
    if (length(find(isnan(flt_pos(:,t,2)))>0))
        mintaft = mint+t;
        flt_pos(:,t,:) = NaN;
    else
        swi = 1;
    end
    else
        if (length(find(isnan(flt_pos(:,t,2)))>0))
            maxtaft = mint+t-2;
            flt_pos(:,t:end,:) = NaN;
            break;
        end
    end
end
flt_pos = flt_pos(:,(mintaft-mint+1):(maxtaft-mint+1),:);
toc
mint = mintaft;
maxt = maxtaft;

%Get Domain limits:
ZGridPos = interp1(zvec,1:length(zvec),flt_pos(:,:,4));
XGridPos = flt_pos(:,:,2);
YGridPos = flt_pos(:,:,3);
minZ = floor(min(min(ZGridPos(ZGridPos<1000 & ZGridPos>0))))-1;
maxZ = length(zvec);
minX = floor(min(min(XGridPos(XGridPos<1000 & XGridPos>0))))-10;
maxX = ceil(max(max(XGridPos(XGridPos<1000 & XGridPos>0))))+10;
minY = floor(min(min(YGridPos(YGridPos<1000 & YGridPos>0))))-10;
maxY = ceil(max(max(YGridPos(YGridPos<1000 & YGridPos>0))))+10;
%Temporary for w'b':
% $$$ minX = floor(min(min(XGridPos(XGridPos<1000 & XGridPos>0))))-35;
% $$$ maxX = ceil(max(max(XGridPos(XGridPos<1000 & XGridPos>0))))+35;
% $$$ minY = floor(min(min(YGridPos(YGridPos<1000 & YGridPos>0))))-35;
% $$$ maxY = ceil(max(max(YGridPos(YGridPos<1000 & YGridPos>0))))+35;

[tmp minT] = min(abs(Htime-ftime(mint)));
[tmp maxT] = min(abs(Htime-ftime(maxt)));
if (minT > maxT)
    temp = minT;
    minT = maxT;
    maxT = temp;
end
domain = [minX minY minZ minT; maxX maxY maxZ maxT];
TimeVector = Htime(minT:maxT)/86400;
tL = length(TimeVector);
fL = length(flt_pos(:,1,1));
%Filtering for all:
% $$$ h_filter = [11 5];
% $$$ v_filter = 5;
% $$$ t_filter = 0;

h_filter = [5 0];
v_filter = 3;
t_filter = 0;

%%H slice PV isopycnals%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%For HCL:
FGF = {'u' 'v' 'w' 'DBDX' 'DBDY' 'N2' 'b_diff' 'something'};
HCL = FGF;
HCL{end} = ['Signlog(-2*(4).*(cdiff(x_rho,(1),''x'').*(4)+cdiff(x_rho,' ...
               '(2),''x'').*(5))-2*(5).*(cdiff(y_rho,(1),''y'').*(4)+cdiff(y_rho,(2),''y'').*(5)),-20,-16)'];
DVADV = FGF;
DVADV{end} = ['Signlog(-2*(4).*(cdiff(x_rho,(3),''x'').*(6))-2*(5).*(cdiff(y_rho,(3),''y'').*(6))+0*(1),-20,-16)'];
HCL{end} = ['-2*(4).*(cdiff(x_rho,(1),''x'').*(4)+cdiff(x_rho,' ...
               '(2),''x'').*(5))-2*(5).*(cdiff(y_rho,(1),''y'').*(4)+cdiff(y_rho,(2),''y'').*(5))'];
DVADV = FGF;
DVADV{end} = ['-2*(4).*(cdiff(x_rho,(3),''x'').*(6))-2*(5).*(cdiff(y_rho,(3),''y'').*(6))+0*(1)'];

FIELD = DVADV;

subplotsh = 1:4;
positions = zeros(4,4);
fig1 = figure
% $$$ set(fig1, 'Position', [112 499 1722 443]);%get(0,'Screensize'));
% $$$ positions(:,3) = 0.2;
% $$$ positions(:,4) = 0.75;
% $$$ positions(:,2) = 0.14;
% $$$ positions(1,1) = 0.05;
% $$$ positions(2,1) = 0.28;
% $$$ positions(3,1) = 0.51;
% $$$ positions(4,1) = 0.74;
% $$$ positions(4,3) = 0.23;
set(fig1, 'Position', [565 3 1276 971]);%get(0,'Screensize'));
positions(:,3) = 0.2039;
positions(:,4) = 0.2420;
positions(:,2) = 0.7158;
positions(1,1) = 0.05;
positions(2,1) = 0.28;
positions(3,1) = 0.51;
positions(4,1) = 0.74;
positions(4,3) = 0.23;
avgden = mean(mean(flt_pos(:,:,7)));
ptimes = [282:(1+1/3):286]*86400;
for i=1:4
timep = ptimes(i);
[temp indp] = min(abs(flt_pos(1,:,1)-timep/86400));
[temp indpH] = min(abs(Htime-timep));
subplot(3,4,subplotsh(i));
pcolPlotNC(hname,0,{'PVv' 'PVh' '(1)+(2)'},[avgden indpH],1);
xlabel('Longitude ($^\circ E$)','Interpreter','latex','FontSize',20);
ylabel('');
title('');
colorbar('off');
set(gca,'YTick',[]);
set(gca,'Xtick',[-142:2:-130]);
set(gca,'pos',positions(i,:));
caxis([-2 10]*1e-9);
%caxis([-1 1]);
axis([-142 -132.5 1 8]);
if (i == 1)
    ylabel('Latitude ($^\circ N$)','Interpreter','latex','FontSize',20);
    set(gca,'YTick',[0:1:9]);
elseif (i == 4 )
    colorbar;
% $$$ cb = colorbar;
% $$$ set(cb,'YTick',[-1:0.5:1]);
% $$$ set(cb,'YTickLabel',{'<-1e-17','-1e-19','||<1e-21','1e-19','>1e-17'});

end
hold on;
plot(flt_pos(:,indp,5),flt_pos(:,indp,6),'.m','MarkerSize',25);
hold off;
end


%%%%%%%%%Pos's,mom's%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CF = pi*6371000/180; 
posx = CF*(flt_pos(:,:,5)+132).*cos(pi/180*flt_pos(:,:,6));
posy = CF*flt_pos(:,:,6);
posz = flt_pos(:,:,4);
COM = [mean(posx,1)' mean(posy,1)' mean(posz,1)'];
hang = atan2(posy-repmat(COM(:,2),[1 fL])',posx-repmat(COM(:,1),[1 ...
                    fL])');
hrad = sqrt((posy-repmat(COM(:,2),[1 fL])').^2+(posx-repmat(COM(:,1),[1 ...
                    fL])').^2);
% $$$ 
% $$$ tvec = 282:4/10:286;
% $$$ tposvec = zeros(length(tvec),1);
% $$$ for t=1:length(tvec)
% $$$ [tmp tposvec(t)] = min(abs(flt_pos(1,:,1)'-tvec(t)));
% $$$ end

tposvec = 2:5:tL;
% $$$ %SD at each radius:
% $$$ angwedges = -pi:pi/2:pi;
% $$$ meanang = (angwedges(2:end)+angwedges(1:(end-1)))/2;
% $$$ 
% $$$ sdrad = zeros(length(angwedges)-1,length(tposvec));
% $$$ for t=1:length(tposvec)
% $$$     rad = hrad(:,tposvec(t));
% $$$     ang = hang(:,tposvec(t));
% $$$ for a = 1:(length(angwedges)-1)
% $$$     sdrad(a,t) = 3*std(rad(ang>angwedges(a) & ang<angwedges(a+1)));
% $$$ end
% $$$ end

separator = 1.2e5;
cenpos = -separator;
%figure;

subplot(3,4,[5 6 7 8]);
for t = 1:length(tposvec)
    cenpos = cenpos+separator;
    plot(hrad(:,tposvec(t)).*cos(hang(:,tposvec(t)))+cenpos,...
         hrad(:,tposvec(t)).*sin(hang(:,tposvec(t))),...
         '.m','MarkerSize',15);
    hold on;
    plot(cenpos,0,'.k','MarkerSize',25);
% $$$     plot([sdrad(:,t).*cos(meanang'); sdrad(1,t).*cos(meanang(1));]+cenpos,[sdrad(:,t).*sin(meanang'); sdrad(1,t).*sin(meanang(1));],'-r', ...
% $$$          'LineWidth',2);
end
daspect([1 1 1]);
axis([-separator/2 length(tposvec)*separator-separator/2 -separator/2 ...
      separator/2]);
xlabel('Day');
ylabel('H distr. (km)');
set(gca,'YTick',[-5 0 5]*1e4);
set(gca,'YTickLabel',{'-50','0','50'});
% $$$ set(gca,'XTick',[0:separator:((length(tposvec)-1)*separator)]);
% $$$ LabelList = cell(length(tposvec),1);
% $$$ for t = 1:length(tposvec)
% $$$     LabelList{t} = num2str(TimeVector(tposvec(t)),4);
% $$$ end

set(gca,'XTick',[0 separator 2*separator-50000 2*separator ...
                 2*separator+50000 3*separator:separator:((length(tposvec)-1)*separator)]);
LabelList = cell(length(tposvec)+2,1);
for t = 1:length(tposvec)
    if (t>3)
        shf = 2;
    elseif (t == 3)
        shf = 1;
    else
        shf = 0;
    end
    LabelList{t+shf} = num2str(TimeVector(tposvec(t)),4);
end

set(gca,'XTickLabel',LabelList);
%set(gca,'XTick',[]);         
set(gca,'Position',[0.0588 0.4202 0.8707 0.2278]);


subplot(3,4,[9 10 11 12]);
FltProfile(flt_pos,hname,dname,{'PVh'},domain,[5 0],0,tposvec);
caxis([-4 4]*1e-9);
%caxis([-1 1]);
text(285.6,-130,'$q_h (s^{-3})$','Interpreter','latex','FontSize',30, ...
     'Color','k');
% $$$ text(285.4,-130,'$HCL (s^{-5})$','Interpreter','latex','FontSize',30, ...
% $$$      'Color','k');
% $$$ text(285.2,-130,'$DVADV (s^{-5})$','Interpreter','latex','FontSize',30, ...
% $$$      'Color','k');
set(gca,'Position',[0.0862 0.0865 0.8331 0.3007]);
set(gca,'XTick',TimeVector(tposvec)); 
set(gca,'XTickLabel',{LabelList{[1 2 4 6:length(LabelList)]}});

%%%%%%%%%Two Panel HCL/DVADV along float
%%%%%%%%%track%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HCLf_h = [5 0];
HCLf_h = [11 0];
HCLf_v = 0;

figure;
set(gcf,'Position',[165 239 1436 652]);

subplot(2,1,1);
[T,D,VarP] = FltProfile(flt_pos,hname,dname,HCL,domain,HCLf_h,HCLf_v, ...
                        tposvec);
h = pcolPlot(interp2(T,2),interp2(D,2),interp2(Signlog(VarP,-20,-16),2));
uistack(h,'bottom');
uistack(h,'up');
caxis([-1 1]);
ylabel('Depth (m)','Interpreter','latex','FontSize',25);
text(285.4,-130,'$HCL (s^{-5})$','Interpreter','latex','FontSize',30, ...
     'Color','k');
set(gca,'Position',[0.13 0.5583 0.7297 0.4080]);
xlabel('');
set(gca,'XTick',TimeVector(tposvec)); 
set(gca,'XTickLabel',[]);%{LabelList{[1 2 4 6:length(LabelList)]}});
cb = colorbar;
set(cb,'YTick',[-1:0.5:1]);
set(cb,'YTickLabel',{'<-1e-16','-1e-18','||<1e-20','1e-18','>1e-16'});

subplot(2,1,2);
[T,D,VarP] = FltProfile(flt_pos,hname,dname,DVADV,domain,HCLf_h,HCLf_v, ...
                        tposvec);
hold on;
h = pcolPlot(interp2(T,2),interp2(D,2),interp2(Signlog(VarP,-20,-16),2));
uistack(h,'bottom');
uistack(h,'up');
caxis([-1 1]);
xlabel('Day','Interpreter','latex','FontSize',25);
ylabel('Depth (m)','Interpreter','latex','FontSize',25);
text(285.2,-130,'$DVADV (s^{-5})$','Interpreter','latex','FontSize',30, ...
     'Color','k');
set(gca,'Position',[0.13 0.1265 0.7297 0.3934]);
cb = colorbar;
set(cb,'YTick',[-1:0.5:1]);
set(cb,'YTickLabel',{'<-1e-16','-1e-18','||<1e-20','1e-18','>1e-16'});
set(gca,'XTick',[282:0.5:286]); 

load('my_colormaps.mat');
colormap(red_blue);

%%%%%%%%%PV terms:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PVv_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname,{'PVv'},domain,...
                           h_filter,v_filter,t_filter);
PVh_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname,{'PVh'},domain,...
                           h_filter,v_filter,t_filter);
PV_timeseries = PVv_timeseries+PVh_timeseries;
fN2_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname,{'N2' 'f' '(1).*f_cor'},domain,...
                           h_filter,v_filter,t_filter);
zetaN2_timeseries = PVv_timeseries-fN2_timeseries;

%5-element line colors vector:
cc = fliplr(jet(7));
cc = jet(7)
cc = cc(7:-1:4,:);
cc(5,:) = [0 0 1];
cc(4,:) = 0;
cc(3,:) = [0 0.5 0];
fig3 = figure;
set(fig3, 'Position', [314 38 1175 936]);
subplot(3,1,1);
%PVv:
plot(TimeVector,mean(PVv_timeseries,1),'-','Color',cc(5,:),'LineWidth',2,'DisplayName','PVv');
hold on;
ebar(TimeVector(8:16:end),mean(PVv_timeseries(:,8:16:end),1), ...
     std(PVv_timeseries(:,8:16:end),1),1/8,'-y','Color',cc(5,:),'LineWidth',2);

%PVh:
plot(TimeVector,mean(PVh_timeseries,1),'-','Color',cc(1,:),'LineWidth',2,'DisplayName','PVh');
ebar(TimeVector(11:16:end),mean(PVh_timeseries(:,11:16:end),1),std(PVh_timeseries(:,11:16:end),1), ...
         1/8,'-y','Color',cc(1,:),'LineWidth',2);
%PV:
plot(TimeVector,mean(PV_timeseries,1),'-','Color',cc(4,:),'LineWidth',2,'DisplayName','PV');
ebar(TimeVector(9:16:end),mean(PV_timeseries(:,9:16:end),1),std(PV_timeseries(:,9:16:end),1), ...
         1/8,'-y','Color',cc(4,:),'LineWidth',2);

%fN2:
plot(TimeVector,mean(fN2_timeseries,1),'-','Color',cc(3,:),'LineWidth',2,'DisplayName','fN^2');
ebar(TimeVector(7:16:end),mean(fN2_timeseries(:,7:16:end),1),std(fN2_timeseries(:,7:16:end),1), ...
         1/8,'-y','Color',cc(3,:),'LineWidth',2);
%zetaN2:
plot(TimeVector,mean(zetaN2_timeseries,1),'-','Color',cc(2,:),'LineWidth',2,'DisplayName','zetaN^2');
ebar(TimeVector(10:16:end),mean(zetaN2_timeseries(:,10:16:end),1),std(zetaN2_timeseries(:,10:16:end),1), ...
         1/8,'-y','Color',cc(2,:),'LineWidth',2);

text(282.2,1.6e-9,'$q_v$','Interpreter','latex','FontSize',30,'Color',cc(5,:));
text(283.25,-1.6e-9,'$q_h$','Interpreter','latex','FontSize',30, ...
     'Color',cc(1,:));
text(285.7,2e-9,'$f N^2$','Interpreter','latex','FontSize',30,'Color',cc(3,:));
text(285.7,-1.9e-9,'$\zeta N^2$','Interpreter','latex','FontSize',30,'Color',cc(2,:));
text(283.43,1.1e-9,'$q$','Interpreter','latex','FontSize',30,'Color',cc(4,:));
plot(custom_tlims,[0 0],'-k','LineWidth',1);
xlabel('Day','Interpreter','latex','FontSize',25);
ylabel('PV ($10^{-9} s^{-3}$)','Interpreter','latex','FontSize',25);
set(gca,'YTick',[-0.2 -0.1 0 0.1 0.2 0.3]*1e-8);
set(gca,'YTickLabel',{'-2','-1','0','1','2','3'});
set(gca,'YLim',[-0.2 0.3]*1e-8);
set(gca,'XLim',[282 custom_tlims(2)]);
%set(gca,'XTick',[]);

%%%%%%%%%Zeta Eqn %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zeta_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname, ...
  {'u' 'v' 'cdiff(x_rho,(2),''x'')-cdiff(y_rho,(1),''y'')'}...
                                 ,domain, h_filter,v_filter,t_filter);
zetaFRIC_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname, ...
  {'u_visc' 'v_visc' 'cdiff(x_rho,(2),''x'')-cdiff(y_rho,(1),''y'')'}...
                                 ,domain, h_filter,v_filter,t_filter);
zetaSTR_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname, ...
  {'u' 'v' 'w' '-(f_cor+cdiff(x_rho,(2),''x'')-cdiff(y_rho,(1),''y'')).*cdiff(zvec,(3),''z'')'}...
                                 ,domain, h_filter,v_filter,t_filter);
zetaTILTch_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname, ...
  {'u' 'v' 'w' 'cdiff(zvec,(1),''z'').*cdiff(y_rho,(3),''y'')-cdiff(zvec,(2),''z'').*cdiff(x_rho,(3),''x'')'}...
                                 ,domain, h_filter,v_filter,t_filter);
%Check: this does not correspond very well to the residual
%TILT. Check the filtering? 

%Time integrate:
zetaBETA_timeseries = -2*7.292e-5*(sin(flt_pos(:,:,6)/180*pi)- ...
                                       repmat(sin(flt_pos(:,1,6)/180*pi),[1 tL]));

zetaFRIC_timeseries = Tint(TimeVector,zetaFRIC_timeseries);
zetaTILTch_timeseries = Tint(TimeVector,zetaTILTch_timeseries);
zetaSTR_timeseries = Tint(TimeVector,zetaSTR_timeseries);
zetaTILT_timeseries = zeta_timeseries-repmat(zeta_timeseries(:,1),[1 tL])-zetaFRIC_timeseries- ...
    zetaSTR_timeseries-zetaBETA_timeseries;
zeta_timeseries = zeta_timeseries-repmat(zeta_timeseries(:,1),[1 tL]);

%Set zero:
% $$$ zetaBETA_timeseries = zetaBETA_timeseries + repmat(zeta_timeseries(:,1),[1 tL]);
% $$$ zetaFRIC_timeseries = zetaFRIC_timeseries + repmat(zeta_timeseries(:,1),[1 tL]);
% $$$ zetaSTR_timeseries = zetaSTR_timeseries + repmat(zeta_timeseries(:,1),[1 tL]);
% $$$ zetaTILT_timeseries = zetaTILT_timeseries + repmat(zeta_timeseries(:,1),[1 tL]);

%Check
% $$$ zetaACCEL_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname, ...
% $$$   {'u_accel' 'v_accel' 'cdiff(x_rho,(2),''x'')-cdiff(y_rho,(1),''y'')'}...
% $$$                                  ,domain, h_filter,v_filter,t_filter);
% $$$ zetaADVCOR_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname, ...
% $$$   {'u_cor' 'v_cor' 'u_adv' 'v_adv' 'cdiff(x_rho,(2)+(4),''x'')-cdiff(y_rho,(1)+(3),''y'')'}...
% $$$                                  ,domain, h_filter,v_filter,t_filter);
% $$$ zeta_timeseries = Tint(TimeVector,zetaACCEL_timeseries);
% $$$ zetaTILT_timeseries = Tint(TimeVector,zetaADVCOR_timeseries-zetaSTR_timeseries);
subplot(3,1,2);
%5-element line colors vector:
cc = fliplr(jet(7));
cc = jet(7)
cc = cc(7:-1:4,:);
cc(5,:) = [0 0 1];
cc(4,:) = 0;
cc(3,:) = [0 0.5 0];
% $$$ fig3 = figure;
% $$$ set(fig3, 'Position', [1 1 1300 650]);
%Zeta:
plot(TimeVector,mean(zeta_timeseries,1),'-','Color',cc(5,:),'LineWidth',2,'DisplayName','PVv');
hold on;
ebar(TimeVector(8:16:end),mean(zeta_timeseries(:,8:16:end),1), ...
     std(zeta_timeseries(:,16:end),1),1/8,'-y','Color',cc(5,:),'LineWidth',2);

%STR:
plot(TimeVector,mean(zetaSTR_timeseries,1),'-','Color',cc(1,:),'LineWidth',2,'DisplayName','PVh');
ebar(TimeVector(11:16:end),mean(zetaSTR_timeseries(:,11:16:end),1),std(zetaSTR_timeseries(:,11:16:end),1), ...
         1/8,'-y','Color',cc(1,:),'LineWidth',2);
%TILT:
plot(TimeVector,mean(zetaTILT_timeseries,1),'-','Color',cc(4,:),'LineWidth',2,'DisplayName','PV');
ebar(TimeVector(9:16:end),mean(zetaTILT_timeseries(:,9:16:end),1),std(zetaTILT_timeseries(:,9:16:end),1), ...
         1/8,'-y','Color',cc(4,:),'LineWidth',2);
%TILTch:
plot(TimeVector,mean(zetaTILTch_timeseries,1),':','Color',cc(4,:),'LineWidth',2,'DisplayName','PV');
ebar(TimeVector(9:16:end),mean(zetaTILTch_timeseries(:,9:16:end),1),std(zetaTILTch_timeseries(:,9:16:end),1), ...
         1/8,'-y','Color',cc(4,:),'LineWidth',2);

%FRIC:
plot(TimeVector,mean(zetaFRIC_timeseries,1),'-','Color',cc(3,:),'LineWidth',2,'DisplayName','fN^2');
ebar(TimeVector(7:16:end),mean(zetaFRIC_timeseries(:,7:16:end),1),std(zetaFRIC_timeseries(:,7:16:end),1), ...
         1/8,'-y','Color',cc(3,:),'LineWidth',2);
%BETA:
plot(TimeVector,mean(zetaBETA_timeseries,1),'-','Color',cc(2,:),'LineWidth',2,'DisplayName','fN^2');
ebar(TimeVector(10:16:end),mean(zetaBETA_timeseries(:,10:16:end),1),std(zetaBETA_timeseries(:,10:16:end),1), ...
         1/8,'-y','Color',cc(2,:),'LineWidth',2);

text(284,1e-5,'$\zeta$','Interpreter','latex','FontSize',30,'Color',cc(5,:));
text(285.5,2.1e-5,'$STR$','Interpreter','latex','FontSize',20,'Color',cc(1,:));
text(285.5,1.3e-5,'$TILT$','Interpreter','latex','FontSize',20,'Color',cc(4,:));
text(283,1.4e-5,'$FRIC$','Interpreter','latex','FontSize',20,'Color',cc(3,:));
text(283,1.4e-5,'$BETA$','Interpreter','latex','FontSize',20,'Color',cc(2,:));
plot(custom_tlims,[0 0],'-k','LineWidth',1);
xlabel('Day','Interpreter','latex','FontSize',25);
ylabel('$\zeta$ ($10^{-5} s^{-1}$)','Interpreter','latex','FontSize',25);
set(gca,'YTick',[-3 -2 -1 0 1 2 3]*1e-5);
set(gca,'YTickLabel',{'-3','-2','-1','0','1','2','3'});
set(gca,'YLim',[-3.7 2.5]*1e-5);
set(gca,'XLim',[282 custom_tlims(2)]);
%set(gca,'XTick',[]);


% $$$ %%%%%%%%%N2 Eqn %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ N2_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname,{'N2'},domain, ...
% $$$                           h_filter,v_filter,t_filter);
% $$$ N2DIAB_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname, {'b_diff' ...
% $$$                     'v' 'cdiff(zvec,(1),''z'')'} ,domain, h_filter, ...
% $$$                               v_filter,t_filter);
% $$$ N2STR_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname, {'N2' ...
% $$$                     'w' '-(1).*cdiff(zvec,(2),''z'')'} ,domain, ...
% $$$                              h_filter,v_filter,t_filter);
% $$$ N2TILTch_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname, {'DBDX' ...
% $$$                     'DBDY' 'u' 'v' '-cdiff(zvec,(3),''z'').*(1)-cdiff(zvec,(4),''z'').*(2)'} ...
% $$$                                 ,domain, h_filter,v_filter,t_filter);
% $$$ %Check: This does not correspond very well to the residual
% $$$ %TILT. -numerical issues. 
% $$$ 
% $$$ %Time integrate:
% $$$ N2DIAB_timeseries = Tint(TimeVector,N2DIAB_timeseries);
% $$$ N2STR_timeseries = Tint(TimeVector,N2STR_timeseries);
% $$$ N2TILTch_timeseries = Tint(TimeVector,N2TILTch_timeseries);
% $$$ N2TILT_timeseries = N2_timeseries-repmat(N2_timeseries(:,1),[1 tL])-N2DIAB_timeseries-N2STR_timeseries;
% $$$ N2_timeseries = N2_timeseries-repmat(N2_timeseries(:,1),[1 tL]);
% $$$ 
% $$$ %Set zero:
% $$$ % $$$ N2DIAB_timeseries = N2DIAB_timeseries + mean(N2_timeseries(:,1));
% $$$ % $$$ N2STR_timeseries = N2STR_timeseries + mean(N2_timeseries(:,1));
% $$$ % $$$ N2TILT_timeseries = N2TILT_timeseries + mean(N2_timeseries(:,1));
% $$$ 
% $$$ %5-element line colors vector:
% $$$ cc = fliplr(jet(7));
% $$$ cc = jet(7)
% $$$ cc = cc(7:-1:4,:);
% $$$ cc(5,:) = [0 0 1];
% $$$ cc(4,:) = 0;
% $$$ cc(3,:) = [0 0.5 0];
% $$$ fig3 = figure;
% $$$ set(fig3, 'Position', [1 1 1300 650]);
% $$$ %N2:
% $$$ plot(TimeVector,mean(N2_timeseries,1),'-','Color',cc(5,:),'LineWidth',2,'DisplayName','PVv');
% $$$ hold on;
% $$$ ebar(TimeVector(8:16:end),mean(N2_timeseries(:,8:16:end),1), ...
% $$$      std(N2_timeseries(:,16:end),1),1/8,'-y','Color',cc(5,:),'LineWidth',2);
% $$$ 
% $$$ %STR:
% $$$ plot(TimeVector,mean(N2STR_timeseries,1),'-','Color',cc(1,:),'LineWidth',2,'DisplayName','PVh');
% $$$ ebar(TimeVector(11:16:end),mean(N2STR_timeseries(:,11:16:end),1),std(N2STR_timeseries(:,11:16:end),1), ...
% $$$          1/8,'-y','Color',cc(1,:),'LineWidth',2);
% $$$ %TILT:
% $$$ plot(TimeVector,mean(N2TILT_timeseries,1),'-','Color',cc(4,:),'LineWidth',2,'DisplayName','PV');
% $$$ ebar(TimeVector(9:16:end),mean(N2TILT_timeseries(:,9:16:end),1),std(N2TILT_timeseries(:,9:16:end),1), ...
% $$$          1/8,'-y','Color',cc(4,:),'LineWidth',2);
% $$$ %TILTch:
% $$$ plot(TimeVector,mean(N2TILTch_timeseries,1),':','Color',cc(4,:),'LineWidth',2,'DisplayName','PV');
% $$$ ebar(TimeVector(9:16:end),mean(N2TILTch_timeseries(:,9:16:end),1),std(N2TILT_timeseries(:,9:16:end),1), ...
% $$$          1/8,'-y','Color',cc(4,:),'LineWidth',2);
% $$$ 
% $$$ %DIAB:
% $$$ plot(TimeVector,mean(N2DIAB_timeseries,1),'-','Color',cc(3,:),'LineWidth',2,'DisplayName','fN^2');
% $$$ ebar(TimeVector(7:16:end),mean(N2DIAB_timeseries(:,7:16:end),1),std(N2DIAB_timeseries(:,7:16:end),1), ...
% $$$          1/8,'-y','Color',cc(3,:),'LineWidth',2);
% $$$ 
% $$$ text(284,0,'$N^2$','Interpreter','latex','FontSize',20,'Color',cc(5,:));
% $$$ text(285.5,0,'$N^2STR$','Interpreter','latex','FontSize',20,'Color',cc(1,:));
% $$$ text(285.5,0,'$N^2TILT$','Interpreter','latex','FontSize',20,'Color',cc(4,:));
% $$$ text(283,0,'$N^2DIAB$','Interpreter','latex','FontSize',20,'Color',cc(3,:));
% $$$ plot(custom_tlims,[0 0],'-k','LineWidth',1);
% $$$ xlabel('Day','Interpreter','latex','FontSize',25);
% $$$ ylabel('$N^2$ ($10^{-4} s^{-2}$)','Interpreter','latex','FontSize',25);
% $$$ set(gca,'YTick',[-1:0.25:1]*1e-4);
% $$$ set(gca,'YTickLabel',{'-1','-0.75','-0.5','-0.25','0','0.25','0.5','0.75','1'});
% $$$ set(gca,'YLim',[-1.2 0.75]*1e-4);
% $$$ set(gca,'XLim',[282 custom_tlims(2)]);
% $$$ %set(gca,'XTick',[]);


%%%%%%%%%Gradhb Eqn %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gradhb_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname,{'DBDX' ...
                    'DBDY' '(1).^2+(2).^2'},domain, h_filter, ...
                              v_filter,t_filter);
gradhbDIAB_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname, ...
                                  {'b_diff' 'DBDX' 'DBDY' ...
                    ['2*(cdiff(x_rho,(1),''x'').*(2)+cdiff(y_rho,' ...
                    '(1),''y'').*(3))']},domain, h_filter, v_filter,t_filter);
gradhbvDIAB_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname, ...
                                  {'b_vdiff' 'DBDX' 'DBDY' ...
                    ['2*(cdiff(x_rho,(1),''x'').*(2)+cdiff(y_rho,' ...
                    '(1),''y'').*(3))']},domain, h_filter, v_filter,t_filter);
gradhbDVADVch_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname, ...
      {'DBDX' 'DBDY' 'w' 'N2' ['-2*(4).*(cdiff(x_rho,(3),''x'').*(1)+cdiff(y_rho,' ...
                    '(3),''y'').*(2))']} ,domain, h_filter, ...
                                   v_filter,t_filter);
gradhbSHR_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname, ...
      {'DBDX' 'DBDY' 'u' 'v' ['-2*(1).*(2).*(cdiff(x_rho,(4),''x'')+cdiff(y_rho,' ...
                    '(3),''y''))']} ,domain, h_filter, ...
                                   v_filter,t_filter);
gradhbSQZ_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname, ...
      {'DBDX' 'DBDY' 'u' 'v' ['-2*(cdiff(x_rho,(3),''x'').*(1).^2+cdiff(y_rho,' ...
                    '(4),''y'').*(2).^2)']} ,domain, h_filter, ...
                                   v_filter,t_filter);
gradhbHCL_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname, ...
      {'DBDX' 'DBDY' 'u' 'v' ['-2*(1).*(2).*(cdiff(x_rho,(4),''x'')+cdiff(y_rho,' ...
                    '(3),''y''))-2*(cdiff(x_rho,(3),''x'').*(1).^2+cdiff(y_rho,' ...
                    '(4),''y'').*(2).^2)']} ,domain, h_filter, ...
                                   v_filter,t_filter);
%Check: This corresponds reasonably well to the residual SQZ. 

%Time integrate:
gradhbhDIAB_timeseries = gradhbDIAB_timeseries-gradhbvDIAB_timeseries;
gradhbDIAB_timeseries = Tint(TimeVector,gradhbDIAB_timeseries);
gradhbvDIAB_timeseries = Tint(TimeVector,gradhbvDIAB_timeseries);
gradhbhDIAB_timeseries = Tint(TimeVector,gradhbhDIAB_timeseries);
gradhbDVADVch_timeseries = Tint(TimeVector,gradhbDVADVch_timeseries);
gradhbSHR_timeseries = Tint(TimeVector,gradhbSHR_timeseries);
gradhbSQZ_timeseries = Tint(TimeVector,gradhbSQZ_timeseries);
gradhbHCL_timeseries = Tint(TimeVector,gradhbHCL_timeseries);
gradhbDVADV_timeseries = gradhb_timeseries-repmat(gradhb_timeseries(:,1),[1 ...
                    tL])-gradhbDIAB_timeseries-gradhbHCL_timeseries;
%gradhb_timeseries = gradhb_timeseries-repmat(gradhb_timeseries(:,1),[1 tL]);

%No time integrals:
% $$$ gradhb_timeseries(:,2:end) = (gradhb_timeseries(:,2:end)- ...
% $$$                                   gradhb_timeseries(:,1:(end-1)))/ ...
% $$$     (TimeVector(2)-TimeVector(1))/86400;
% $$$ gradhb_timeseries(:,1) = 0;
% $$$ gradhbCON_timeseries = gradhb_timeseries-repmat(gradhb_timeseries(:,1),[1 ...
% $$$                     tL])-gradhbDIAB_timeseries-gradhbDVADV_timeseries;
% $$$ gradhbSHR_timeseries = gradhbCON_timeseries-gradhbSQZ_timeseries;


%5-element line colors vector:
cc = fliplr(jet(7));
cc = jet(7)
cc = cc(7:-1:4,:);
cc(5,:) = [0 0 1];
cc(4,:) = 0;
cc(3,:) = [0 0.5 0];
% $$$ fig3 = figure;
% $$$ set(fig3, 'Position', [1 1 1300 650]);
subplot(3,1,3);
%gradhb:
plot(TimeVector,mean(gradhb_timeseries,1),'-','Color',cc(5,:),'LineWidth',2);
hold on;
ebar(TimeVector(8:8:end),mean(gradhb_timeseries(:,8:8:end),1), ...
     std(gradhb_timeseries(:,8:8:end),1),1/8,'-y','Color',cc(5,:),'LineWidth',2);

%DVADV:
plot(TimeVector,mean(gradhbDVADV_timeseries,1),'-','Color',cc(1,:),'LineWidth',2);
ebar(TimeVector(11:8:end),mean(gradhbDVADV_timeseries(:,11:8:end),1),std(gradhbDVADV_timeseries(:,11:8:end),1), ...
         1/8,'-y','Color',cc(1,:),'LineWidth',2);
%DVADVch:
plot(TimeVector,mean(gradhbDVADVch_timeseries,1),':','Color',cc(1,:),'LineWidth',2);
ebar(TimeVector(9:8:end),mean(gradhbDVADVch_timeseries(:,9:8:end),1),std(gradhbDVADVch_timeseries(:,9:8:end),1), ...
         1/8,'-y','Color',cc(1,:),'LineWidth',2);
%SQZ:
plot(TimeVector,mean(gradhbSQZ_timeseries,1),'-','Color',cc(4,:),'LineWidth',2);
ebar(TimeVector(9:8:end),mean(gradhbSQZ_timeseries(:,9:8:end),1),std(gradhbSQZ_timeseries(:,9:8:end),1), ...
         1/8,'-y','Color',cc(4,:),'LineWidth',2);

%DIAB:
plot(TimeVector,mean(gradhbDIAB_timeseries,1),'-','Color',cc(3,:),'LineWidth',2);
ebar(TimeVector(7:8:end),mean(gradhbDIAB_timeseries(:,7:8:end),1),std(gradhbDIAB_timeseries(:,7:8:end),1), ...
         1/8,'-y','Color',cc(3,:),'LineWidth',2);

%DIABv:
plot(TimeVector,mean(gradhbvDIAB_timeseries,1),'--','Color',cc(3,:),'LineWidth',2);
% $$$ ebar(TimeVector(7:8:end),mean(gradhbvDIAB_timeseries(:,7:8:end),1),std(gradhbvDIAB_timeseries(:,7:8:end),1), ...
% $$$          1/8,'-y','Color',cc(3,:),'LineWidth',2);

%DIABh:
plot(TimeVector,mean(gradhbhDIAB_timeseries,1),':','Color',cc(3,:),'LineWidth',2);
% $$$ ebar(TimeVector(7:8:end),mean(gradhbhDIAB_timeseries(:,7:8:end),1),std(gradhbhDIAB_timeseries(:,7:8:end),1), ...
% $$$          1/8,'-y','Color',cc(3,:),'LineWidth',2);

%SHR:
plot(TimeVector,mean(gradhbSHR_timeseries,1),'-','Color',cc(2,:),'LineWidth',2);
ebar(TimeVector(10:8:end),mean(gradhbSHR_timeseries(:,10:8:end),1),std(gradhbSHR_timeseries(:,10:8:end),1), ...
         1/8,'-y','Color',cc(2,:),'LineWidth',2);

%HCL:
plot(TimeVector,mean(gradhbHCL_timeseries,1),'-','Color',[0.789 0 0.789],'LineWidth',2);
ebar(TimeVector(10:8:end),mean(gradhbHCL_timeseries(:,10:8:end),1),std(gradhbHCL_timeseries(:,10:8:end),1), ...
         1/8,'-y','Color',[0.789 0 0.789],'LineWidth',2);

%All Dyn:
gradhbDYN_timeseries = gradhbSHR_timeseries+gradhbSQZ_timeseries+gradhbDVADV_timeseries;
plot(TimeVector,mean(gradhbDYN_timeseries,1),'-','Color',cc(2,:),'LineWidth',2);
ebar(TimeVector(10:8:end),mean(gradhbDYN_timeseries(:,10:8:end),1),std(gradhbSHR_timeseries(:,10:8:end),1), ...
         1/8,'-y','Color',cc(2,:),'LineWidth',2);

text(284,0,'$|\nabla_h b|^2$','Interpreter','latex','FontSize',20,'Color',cc(5,:));
text(285.5,0,'$DVADV$','Interpreter','latex','FontSize',20,'Color',cc(1,:));
text(285.5,0,'$SQZ$','Interpreter','latex','FontSize',20,'Color',cc(4,:));
text(283,0,'$DIAB$','Interpreter','latex','FontSize',20,'Color',cc(3,:));
text(283,0,'$SHR$','Interpreter','latex','FontSize',20,'Color',cc(2,:));
text(283,0,'$DYN$','Interpreter','latex','FontSize',20,'Color',cc(2,:));
plot(custom_tlims,[0 0],'-k','LineWidth',1);
xlabel('Day','Interpreter','latex','FontSize',25);
ylabel('$|\nabla_h b|^2$ ($10^{-13} s^{-4}$)','Interpreter','latex','FontSize',25);
set(gca,'YTick',[-0.5:0.25:0.75]*1e-13);
set(gca,'YTickLabel',{'-0.5','-0.25','0','0.25','0.5','0.75'});
set(gca,'YLim',[-0.5 0.75]*1e-13);
set(gca,'XLim',[282 custom_tlims(2)]);
%set(gca,'XTick',[]);


% $$$ %%%%%%%%%omgh Eqn %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ omgh_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname,{'u' ...
% $$$                     'v' 'cdiff(zvec,(2),''z'').^2+cdiff(zvec,(1),''z'').^2'},domain, h_filter, ...
% $$$                               v_filter,t_filter);
% $$$ omghDDT_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname,{'u' ...
% $$$                     'v' 'u_accel' 'v_accel' 'u_adv' 'v_adv' ...
% $$$                     ['-2*cdiff(zvec,(2),''z'').*(-cdiff(zvec,(4)-(6),' ...
% $$$                     '''z''))+2*cdiff(zvec,(1),''z'').*cdiff(zvec,' ...
% $$$                     '(3)-(5),''z'')']},domain, h_filter, v_filter, ...
% $$$                                t_filter);
% $$$ %Check: time integrating this does not give the same results as
% $$$ %omgh_timeseries! But since it gives a negative results I don't
% $$$ %trust it as much as the omgh_timeseries.
% $$$ 
% $$$ omghBTOR_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname,{'u' ...
% $$$                     'v' 'u_prsgrd' 'v_prsgrd' ...
% $$$                     ['-2*cdiff(zvec,(2),''z'').*(-cdiff(zvec,(4),' ...
% $$$                     '''z''))+2*cdiff(zvec,(1),''z'').*cdiff(zvec,' ...
% $$$                     '(3),''z'')']},domain, h_filter, v_filter,t_filter);
% $$$ omghFRIC_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname,{'u' ...
% $$$                     'v' 'u_visc' 'v_visc' ...
% $$$                     ['-2*cdiff(zvec,(2),''z'').*(-cdiff(zvec,(4),' ...
% $$$                     '''z''))+2*cdiff(zvec,(1),''z'').*cdiff(zvec,' ...
% $$$                     '(3),''z'')']},domain, h_filter, v_filter,t_filter);
% $$$ % $$$ omghCOR_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname,{'u' ...
% $$$ % $$$                     'v' 'u_cor' 'v_cor' ...
% $$$ % $$$                     ['-2*cdiff(zvec,(2),''z'').*(-cdiff(zvec,(4),' ...
% $$$ % $$$                     '''z''))+2*cdiff(zvec,(1),''z'').*cdiff(zvec,' ...
% $$$ % $$$                     '(3),''z'')']},domain, h_filter, v_filter,t_filter);
% $$$ %Check: this is zero.
% $$$ 
% $$$ omghFRICv_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname, ...
% $$$                                   {'u_vvisc' 'v_vvisc' 'u' 'v' ...
% $$$                     ['2*cdiff(zvec,(4),''z'').*cdiff(zvec,(2),''z'')+'...
% $$$                      '2*cdiff(zvec,(3),''z'').*cdiff(zvec,(1),''z'')']},domain, h_filter, v_filter,t_filter);
% $$$ 
% $$$ % $$$ omghBTOR_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname, ...
% $$$ % $$$                                   {'DBDX' 'DBDY' 'u' 'v' ...
% $$$ % $$$                     ['-2*cdiff(zvec,(4),''z'').*(2)-'...
% $$$ % $$$                      '2*cdiff(zvec,(3),''z'').*(1)']},domain, h_filter, v_filter,t_filter);
% $$$ %Check: this is same as prsgrd method.
% $$$ 
% $$$ omghSQH_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname, {'u' ...
% $$$                     'v' ['2*cdiff(zvec,(2),''z'').^2.*cdiff(x_rho,(1),''x'')+'...
% $$$                     '2*cdiff(zvec,(1),''z'').^2.*cdiff(y_rho,(2),''y'')']} ,domain, h_filter, ...
% $$$                                v_filter,t_filter);
% $$$ omghTILT_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname, {'u' ...
% $$$                     'v' ['-2*cdiff(zvec,(2),''z'').*cdiff(zvec,(1),''z'').*('...
% $$$                     'cdiff(y_rho,(1),''y'')+cdiff(x_rho,(2),''x''))']} ,domain, h_filter, ...
% $$$                                v_filter,t_filter);
% $$$ 
% $$$ omghTILTEk_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname, ...
% $$$                                   {'u' 'v' 'u_vvisc' 'v_vvisc' ...
% $$$                     ['-2*cdiff(zvec,-(3)./f_cor,''z'').*' ...
% $$$                     'cdiff(zvec,(4)./f_cor,''z'').*(' ...
% $$$                     'cdiff(y_rho,(1),''y'')+cdiff(x_rho,(2),''x''))']} ...
% $$$     ,domain, h_filter, v_filter,t_filter);
% $$$ 
% $$$ %Time integrate:
% $$$ omghFRICh_timeseries = omghFRIC_timeseries-omghFRICv_timeseries;
% $$$ omghFRIC_timeseries = Tint(TimeVector,omghFRIC_timeseries);
% $$$ omghFRICv_timeseries = Tint(TimeVector,omghFRICv_timeseries);
% $$$ omghFRICh_timeseries = Tint(TimeVector,omghFRICh_timeseries);
% $$$ omghBTOR_timeseries = Tint(TimeVector,omghBTOR_timeseries);
% $$$ %omghSQH_timeseries = Tint(TimeVector,omghSQH_timeseries);
% $$$ %omghDDT_timeseries = Tint(TimeVector,omghDDT_timeseries);
% $$$ omghTILTEk_timeseries = Tint(TimeVector,omghTILTEk_timeseries);
% $$$ omghTILT_timeseries = Tint(TimeVector,omghTILT_timeseries);
% $$$ omghSQH_timeseries = omgh_timeseries-repmat(omgh_timeseries(:,1),[1 ...
% $$$                     tL])-omghFRIC_timeseries-omghBTOR_timeseries- ...
% $$$     omghTILT_timeseries;
% $$$ %omgh_timeseries = omgh_timeseries-repmat(omgh_timeseries(:,1),[1 tL]);
% $$$ 
% $$$ %5-element line colors vector:
% $$$ cc = fliplr(jet(7));
% $$$ cc = jet(7)
% $$$ cc = cc(7:-1:4,:);
% $$$ cc(5,:) = [0 0 1];
% $$$ cc(4,:) = 0;
% $$$ cc(3,:) = [0 0.5 0];
% $$$ fig3 = figure;
% $$$ set(fig3, 'Position', [1 1 1300 650]);
% $$$ %omgh:
% $$$ plot(TimeVector,mean(omgh_timeseries,1),'-','Color',cc(5,:),'LineWidth',2);
% $$$ hold on;
% $$$ ebar(TimeVector(8:16:end),mean(omgh_timeseries(:,8:16:end),1), ...
% $$$      std(omgh_timeseries(:,16:end),1),1/8,'-y','Color',cc(5,:),'LineWidth',2);
% $$$ 
% $$$ %BTOR:
% $$$ plot(TimeVector,mean(omghBTOR_timeseries,1),'-','Color',cc(1,:),'LineWidth',2);
% $$$ ebar(TimeVector(11:16:end),mean(omghBTOR_timeseries(:,11:16:end),1),std(omghBTOR_timeseries(:,11:16:end),1), ...
% $$$          1/8,'-y','Color',cc(1,:),'LineWidth',2);
% $$$ %SQH:
% $$$ plot(TimeVector,mean(omghSQH_timeseries,1),'-','Color',cc(4,:),'LineWidth',2);
% $$$ ebar(TimeVector(9:16:end),mean(omghSQH_timeseries(:,9:16:end),1),std(omghSQH_timeseries(:,9:16:end),1), ...
% $$$          1/8,'-y','Color',cc(4,:),'LineWidth',2);
% $$$ 
% $$$ %FRIC:
% $$$ plot(TimeVector,mean(omghFRIC_timeseries,1),'-','Color',cc(3,:),'LineWidth',2);
% $$$ ebar(TimeVector(7:16:end),mean(omghFRIC_timeseries(:,7:16:end),1),std(omghFRIC_timeseries(:,7:16:end),1), ...
% $$$          1/8,'-y','Color',cc(3,:),'LineWidth',2);
% $$$ 
% $$$ %FRICv:
% $$$ plot(TimeVector,mean(omghFRICv_timeseries,1),'--','Color',cc(3,:),'LineWidth',2);
% $$$ % $$$ ebar(TimeVector(7:16:end),mean(omghFRICv_timeseries(:,7:16:end),1),std(omghFRICv_timeseries(:,7:16:end),1), ...
% $$$ % $$$          1/8,'-y','Color',cc(3,:),'LineWidth',2);
% $$$ 
% $$$ %FRICh:
% $$$ % $$$ plot(TimeVector,mean(omghFRICh_timeseries,1),':','Color',cc(3,:),'LineWidth',2);
% $$$ % $$$ ebar(TimeVector(7:16:end),mean(omghFRICh_timeseries(:,7:16:end),1),std(omghFRICh_timeseries(:,7:16:end),1), ...
% $$$ % $$$          1/8,'-y','Color',cc(3,:),'LineWidth',2);
% $$$ 
% $$$ %TILT:
% $$$ plot(TimeVector,mean(omghTILT_timeseries,1),'-','Color',cc(2,:),'LineWidth',2);
% $$$ ebar(TimeVector(10:16:end),mean(omghTILT_timeseries(:,10:16:end),1),std(omghTILT_timeseries(:,10:16:end),1), ...
% $$$          1/8,'-y','Color',cc(2,:),'LineWidth',2);
% $$$ plot(TimeVector,mean(omghTILTEk_timeseries,1),':','Color',cc(2,:),'LineWidth',2);
% $$$ ebar(TimeVector(10:16:end),mean(omghTILTEk_timeseries(:,10:16:end),1),std(omghTILTEk_timeseries(:,10:16:end),1), ...
% $$$          1/8,'-y','Color',cc(2,:),'LineWidth',2);
% $$$ 
% $$$ % $$$ %SQHch:
% $$$ % $$$ plot(TimeVector,mean(omghSQHch2_timeseries,1),':','Color',cc(4,:),'LineWidth',2);
% $$$ % $$$ ebar(TimeVector(10:16:end),mean(omghSQHch2_timeseries(:,10:16:end),1),std(omghSQHch2_timeseries(:,10:16:end),1), ...
% $$$ % $$$          1/8,'-y','Color',cc(4,:),'LineWidth',2);
% $$$ % $$$ 
% $$$ % $$$ %DDT:
% $$$ % $$$ plot(TimeVector,mean(omghDDT_timeseries,1),':','Color',cc(5,:),'LineWidth',2);
% $$$ % $$$ ebar(TimeVector(10:16:end),mean(omghDDT_timeseries(:,10:16:end),1),std(omghDDT_timeseries(:,10:16:end),1), ...
% $$$ % $$$          1/8,'-y','Color',cc(5,:),'LineWidth',2);
% $$$ % $$$ 
% $$$ % $$$ %DDT:
% $$$ % $$$ omghDDTn_timeseries = omghDDT_timeseries+repmat(omgh_timeseries(:,1),[1 tL]);
% $$$ % $$$ plot(TimeVector,mean(omghDDTn_timeseries,1),':','Color',cc(5,:),'LineWidth',2);
% $$$ % $$$ ebar(TimeVector(10:16:end),mean(omghDDTn_timeseries(:,10:16:end),1),std(omghDDTn_timeseries(:,10:16:end),1), ...
% $$$ % $$$          1/8,'-y','Color',cc(5,:),'LineWidth',2);
% $$$ 
% $$$ text(284,0,'$|\omega_h|^2$','Interpreter','latex','FontSize',20,'Color',cc(5,:));
% $$$ text(285.5,0,'$BTOR$','Interpreter','latex','FontSize',20,'Color',cc(1,:));
% $$$ text(285.5,0,'$SQH$','Interpreter','latex','FontSize',20,'Color',cc(4,:));
% $$$ text(283,0,'$FRIC$','Interpreter','latex','FontSize',20,'Color',cc(3,:));
% $$$ text(283,0,'$FRIC_v$','Interpreter','latex','FontSize',20,'Color',cc(3,:));
% $$$ text(283,0,'$TILT$','Interpreter','latex','FontSize',20,'Color',cc(2,:));
% $$$ plot(custom_tlims,[0 0],'-k','LineWidth',1);
% $$$ xlabel('Day','Interpreter','latex','FontSize',25);
% $$$ ylabel('$|\omega_h|^2$ ($10^{-3} s^{-2}$)','Interpreter','latex','FontSize',25);
% $$$ set(gca,'YTick',[-2:0.5:1.5]*1e-3);
% $$$ set(gca,'YTickLabel',{'-2','-1.5','-1','-0.5','0','0.5','1','1.5'});
% $$$ set(gca,'YLim',[-2 1.5]*1e-3);
% $$$ set(gca,'XLim',[282 custom_tlims(2)]);
% $$$ %set(gca,'XTick',[]);

%%%Viscosity along tracks:
% $$$ fig3 = figure;
% $$$ set(fig3, 'Position', [347 4 1176 970]);
% $$$ subplot(4,1,1);
% $$$ FltProfile(flt_pos,hname,dname,{'b_vdiff'},domain,[0 0],0,tposvec);
% $$$ caxis([-0.5 0.5]*1e-8);
% $$$ text(285.6,-130,'$\mathcal{H}_v (ms^{-3})$','Interpreter','latex','FontSize',30, ...
% $$$      'Color','k');
% $$$ subplot(4,1,2);
% $$$ FltProfile(flt_pos,hname,dname,{'b_diff' 'b_vdiff' '(1)-(2)'},domain,[0 0],0,tposvec);
% $$$ caxis([-0.5 0.5]*1e-8);
% $$$ text(285.6,-130,'$\mathcal{H}_h (ms^{-3})$','Interpreter','latex','FontSize',30, ...
% $$$      'Color','k');
% $$$ subplot(4,1,3);
% $$$ FltProfile(flt_pos,hname,dname,{'u_vvisc' 'v_vvisc' 'sqrt((1).^2+(2).^2)'},domain,[0 0],0,tposvec);
% $$$ caxis([0 2]*1e-6);
% $$$ text(285.6,-130,'$|F|_v (ms^{-2})$','Interpreter','latex','FontSize',30, ...
% $$$      'Color','k');
% $$$ subplot(4,1,4);
% $$$ FltProfile(flt_pos,hname,dname,{'u_visc' 'v_visc' 'u_vvisc' 'v_vvisc' 'sqrt(((1)-(3)).^2+((2)-(4)).^2)'},domain,[0 0],0,tposvec);
% $$$ caxis([0 2]*1e-6);
% $$$ text(285.6,-130,'$|F|_h (ms^{-2})$','Interpreter','latex','FontSize',30, ...
% $$$      'Color','k');

%%%%%DBDX + DBDY:
cc = fliplr(jet(7));
cc = jet(7)
cc = cc(7:-1:4,:);
cc(5,:) = [0 0 1];
cc(4,:) = 0;
cc(3,:) = [0 0.5 0];
fig3 = figure;
set(fig3, 'Position', [314 38 1175 600]);

DBDX_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname,{'DBDX'}, ...
                            domain,h_filter,v_filter,t_filter);
DBDY_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname,{'DBDY'}, ...
                            domain,h_filter,v_filter,t_filter);
gradhbmag_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname,{'DBDX' ...
                    'DBDY' 'sqrt((1).^2+(2).^2)'}, domain,h_filter, ...
                              v_filter,t_filter);
gradhbmagsq_timeseries = sqrt(FltInterp(flt_pos(:,:,2:4),hname,dname,{'DBDX' ...
                    'DBDY' '(1).^2+(2).^2'}, domain,h_filter, ...
                              v_filter,t_filter));

%DBDX:
plot(TimeVector,mean(DBDX_timeseries,1),'-','Color',cc(1,:),'LineWidth',2);
hold on;
ebar(TimeVector(7:8:end),mean(DBDX_timeseries(:,7:8:end),1), ...
     std(DBDX_timeseries(:,7:8:end),1),1/8,'-y','Color',cc(1,:),'LineWidth',2);

%DBDY:
plot(TimeVector,mean(DBDY_timeseries,1),'-','Color',cc(2,:),'LineWidth',2);
hold on;
ebar(TimeVector(6:8:end),mean(DBDY_timeseries(:,6:8:end),1), ...
     std(DBDY_timeseries(:,6:8:end),1),1/8,'-y','Color',cc(2,:),'LineWidth',2);

gradhbmagchk_timeseries = sqrt(DBDX_timeseries.^2+DBDY_timeseries.^2);

%Magnitude:
plot(TimeVector,mean(gradhbmag_timeseries,1),'-','Color',cc(3,:),'LineWidth',2);
hold on;
ebar(TimeVector(5:8:end),mean(gradhbmag_timeseries(:,5:8:end),1), ...
     std(gradhbmag_timeseries(:,5:8:end),1),1/8,'-y','Color',cc(3,:),'LineWidth',2);

%Magnitude check:
plot(TimeVector,mean(gradhbmagchk_timeseries,1),'-','Color',cc(4,:),'LineWidth',2);
hold on;
ebar(TimeVector(4:8:end),mean(gradhbmagchk_timeseries(:,4:8:end),1), ...
     std(gradhbmag_timeseries(:,4:8:end),1),1/8,'-y','Color',cc(4,:),'LineWidth',2);


%Magnitude check2:
plot(TimeVector,mean(gradhbmagsq_timeseries,1),'-','Color',cc(5,:),'LineWidth',2);
hold on;
ebar(TimeVector(3:8:end),mean(gradhbmagsq_timeseries(:,3:8:end),1), ...
     std(gradhbmag_timeseries(:,3:8:end),1),1/8,'-y','Color',cc(5,:),'LineWidth',2);

text(284,0,'$|\nabla_h b|$','Interpreter','latex','FontSize',20,'Color',cc(3,:));
text(284,0,'$\frac{\partial b}{\partial x}$','Interpreter','latex','FontSize',20,'Color',cc(1,:));
text(284,0,'$\frac{\partial b}{\partial y}$','Interpreter','latex','FontSize',20,'Color',cc(2,:));
text(284,0,'$\sqrt{\left(\frac{\partial b}{\partial x}\right)^2+\left(\frac{\partial b}{\partial y}\right)^2}$','Interpreter','latex','FontSize',20,'Color',cc(4,:));
text(284,0,'$\sqrt{|\nabla_h b|^2}$','Interpreter','latex','FontSize',20,'Color',cc(5,:));
plot(custom_tlims,[0 0],'-k','LineWidth',1);
xlabel('Day','Interpreter','latex','FontSize',25);
ylabel('$|\nabla_h b|$ ($10^{-7} s^{-2}$)','Interpreter','latex','FontSize',25);
set(gca,'YTick',[-0.5:0.5:2]*1e-7);
set(gca,'YTickLabel',{'-0.5','0','0.5','1','1.5','2'});
set(gca,'YLim',[-0.5 2]*1e-7);
set(gca,'XLim',[282 custom_tlims(2)]);



%%%%%%%%%omghx and omghy%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc = fliplr(jet(7));
cc = jet(7)
cc = cc(7:-1:4,:);
cc(5,:) = [0 0 1];
cc(4,:) = 0;
cc(3,:) = [0 0.5 0];
fig3 = figure;
set(fig3, 'Position', [314 38 1175 600]);

omghx_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname,{'u' ...
                    'v' '-cdiff(zvec,(2),''z'')+0*(1)'},domain, h_filter, ...
                              v_filter,t_filter);
omghy_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname,{'u' ...
                    'v' 'cdiff(zvec,(1),''z'')'},domain, h_filter, ...
                              v_filter,t_filter);
omghmag_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname,{'u' ...
                    'v' 'sqrt(cdiff(zvec,(2),''z'').^2+cdiff(zvec,(1),''z'').^2)'},domain, h_filter, ...
                              v_filter,t_filter);
omghmagsq_timeseries = sqrt(FltInterp(flt_pos(:,:,2:4),hname,dname,{'u' ...
                    'v' 'cdiff(zvec,(2),''z'').^2+cdiff(zvec,(1),''z'').^2'},domain, h_filter, ...
                              v_filter,t_filter));

%omghx:
plot(TimeVector,mean(omghx_timeseries,1),'-','Color',cc(1,:),'LineWidth',2);
hold on;
ebar(TimeVector(7:8:end),mean(omghx_timeseries(:,7:8:end),1), ...
     std(omghx_timeseries(:,7:8:end),1),1/8,'-y','Color',cc(1,:),'LineWidth',2);

%omghy:
plot(TimeVector,mean(omghy_timeseries,1),'-','Color',cc(2,:),'LineWidth',2);
hold on;
ebar(TimeVector(6:8:end),mean(omghy_timeseries(:,6:8:end),1), ...
     std(omghy_timeseries(:,6:8:end),1),1/8,'-y','Color',cc(2,:),'LineWidth',2);

omghmagchk_timeseries = sqrt(omghx_timeseries.^2+omghy_timeseries.^2);

%Magnitude:
plot(TimeVector,mean(omghmag_timeseries,1),'-','Color',cc(3,:),'LineWidth',2);
hold on;
ebar(TimeVector(5:8:end),mean(omghmag_timeseries(:,5:8:end),1), ...
     std(omghmag_timeseries(:,5:8:end),1),1/8,'-y','Color',cc(3,:),'LineWidth',2);

%Magnitude check:
plot(TimeVector,mean(omghmagchk_timeseries,1),'-','Color',cc(4,:),'LineWidth',2);
hold on;
ebar(TimeVector(4:8:end),mean(omghmagchk_timeseries(:,4:8:end),1), ...
     std(omghmag_timeseries(:,4:8:end),1),1/8,'-y','Color',cc(4,:),'LineWidth',2);

%Magnitude check2:
plot(TimeVector,mean(omghmagsq_timeseries,1),'-','Color',cc(5,:),'LineWidth',2);
hold on;
ebar(TimeVector(3:8:end),mean(omghmagsq_timeseries(:,3:8:end),1), ...
     std(omghmag_timeseries(:,3:8:end),1),1/8,'-y','Color',cc(5,:),'LineWidth',2);

%Thermal wind buoyancy gradients:
omghyTW_timeseries = -DBDY_timeseries./(2*7.292*1e-5*sin(pi/180*flt_pos(:,:,6)));
omghxTW_timeseries = -DBDX_timeseries./(2*7.292*1e-5*sin(pi/180*flt_pos(:,:,6)));
plot(TimeVector,mean(omghxTW_timeseries,1),'--','Color',cc(1,:),'LineWidth',2);
hold on;
ebar(TimeVector(7:8:end),mean(omghxTW_timeseries(:,7:8:end),1), ...
     std(omghxTW_timeseries(:,7:8:end),1),1/8,'-y','Color',cc(1,:),'LineWidth',2);
plot(TimeVector,mean(omghyTW_timeseries,1),'--','Color',cc(2,:),'LineWidth',2);
hold on;
ebar(TimeVector(6:8:end),mean(omghyTW_timeseries(:,6:8:end),1), ...
     std(omghyTW_timeseries(:,6:8:end),1),1/8,'-y','Color',cc(2,:),'LineWidth',2);

%component in direction of buoyancy gradient:
bhatx = DBDX_timeseries./sqrt(DBDX_timeseries.^2+ ...
                              DBDY_timeseries.^2);
bhaty = DBDY_timeseries./sqrt(DBDX_timeseries.^2+ ...
                              DBDY_timeseries.^2);
omghDBmag = (omghx_timeseries.*bhatx+omghy_timeseries.*bhaty);

%Thermal wind buoyancy gradients:
plot(TimeVector,mean(omghDBmag,1),'--','Color',cc(3,:),'LineWidth',2);
hold on;
ebar(TimeVector(6:8:end),mean(omghDBmag(:,6:8:end),1), ...
     std(omghyTW_timeseries(:,6:8:end),1),1/8,'-y','Color',cc(3,:),'LineWidth',2);


text(284,0,'$|\omega_h|$','Interpreter','latex','FontSize',20,'Color',cc(3,:));
text(284,0,'$\omega_x$','Interpreter','latex','FontSize',20,'Color',cc(1,:));
text(284,0,'$\omega_y$','Interpreter','latex','FontSize',20,'Color',cc(2,:));
text(284,0,'$-\frac{1}{f}\frac{\partial b}{\partial x}$','Interpreter','latex','FontSize',20,'Color',cc(1,:));
text(284,0,'$-\frac{1}{f}\frac{\partial b}{\partial y}$','Interpreter','latex','FontSize',20,'Color',cc(2,:));
text(284,0,'$\sqrt{\left(\omega_x\right)^2+\left(\omega_y\right)^2}$','Interpreter','latex','FontSize',20,'Color',cc(4,:));
text(284,0,'$\sqrt{|\omega_h|^2}$','Interpreter','latex','FontSize',20,'Color',cc(5,:));
text(284,0,'$\frac{\omega_h \cdot \nabla_h b}{|\nabla_h b|}$', ...
     'Interpreter','latex','FontSize',20,'Color',cc(3,:));
plot(custom_tlims,[0 0],'-k','LineWidth',1);
xlabel('Day','Interpreter','latex','FontSize',25);
ylabel('$|\omega_h|$ ($10^{-7} s^{-2}$)','Interpreter','latex','FontSize',25);
set(gca,'YTick',[-4:2:4]*1e-2);
set(gca,'YTickLabel',{'-4','-2','0','2','4'});
set(gca,'YLim',[-4 4]*1e-2);
set(gca,'XLim',[282 custom_tlims(2)]);

% $$$ %Direction plots:
% $$$ fbar = 1.27e-5;
% $$$ DBDX = mean(DBDX_timeseries,1)/fbar;
% $$$ DBDY = mean(DBDY_timeseries,1)/fbar;
% $$$ omghx = mean(omghx_timeseries,1);
% $$$ omghy = mean(omghy_timeseries,1);
% $$$ omghhatx = omghx_timeseries./sqrt(omghx_timeseries.^2+ omghy_timeseries.^2);
% $$$ omghhaty = omghy_timeseries./sqrt(omghx_timeseries.^2+ omghy_timeseries.^2);
% $$$ DBmag = (DBDX_timeseries.*omghhatx+DBDY_timeseries.*omghhaty)/fbar;
% $$$ DBDBx = mean(DBmag.*omghhatx,1);
% $$$ DBDBy = mean(DBmag.*omghhaty,1);
% $$$ 
% $$$ fig3 = figure;
% $$$ set(fig3, 'Position', [314 38 1175 600]);
% $$$ tposvec = 2:4:tL;
% $$$ separator = 1.2e5;
% $$$ cenpos = -separator;
% $$$ scale = 0.4e5/mean(sqrt(DBDX.^2+DBDY.^2));
% $$$ for t = 1:length(tposvec)
% $$$     cenpos = cenpos+separator;
% $$$     plot([scale*DBDX(tposvec(t)) 0]+cenpos,[scale*DBDY(tposvec(t)) 0], ...
% $$$           '-','LineWidth',2,'Color',cc(3,:));
% $$$     hold on;
% $$$     plot([scale*omghx(tposvec(t)) 0]+cenpos,[scale*omghy(tposvec(t)) 0], ...
% $$$           '-','LineWidth',2,'Color',cc(2,:));
% $$$     plot([scale*DBDBx(tposvec(t)) 0]+cenpos,[scale*DBDBy(tposvec(t)) 0], ...
% $$$           '-','LineWidth',2,'Color',cc(1,:));
% $$$ end
% $$$ daspect([1 1 1]);
% $$$ axis([-separator/2 length(tposvec)*separator-separator/2 -separator/2 ...
% $$$       separator/2]);
% $$$ xlabel('Day','Interpreter','latex','FontSize',25);
% $$$ ylabel('$\nabla_h b$','Interpreter','latex','FontSize',25);
% $$$ set(gca,'YTick',[]);
% $$$ set(gca,'YTickLabel',[]);
% $$$ set(gca,'XTick',[0 separator 2*separator-50000 2*separator ...
% $$$                  2*separator+50000 3*separator:separator:((length(tposvec)-1)*separator)]);
% $$$ LabelList = cell(length(tposvec)+2,1);
% $$$ for t = 1:length(tposvec)
% $$$     if (t>3)
% $$$         shf = 2;
% $$$     elseif (t == 3)
% $$$         shf = 1;
% $$$     else
% $$$         shf = 0;
% $$$     end
% $$$     LabelList{t+shf} = num2str(TimeVector(tposvec(t)),4);
% $$$ end
% $$$ set(gca,'XTickLabel',LabelList);
% $$$ 

%%%%%%%%%w'b' and HCL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For HCL:
FGF = {'u' 'v' 'w' 'DBDX' 'DBDY' 'N2' 'b_diff' 'something'};
HCL = FGF;
HCL{end} = ['-2*(4).*(cdiff(x_rho,(1),''x'').*(4)+cdiff(x_rho,' ...
               '(2),''x'').*(5))-2*(5).*(cdiff(y_rho,(1),''y'').*(4)+cdiff(y_rho,(2),''y'').*(5))'];

filtkm = 75;
filtpts = floor(filtkm/3);
if (mod(filtpts,2)==0)
    filtpts=filtpts+1;
end
wbVarOp = {'w' 'rho' ['((1)-filter_field((1),' num2str(filtpts) ',''-s'')).*(-9.81/1025*((2)+1000)-' ...
                    'filter_field(-9.81/1025*((2)+1000),' num2str(filtpts) ',' ...
                    '''-s''))']};

cc = fliplr(jet(7));
cc = jet(7)
cc = cc(7:-1:4,:);
cc(5,:) = [0 0 1];
cc(4,:) = 0;
cc(3,:) = [0 0.5 0];
fig3 = figure;
set(fig3, 'Position', [314 38 1175 600]);

wb_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname,wbVarOp,domain, h_filter, v_filter,t_filter);

%wb:
plot(TimeVector,mean(wb_timeseries,1),'-','Color',cc(1,:),'LineWidth',2);
hold on;
ebar(TimeVector(7:8:end),mean(wb_timeseries(:,7:8:end),1), ...
     std(wb_timeseries(:,7:8:end),1),1/8,'-y','Color',cc(1,:),'LineWidth',2);

vtext(284,0,'$w''b''$','Interpreter','latex','FontSize',20,'Color',cc(1,:));
plot(custom_tlims,[0 0],'-k','LineWidth',1);
xlabel('Day','Interpreter','latex','FontSize',25);
ylabel('$w''b''$ ($10^{-7} m^{2}s^{-3}$)','Interpreter','latex','FontSize',25);
set(gca,'YTick',[-4:2:4]*1e-7);
set(gca,'YTickLabel',{'-4','-2','0','2','4'});
set(gca,'YLim',[-4 4]*1e-7);
set(gca,'XLim',[282 custom_tlims(2)]);


fig3 = figure;
set(fig3, 'Position', [314 38 1175 600]);
HCL_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname,HCL,domain, h_filter, v_filter,t_filter);
%HCL:
plot(TimeVector,mean(HCL_timeseries,1),'-','Color',cc(1,:),'LineWidth',2);
hold on;
ebar(TimeVector(7:8:end),mean(HCL_timeseries(:,7:8:end),1), ...
     std(HCL_timeseries(:,7:8:end),1),1/8,'-y','Color',cc(1,:),'LineWidth',2);

text(284,0,'$HCL$','Interpreter','latex','FontSize',20,'Color',cc(1,:));
plot(custom_tlims,[0 0],'-k','LineWidth',1);
xlabel('Day','Interpreter','latex','FontSize',25);
ylabel('$w''b''$ ($10^{-19} m^{2}s^{-3}$)','Interpreter','latex','FontSize',25);
set(gca,'YTick',[-10:5:10]*1e-19);
set(gca,'YTickLabel',{'-10','-5','0','5','10'});
set(gca,'YLim',[-10 10]*1e-19);
set(gca,'XLim',[282 custom_tlims(2)]);



%AkT:
hisname = '/mnt/Data1/ryan/UW_ROMS_EQ40/run6k1/ocean_his.nc';
slice = {domain(:,1),domain(:,2),domain(:,3),domain(:,4)}
AKt_timeseries = FltInterpORIG(flt_pos(:,:,2:4),hisname,dname,{'AKt'},domain,h_filter,v_filter,t_filter);

fig3 = figure;
set(fig3, 'Position', [314 38 1175 600]);
plot(TimeVector,mean(AKt_timeseries,1),'-','Color',cc(1,:),'LineWidth',2);
hold on;
ebar(TimeVector(7:8:end),mean(AKt_timeseries(:,7:8:end),1), ...
     std(AKt_timeseries(:,7:8:end),1),1/8,'-y','Color',cc(1,:),'LineWidth',2);
text(284,4e-5,'$AKt$','Interpreter','latex','FontSize',20,'Color',cc(1,:));
plot(custom_tlims,[0 0],'-k','LineWidth',1);
xlabel('Day','Interpreter','latex','FontSize',25);
ylabel('$w''b''$ ($10^{-19} m^{2}s^{-3}$)','Interpreter','latex','FontSize',25);
% $$$ set(gca,'YTick',[3:0.5:5.5]*1e-5);
% $$$ set(gca,'YTickLabel',{'3','3.5','4','4.5','5','5.5'});
% $$$ set(gca,'YLim',[3 5.5]*1e-5);
set(gca,'XLim',[282 custom_tlims(2)]);


%N2:
N2_timeseries = FltInterpORIG(flt_pos(:,:,2:4),hisname,dname,{'rho' ...
                   'cdiff(zvec,(1),''z'')'},domain,h_filter,v_filter,t_filter);

fig3 = figure;
set(fig3, 'Position', [314 38 1175 600]);
plot(TimeVector,mean(N2_timeseries,1),'-','Color',cc(1,:),'LineWidth',2);
hold on;
ebar(TimeVector(7:8:end),mean(N2_timeseries(:,7:8:end),1), ...
     std(N2_timeseries(:,7:8:end),1),1/8,'-y','Color',cc(1,:),'LineWidth',2);
text(284,4e-10,'$N2$','Interpreter','latex','FontSize',20,'Color',cc(1,:));
plot(custom_tlims,[0 0],'-k','LineWidth',1);
xlabel('Day','Interpreter','latex','FontSize',25);
ylabel('$w''b''$ ($10^{-19} m^{2}s^{-3}$)','Interpreter','latex','FontSize',25);
% $$$ set(gca,'YTick',[3:0.5:5.5]*1e-5);
% $$$ set(gca,'YTickLabel',{'3','3.5','4','4.5','5','5.5'});
% $$$ set(gca,'YLim',[3 5.5]*1e-5);
set(gca,'XLim',[282 custom_tlims(2)]);

%Ri:
Ri_timeseries = FltInterp(flt_pos(:,:,2:4),hname,dname,{'u' ...
                    'v' 'N2' '(3)./(cdiff(zvec,(1),''z'').^2+cdiff(zvec,(2),''z'').^2)'},domain,h_filter,v_filter,t_filter);

fig3 = figure;
set(fig3, 'Position', [314 38 1175 600]);
plot(TimeVector,mean(Ri_timeseries,1),'-','Color',cc(1,:),'LineWidth',2);
hold on;
ebar(TimeVector(7:8:end),mean(Ri_timeseries(:,7:8:end),1), ...
     std(Ri_timeseries(:,7:8:end),1),1/8,'-y','Color',cc(1,:),'LineWidth',2);
text(284,5,'$Ri$','Interpreter','latex','FontSize',20,'Color',cc(1,:));
plot(custom_tlims,[0 0],'-k','LineWidth',1);
xlabel('Day','Interpreter','latex','FontSize',25);
ylabel('$w''b''$ ($10^{-19} m^{2}s^{-3}$)','Interpreter','latex','FontSize',25);
set(gca,'YTick',[0:2.5:10]);
set(gca,'YTickLabel',{'0','2.5','5','7.5','10'});
set(gca,'YLim',[0 10]);
set(gca,'XLim',[282 custom_tlims(2)]);
