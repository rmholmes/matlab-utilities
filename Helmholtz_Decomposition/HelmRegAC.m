close all;
clear all;
load filenames;

%Testing:
fname = BASF_his;
time = 100;
box = [-160 -140 -20 25];
zlvl = 40;

%Get SSH and real velocity:
lon = ncread(fname,'lon_rho');
lat = ncread(fname,'lat_rho');
f = ncread(fname,'f');
[tmp ltmn] = min(abs(lat(1,:)-box(3)));
[tmp ltmx] = min(abs(lat(1,:)-box(4)));
[tmp lnmn] = min(abs(lon(:,round((ltmn+ltmx)/2))-box(1)));
[tmp lnmx] = min(abs(lon(:,round((ltmn+ltmx)/2))-box(2)));

SSH = GetVar(fname,0,{'zeta'},{[lnmn lnmx],[ltmn ltmx],[zlvl zlvl],[time ...
                    time]});
u = GetVar(fname,0,{'u'},{[lnmn-1 lnmx+1],[ltmn-1 ltmx+1],[zlvl zlvl],[time ...
                    time]});
v = GetVar(fname,0,{'v'},{[lnmn-1 lnmx+1],[ltmn-1 ltmx+1],[zlvl zlvl],[time ...
                    time]});
u = u(2:(end-1),2:(end-1));
v = v(2:(end-1),2:(end-1));
lon = lon(lnmn:lnmx,ltmn:ltmx);
lat = lat(lnmn:lnmx,ltmn:ltmx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Get x,y:
CF = pi*6371000/180; %conversion factor from lat/lon to distance.
x = CF*(lon+132).*cos(pi/180*lat);
y = CF*lat;

%grid spacings:
dx = min(min(x(2:end,:)-x(1:(end-1),:)));
dy = min(min(y(:,2:end)-y(:,1:(end-1))));

%Generate regular grid which is all interpolated:
xmin = max(min(x,[],1),[],2)+dx;
xmax = min(max(x,[],1),[],2)-dx;
ymin = max(min(y,[],2),[],1)+dy;
ymax = min(max(y,[],2),[],1)-dy;
xmax = floor((xmax-xmin)/dx)*dx+xmin;
ymax = floor((ymax-ymin)/dy)*dy+ymin;

[X_r,Y_r] = ndgrid(xmin:dx:xmax,ymin:dy:ymax);
[X_u,Y_u] = ndgrid((xmin-dx/2):dx:(xmax+dx/2),ymin:dy:ymax);
[X_v,Y_v] = ndgrid(xmin:dx:xmax,(ymin-dy/2):dy:(ymax+dy/2));
[X_p,Y_p] = ndgrid((xmin-dx/2):dx:(xmax+dx/2),(ymin-dy/2):dy:(ymax+ ...
                                                  dy/2));
LAT_r = Y_r/CF;
LON_r = X_r/CF./cos(pi/180*LAT_r)-132;
LAT_u = Y_u/CF;
LON_u = X_u/CF./cos(pi/180*LAT_u)-132;
LAT_v = Y_v/CF;
LON_v = X_v/CF./cos(pi/180*LAT_v)-132;
LAT_p = Y_p/CF;
LON_p = X_p/CF./cos(pi/180*LAT_p)-132;
[xL,yL] = size(X_r);

%Interpolate velocities:
U = interp2(lon',lat',u',LON_u,LAT_u,'spline');
V = interp2(lon',lat',v',LON_v,LAT_v,'spline');

%Calculate divergence and curl:
DIV = Dsimp2(dx,U,'x')+Dsimp2(dy,V,'y');
tmp = Dsimp2(dx,V,'x');
tmp2 = Dsimp2(dy,U,'y');
CURL = tmp(:,2:(end-1))-tmp2(2:(end-1),:);

%Solve for phi first:
PHI = zeros(xL,yL);

%Construct Dirichlett Laplacian Matrix:
one = ones((xL-2)*(yL-2),1);
D = spdiags([1/dx/dx.*one (-2*(1/dx/dx+1/dy/dy)).*one ...
             1/dx/dx.*one],[-1 0 1],xL-2,xL-2);
M = spdiags([1/dy/dy.*one],[0],xL-2,xL-2);

LAP = zeros(2,2);
LAP = sparse(LAP);
for i=1:(yL-2)
    for k=1:(yL-2)
    	if (i == k)
            LAP(((i-1)*(xL-2)+1):(i*(xL-2)),((i-1)*(xL-2)+1):(i*(xL-2))) = D;
        elseif( i == k-1 | i == k+1)
            LAP(((i-1)*(xL-2)+1):(i*(xL-2)),((k-1)*(xL-2)+1):(k*(xL-2))) = M;
        end
    end
end

%RHS:
RHSPHI = zeros((xL-2)*(yL-2),1);
cnt = 1;
for j=2:(yL-1)
    RHSPHI(((j-2)*(xL-2)+1):((j-1)*(xL-2))) = DIV(2:(xL-1),j);
end

%Solve:
PHIlist = LAP\RHSPHI;

%Reorder:
for j=2:(yL-1)
    PHI(2:(xL-1),j) = PHIlist(((j-2)*(xL-2)+1):((j-1)*(xL-2)));
end

%Calculate Divergent and rotational velocities:
Udiv = Dsimp2(dx,PHI,'x');
Vdiv = Dsimp2(dy,PHI,'y');
Urot = U(2:(end-1),:)-Udiv;
Vrot = V(:,2:(end-1))-Vdiv;
Udiv = Udiv(:,2:(end-1));
Vdiv = Vdiv(2:(end-1),:);
Urot = Urot(:,2:(end-1));
Vrot = Vrot(2:(end-1),:);

%Check Divergence and Curl:
Divdiv = Dsimp2(dx,Udiv,'x')+Dsimp2(dy,Vdiv,'y');
Divrot = Dsimp2(dx,Urot,'x')+Dsimp2(dy,Vrot,'y');
tmp = Dsimp2(dx,Vdiv,'x');
tmp2 = Dsimp2(dy,Udiv,'y');
Curdiv = tmp(:,2:(end-1))-tmp2(2:(end-1),:);
tmp = Dsimp2(dx,Vrot,'x');
tmp2 = Dsimp2(dy,Urot,'y');
Currot = tmp(:,2:(end-1))-tmp2(2:(end-1),:);

LON_rM = LON_r(2:(end-1),2:(end-1));
LAT_rM = LAT_r(2:(end-1),2:(end-1));
LON_pM = LON_p(2:(end-1),2:(end-1));
LAT_pM = LAT_p(2:(end-1),2:(end-1));
LON_pM2 = LON_p(3:(end-2),3:(end-2));
LAT_pM2 = LAT_p(3:(end-2),3:(end-2));

% $$$ caxs = [-1 1]*1e-5;
% $$$ figure;
% $$$ set(gcf,'Position',[172 203 1681 741]);
% $$$ subplot(2,3,1);
% $$$ pcolPlot(LON_r,LAT_r,DIV);
% $$$ caxis(caxs);
% $$$ subplot(2,3,4);
% $$$ pcolPlot(LON_pM,LAT_pM,CURL);
% $$$ caxis(caxs);
% $$$ subplot(2,3,2);
% $$$ pcolPlot(LON_rM,LAT_rM,Divdiv);
% $$$ caxis(caxs);
% $$$ subplot(2,3,5);
% $$$ pcolPlot(LON_pM2,LAT_pM2,Curdiv);
% $$$ caxis(caxs);
% $$$ subplot(2,3,3);
% $$$ pcolPlot(LON_rM,LAT_rM,Divrot);
% $$$ caxis(caxs);
% $$$ subplot(2,3,6);
% $$$ pcolPlot(LON_pM2,LAT_pM2,Currot);
% $$$ caxis(caxs);

PSI = zeros(xL,yL);

%Now solve for psi on a smaller grid:
[xL,yL] = size(Currot);
Urot = Urot(2:(end-1),:);
Vrot = Vrot(:,2:(end-1));

%Construct Dirichlett Laplacian Matrix:
one = ones((xL-2)*(yL-2),1);
D = spdiags([1/dx/dx.*one (-2*(1/dx/dx+1/dy/dy)).*one ...
             1/dx/dx.*one],[-1 0 1],xL-2,xL-2);
M = spdiags([1/dy/dy.*one],[0],xL-2,xL-2);

LAP = zeros(2,2);
LAP = sparse(LAP);
for i=1:(yL-2)
    for k=1:(yL-2)
    	if (i == k)
            LAP(((i-1)*(xL-2)+1):(i*(xL-2)),((i-1)*(xL-2)+1):(i*(xL-2))) = D;
        elseif( i == k-1 | i == k+1)
            LAP(((i-1)*(xL-2)+1):(i*(xL-2)),((k-1)*(xL-2)+1):(k*(xL-2))) = M;
        end
    end
end

%Perimeter integral:
PSI = zeros(xL,yL);
for j=2:yL %Left edge
    PSI(1,j) = dy*Urot(1,j-1)+PSI(1,j-1);
end
for i=2:xL %Top edge
    PSI(i,yL) = -dx*Vrot(i-1,yL)+PSI(i-1,yL);
end
for j=(yL-1):-1:1 %Right edge
    PSI(xL,j) = -Urot(xL,j)*dy+PSI(xL,j+1);
end
for i=(xL-1):-1:1 %Bottom edge
    PSI(i,1) = Vrot(i,1)*dx+PSI(i+1,1);
end

%Make RHS matrices:
%PSI:
RPSI = zeros((xL-2)*(yL-2),1);
cnt = 1;
for j=2:(yL-1)
    RPSI(((j-2)*(xL-2)+1):((j-1)*(xL-2))) = -Currot(2:(xL-1),j);
end
RPSI(1:(xL-2)) = RPSI(1:(xL-2)) - PSI(2:(xL-1),1)/dy/dy;
RPSI(((xL-2)*(yL-3)+1):((xL-2)*(yL-2))) = RPSI(((xL-2)*(yL-3)+1):((xL-2)*(yL-2)))-PSI(2:(xL-1),yL)/dy/dy;
RPSI(1:(xL-2):((xL-2)*(yL-3)+1)) = RPSI(1:(xL-2):((xL-2)*(yL-3)+1))-PSI(1,2:(yL-1))'/dx/dx;
RPSI((xL-2):(xL-2):((xL-2)*(yL-2))) = RPSI((xL-2):(xL-2):((xL-2)*(yL-2)))-PSI(xL,2:(yL-1))'/dx/dx;

%Solve:
PSIlist = LAP\RPSI;

%Reorder:
for j=2:(yL-1)
    PSI(2:(xL-1),j) = PSIlist(((j-2)*(xL-2)+1):((j-1)*(xL-2)));
end

%Organise out's:
PHI = PHI(2:(end-1),2:(end-1));
LON_r = LON_r(2:(end-1),2:(end-1));
LAT_r = LAT_r(2:(end-1),2:(end-1));
LON_p = LON_p(3:(end-2),3:(end-2));
LAT_p = LAT_p(3:(end-2),3:(end-2));

LON_u = LON_u(3:(end-2),3:(end-2));
LAT_u = LAT_u(3:(end-2),3:(end-2));
LON_v = LON_v(3:(end-2),3:(end-2));
LAT_v = LAT_v(3:(end-2),3:(end-2));

%Checks:
tmp = Dsimp2(dx,PHI,'x');
Udiv = tmp(:,2:(end-1));
Urot = Dsimp2(dy,PSI,'y');
Uc = Udiv+Urot;
tmp = Dsimp2(dy,PHI,'y');
Vdiv = tmp(2:(end-1),:);
Vrot = -Dsimp2(dx,PSI,'x');
Vc = Vdiv+Vrot;

%Plot velocities:
caxs = [-1.5 1.5];
figure;
set(gcf,'Position',[296 187 1219 741]);
subplot(2,2,1);
pcolPlot(lon,lat,u);
caxis(caxs);
axis(box);title('u');
subplot(2,2,2);
pcolPlot(lon,lat,v);
caxis(caxs);
axis(box);title('v');
subplot(2,2,3);
pcolPlot(LON_u,LAT_u,Uc);
caxis(caxs);
axis(box);title('$Dx\phi+Dy\psi$');
subplot(2,2,4);
pcolPlot(LON_v,LAT_v,Vc);
caxis(caxs);
axis(box);title('$Dy\phi-Dx\psi$');
figure;
set(gcf,'Position',[296 187 1219 741]);
subplot(2,2,1);
pcolPlot(LON_u,LAT_u,Uc);
caxis(caxs);
axis(box);title('u');
subplot(2,2,2);
pcolPlot(LON_v,LAT_v,Vc);
caxis(caxs);
axis(box);title('v');
subplot(2,2,3);
pcolPlot(lon,lat,u);
caxis(caxs);
axis(box);title('$Dx\phi+Dy\psi$');
subplot(2,2,4);
pcolPlot(lon,lat,v);
caxis(caxs);
axis(box);title('$Dy\phi-Dx\psi$');

figure;
set(gcf,'Position',[296 187 1219 741]);
subplot(2,2,1);
pcolPlot(LON_u,LAT_u,Udiv);
caxis(caxs);
axis(box);title('Udiv');
subplot(2,2,2);
pcolPlot(LON_v,LAT_v,Vdiv);
caxis(caxs);
axis(box);title('Vdiv');
subplot(2,2,3);
pcolPlot(LON_u,LAT_u,Urot);
caxis(caxs);
axis(box);title('Urot');
subplot(2,2,4);
pcolPlot(LON_v,LAT_v,Vrot);
caxis(caxs);
axis(box);title('Vrot');

%Check Divs and Curls:
Divdiv = Dsimp2(dx,Udiv,'x')+Dsimp2(dy,Vdiv,'y');
Divrot = Dsimp2(dx,Urot,'x')+Dsimp2(dy,Vrot,'y');
tmp = Dsimp2(dx,Vdiv,'x');
tmp2 = Dsimp2(dy,Udiv,'y');
Curdiv = tmp(:,2:(end-1))-tmp2(2:(end-1),:);
tmp = Dsimp2(dx,Vrot,'x');
tmp2 = Dsimp2(dy,Urot,'y');
Currot = tmp(:,2:(end-1))-tmp2(2:(end-1),:);

LON_rM = LON_r(2:(end-1),2:(end-1));
LAT_rM = LAT_r(2:(end-1),2:(end-1));
LON_pM = LON_p(2:(end-1),2:(end-1));
LAT_pM = LAT_p(2:(end-1),2:(end-1));
LON_pM2 = LON_p(3:(end-2),3:(end-2));
LAT_pM2 = LAT_p(3:(end-2),3:(end-2));

caxs = [-1 1]*1e-5;
figure;
set(gcf,'Position',[172 203 1681 741]);
subplot(2,3,1);
pcolPlot(LON_r,LAT_r,DIV(2:(end-1),2:(end-1)));
caxis(caxs);
axis(box);title('Orig. Divergence');
subplot(2,3,4);
pcolPlot(LON_p,LAT_p,CURL(2:(end-1),2:(end-1)));
caxis(caxs);
axis(box);title('Orig. Curl');
subplot(2,3,2);
pcolPlot(LON_rM,LAT_rM,Divdiv);
caxis(caxs);
axis(box);title('Divergence U\_div');
subplot(2,3,5);
pcolPlot(LON_pM2,LAT_pM2,Curdiv(2:(end-1),2:(end-1)));
caxis(caxs);
axis(box);title('Curl U\_div');
subplot(2,3,3);
pcolPlot(LON_rM,LAT_rM,Divrot);
caxis(caxs);
axis(box);title('Divergence U\_rot');
subplot(2,3,6);
pcolPlot(LON_pM2,LAT_pM2,Currot(2:(end-1),2:(end-1)));
caxis(caxs);
axis(box);title('Curl U\_rot');
