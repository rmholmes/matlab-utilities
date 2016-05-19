
%This function calculates the streamfunction \psi and potential \phi
%in the Helmholtz decomposition of a 2D vector field (u,v) defined on
%the longitude/latitude grid lon,lat.
%
% u = \nabla\psi \times k + \nabla \phi
%
% \nabla^2 \psi = -\nabla\times u
%
% \nabla^2 \phi = \nabla\cdot u
%
% The resolution used on the interpolated rectangular grid is
% determined by a multiplyer (parameter below) times the resolution
% in m's of the center of the grid.
%
%NOTE: THIS DECOMPOSITION I don't think it is unique when
%considering open boundary conditions! I have chosen to calculate
%the divergent part first and then take this off the result and
%calculate the streamfunction. This appears to work pretty well but
%may be unfairly biasing to psi or something?

close all;
clear all;
load filenames;

%Testing:
fname = BASF_his;
time = 300;
box = [-180 -120 -20 25];
zlvl = 50;

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

%Get size:
[xL,yL] = size(lon);

%Get x,y:
CF = pi*6371000/180; %conversion factor from lat/lon to distance.
x = CF*(lon+132).*cos(pi/180*lat);
y = CF*lat;

%Generate regular grid which is all interpolated:
dx = 10000;
dy = 20000;
xmin = max(min(x,[],1),[],2)+dx;
xmax = min(max(x,[],1),[],2)-dx;
ymin = max(min(y,[],2),[],1)+dy;
ymax = min(max(y,[],2),[],1)-dy;

[X,Y] = ndgrid(xmin:dx:xmax,ymin:dy:ymax);
LAT = Y/CF;
LON = X/CF./cos(pi/180*LAT)-132;
[xL,yL] = size(X);

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

%Interpolate variables:
U = interp2(lon',lat',u',LON,LAT,'spline');
V = interp2(lon',lat',v',LON,LAT,'spline');
DIV = zeros(xL-2,yL-2);
CURL = zeros(xL-2,yL-2);
DIV = Dsimp(dx,U,'x')+Dsimp(dy,V,'y');
CURL = Dsimp(dx,V,'x')-Dsimp(dy,U,'y');

%Solve for Phi first:
PHI = zeros(xL,yL);

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
Udiv = Dsimp(dx,PHI,'x');
Vdiv = Dsimp(dy,PHI,'y');
Urot = U-Udiv;
Vrot = V-Vdiv;

Urot = Dsimp(dy,PHI,'y');
Vrot = -Dsimp(dx,PHI,'x');
Udiv = U-Urot;
Vdiv = V-Vrot;

%Check Divergence and Curl:
Divdiv = Dsimp(dx,Udiv,'x')+Dsimp(dy,Vdiv,'y');
Divrot = Dsimp(dx,Urot,'x')+Dsimp(dy,Vrot,'y');
Rotdiv = Dsimp(dx,Vdiv,'x')-Dsimp(dy,Udiv,'y');
Rotrot = Dsimp(dx,Vrot,'x')-Dsimp(dy,Urot,'y');

caxs = [-1 1]*1e-5;
figure;
set(gcf,'Position',[172 203 1681 741]);
subplot(2,3,1);
pcolPlot(X,Y,DIV);
caxis(caxs);
subplot(2,3,4);
pcolPlot(X,Y,CURL);
caxis(caxs);
subplot(2,3,2);
pcolPlot(X,Y,Divdiv);
caxis(caxs);
subplot(2,3,5);
pcolPlot(X,Y,Rotdiv);
caxis(caxs);
subplot(2,3,3);
pcolPlot(X,Y,Divrot);
caxis(caxs);
subplot(2,3,6);
pcolPlot(X,Y,Rotrot);
caxis(caxs);

caxs = [-1.5 1.5];
figure;
set(gcf,'Position',[172 203 1681 741]);
subplot(2,3,1);
pcolPlot(X,Y,U);
caxis(caxs);
subplot(2,3,4);
pcolPlot(X,Y,V);
caxis(caxs);
subplot(2,3,2);
pcolPlot(X,Y,Udiv);
caxis(caxs);
subplot(2,3,5);
pcolPlot(X,Y,Vdiv);
caxis(caxs);
subplot(2,3,3);
pcolPlot(X,Y,Urot);
caxis(caxs);
subplot(2,3,6);
pcolPlot(X,Y,Vrot);
caxis(caxs);