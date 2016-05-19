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

% $$$ %Now solve for psi on a smaller grid:
% $$$ X = X(2:(end-1),2:(end-1));
% $$$ Y = Y(2:(end-1),2:(end-1));
% $$$ Urot = Urot(2:(end-1),2:(end-1));
% $$$ Vrot = Vrot(2:(end-1),2:(end-1));
% $$$ CURL = Dsimp(dx,Vrot,'x')-Dsimp(dy,Urot,'y');
% $$$ [xL,yL] = size(X);
% $$$ 
% $$$ %Construct Dirichlett Laplacian Matrix:
% $$$ one = ones((xL-2)*(yL-2),1);
% $$$ D = spdiags([1/dx/dx.*one (-2*(1/dx/dx+1/dy/dy)).*one ...
% $$$              1/dx/dx.*one],[-1 0 1],xL-2,xL-2);
% $$$ M = spdiags([1/dy/dy.*one],[0],xL-2,xL-2);
% $$$ 
% $$$ LAP = zeros(2,2);
% $$$ LAP = sparse(LAP);
% $$$ for i=1:(yL-2)
% $$$     for k=1:(yL-2)
% $$$     	if (i == k)
% $$$             LAP(((i-1)*(xL-2)+1):(i*(xL-2)),((i-1)*(xL-2)+1):(i*(xL-2))) = D;
% $$$         elseif( i == k-1 | i == k+1)
% $$$             LAP(((i-1)*(xL-2)+1):(i*(xL-2)),((k-1)*(xL-2)+1):(k*(xL-2))) = M;
% $$$         end
% $$$     end
% $$$ end
% $$$ 
% $$$ %Perimeter integral:
% $$$ PSI = zeros(xL,yL);
% $$$ for j=2:yL %Left edge
% $$$     PSI(1,j) = dy*Urot(1,j-1)+PSI(1,j-1);
% $$$ end
% $$$ for i=2:xL %Top edge
% $$$     PSI(i,yL) = -dx*Vrot(i-1,yL)+PSI(i-1,yL);
% $$$ end
% $$$ for j=(yL-1):-1:1 %Right edge
% $$$     PSI(xL,j) = -Urot(xL,j)*dy+PSI(xL,j+1);
% $$$ end
% $$$ for i=(xL-1):-1:1 %Bottom edge
% $$$     PSI(i,1) = Vrot(i,1)*dx+PSI(i+1,1);
% $$$ end
% $$$ 
% $$$ %Finally, make the sizes correct:
% $$$ PHI = PHI(2:(end-1),2:(end-1));
% $$$ LAT = Y/CF;
% $$$ LON = X/CF./cos(pi/180*LAT)-132;
% $$$ 
% $$$ Udiv = Dsimp(dx,PHI,'x');
% $$$ Vdiv = Dsimp(dy,PHI,'y');
% $$$ Urot = Dsimp(dy,PSI,'y');
% $$$ Vrot = -Dsimp(dx,PSI,'x');
% $$$ Uc = Udiv+Urot;
% $$$ Vc = Udiv+Urot;

%Check Divergence and Curl:
Divdiv = Dsimp(dx,Udiv,'x')+Dsimp(dy,Vdiv,'y');
Divrot = Dsimp(dx,Urot,'x')+Dsimp(dy,Vrot,'y');
Rotdiv = Dsimp(dx,Vdiv,'x')-Dsimp(dy,Udiv,'y');
Rotrot = Dsimp(dx,Vrot,'x')-Dsimp(dy,Urot,'y');

caxs = [-1 1]*1e-5;
figure;
set(gcf,'Position',[172 203 1681 741]);
subplot(2,3,1);
pcolPlot(LON,LAT,DIV);
caxis(caxs);
subplot(2,3,4);
pcolPlot(LON,LAT,CURL);
caxis(caxs);
subplot(2,3,2);
pcolPlot(LON,LAT,Divdiv);
caxis(caxs);
subplot(2,3,5);
pcolPlot(LON,LAT,Rotdiv);
caxis(caxs);
subplot(2,3,3);
pcolPlot(LON,LAT,Divrot);
caxis(caxs);
subplot(2,3,6);
pcolPlot(LON,LAT,DIV-Divdiv);
caxis(caxs);

caxs = [-1.5 1.5];
figure;
set(gcf,'Position',[172 203 1681 741]);
subplot(2,3,1);
pcolPlot(lon,lat,u);
caxis(caxs);
subplot(2,3,4);
pcolPlot(lon,lat,v);
caxis(caxs);
subplot(2,3,2);
pcolPlot(LON,LAT,Udiv);
caxis(caxs);
subplot(2,3,5);
pcolPlot(LON,LAT,Vdiv);
caxis(caxs);
subplot(2,3,3);
pcolPlot(LON,LAT,Urot);
caxis(caxs);
subplot(2,3,6);
pcolPlot(LON,LAT,Vrot);
caxis(caxs);