
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

function [PHI,PSI,LON,LAT,dx,dy] = HelmDReg(lon,lat,u,v)

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
% $$$ DIV = zeros(xL-2,yL-2);
% $$$ CURL = zeros(xL-2,yL-2);
DIV = Dsimp(dx,U,'x')+Dsimp(dy,V,'y');
% $$$ CURL = Dsimp(dx,V,'x')-Dsimp(dy,U,'y');

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

%Now solve for psi on a smaller grid:
X = X(2:(end-1),2:(end-1));
Y = Y(2:(end-1),2:(end-1));
Urot = Urot(2:(end-1),2:(end-1));
Vrot = Vrot(2:(end-1),2:(end-1));
CURL = Dsimp(dx,Vrot,'x')-Dsimp(dy,Urot,'y');
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

%Finally, make the sizes correct:
PHI = PHI(2:(end-1),2:(end-1));
LAT = Y/CF;
LON = X/CF./cos(pi/180*LAT)-132;
end
