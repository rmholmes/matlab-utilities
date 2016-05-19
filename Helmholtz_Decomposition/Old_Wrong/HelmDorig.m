
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

function [psi,phi] = HelmDorig(lon,lat,u,v)

%%%%%%%%%%Parameters:
xres_multiplyer = 1;
yres_multiplyer = 1;

%%%%%%%%%%%%%%%%%%%%%%

%Get size:
[xL,yL] = size(lon);

%Get x,y:
CF = pi*6371000/180; %conversion factor from lat/lon to distance.
x = CF*(lon+132).*cos(pi/180*lat);
y = CF*lat;

%Get Resolutions:
xres = xres_multiplyer*abs(x(round(xL/2),round(yL/2))-x(round(xL/2)+1,round(yL/2)));
yres = yres_multiplyer*abs(y(round(xL/2),round(yL/2))-y(round(xL/2),round(yL/2)+1));

%Create regular grid:
[Y,X] = meshgrid(min(min(y)):yres:max(max(y)),min(min(x)):xres:max(max(x)));
[m,n] = size(X);
Lat = Y/CF;
Lon = X/CF./cos(pi/180*Lat)-132;

%Construct Dirichlett Laplacian Matrix:
one = ones((m-2)*(n-2),1);
D = spdiags([1/xres/xres.*one (-2*(1/xres/xres+1/yres/yres)).*one ...
             1/xres/xres.*one],[-1 0 1],m-2,m-2);
M = spdiags([1/yres/yres.*one],[0],m-2,m-2);

LAP = zeros(2,2);
LAP = sparse(LAP);
for i=1:(n-2)
    for k=1:(n-2)
    	if (i == k)
            LAP(((i-1)*(m-2)+1):(i*(m-2)),((i-1)*(m-2)+1):(i*(m-2))) = D;
        elseif( i == k-1 | i == k+1)
            LAP(((i-1)*(m-2)+1):(i*(m-2)),((k-1)*(m-2)+1):(k*(m-2))) = M;
        end
    end
end

%Interpolate variables:
U = interp2(lon',lat',u',Lon,Lat,'spline');
V = interp2(lon',lat',v',Lon,Lat,'spline');
Div = cdiff(X,U,'x')+cdiff(Y,V,'y');

%%%Solve for PHI first:
%Get Dirichlett boundary conditions with perimeter integral:
PHI = zeros(m,n);
for j=2:n %Left edge
    PHI(1,j) = yres*V(1,j-1)+PHI(1,j-1);
end
for i=2:m %Top edge
    PHI(i,n) = xres*U(i-1,n)+PHI(i-1,n);
end
for j=(n-1):-1:1 %Right edge
    PHI(m,j) = -V(m,j)*yres+PHI(m,j+1);
end
for i=(m-1):-1:2 %Bottom edge
    PHI(i,1) = -U(i,1)*xres+PHI(i+1,1);
end

%Make RHS matrices:
%PHI:
RPHI = zeros((m-2)*(n-2),1);
cnt = 1;
for j=2:(n-1)
    RPHI(((j-2)*(m-2)+1):((j-1)*(m-2))) = Div(2:(m-1),j);
end
RPHI(1:(m-2)) = RPHI(1:(m-2)) - PHI(2:(m-1),1)/yres/yres;
RPHI(((m-2)*(n-3)+1):((m-2)*(n-2))) = RPHI(((m-2)*(n-3)+1):((m-2)*(n-2)))-PHI(2:(m-1),n)/yres/yres;
RPHI(1:(m-2):((m-2)*(n-3)+1)) = RPHI(1:(m-2):((m-2)*(n-3)+1))-PHI(1,2:(n-1))'/xres/xres;
RPHI((m-2):(m-2):((m-2)*(n-2))) = RPHI((m-2):(m-2):((m-2)*(n- 2)))-PHI(m,2:(n-1))'/xres/xres;

%Solve:
PHIlist = LAP\RPHI;

%Reorder:
for j=2:(n-1)
    PHI(2:(m-1),j) = PHIlist(((j-2)*(m-2)+1):((j-1)*(m-2)));
end

%Interpolate back to orginal grid:
phi = interp2(X',Y',PHI',x,y,'spline');

%Take this part off the velocity:
U = U - cdiff(X,PHI,'x');
V = V - cdiff(Y,PHI,'y');
Curl = cdiff(X,V,'x')-cdiff(Y,U,'y');

%%%Now solve for PSI:
%Get Dirichlett boundary conditions with perimeter integral:
PSI = zeros(m,n);
for j=2:n %Left edge
    PSI(1,j) = yres*U(1,j-1)+PSI(1,j-1);
end
for i=2:m %Top edge
    PSI(i,n) = -xres*V(i-1,n)+PSI(i-1,n);
end
for j=(n-1):-1:1 %Right edge
    PSI(m,j) = -U(m,j)*yres+PSI(m,j+1);
end
for i=(m-1):-1:2 %Bottom edge
    PSI(i,1) = V(i,1)*xres+PSI(i+1,1);
end

%Make RHS matrices:
%PSI:
RPSI = zeros((m-2)*(n-2),1);
cnt = 1;
for j=2:(n-1)
    RPSI(((j-2)*(m-2)+1):((j-1)*(m-2))) = -Curl(2:(m-1),j);
end
RPSI(1:(m-2)) = RPSI(1:(m-2)) - PSI(2:(m-1),1)/yres/yres;
RPSI(((m-2)*(n-3)+1):((m-2)*(n-2))) = RPSI(((m-2)*(n-3)+1):((m-2)*(n-2)))-PSI(2:(m-1),n)/yres/yres;
RPSI(1:(m-2):((m-2)*(n-3)+1)) = RPSI(1:(m-2):((m-2)*(n-3)+1))-PSI(1,2:(n-1))'/xres/xres;
RPSI((m-2):(m-2):((m-2)*(n-2))) = RPSI((m-2):(m-2):((m-2)*(n-2)))-PSI(m,2:(n-1))'/xres/xres;

%Solve:
PSIlist = LAP\RPSI;

%Reorder:
for j=2:(n-1)
    PSI(2:(m-1),j) = PSIlist(((j-2)*(m-2)+1):((j-1)*(m-2)));
end

%Interpolate back to orginal grid:
psi = interp2(X',Y',PSI',x,y,'spline');

% $$$ 
% $$$ %%%%%%%%%%%BOTH AT SAME TIME
% $$$ %Get Dirichlett boundary conditions with perimeter integral:
% $$$ PHI = zeros(m,n);
% $$$ PSI = zeros(m,n);
% $$$ for j=2:n %Left edge
% $$$     PHI(1,j) = yres*V(1,j-1)+PHI(1,j-1);
% $$$     PSI(1,j) = yres*U(1,j-1)+PSI(1,j-1);
% $$$ end
% $$$ for i=2:m %Top edge
% $$$     PHI(i,n) = xres*U(i-1,n)+PHI(i-1,n);
% $$$     PSI(i,n) = -xres*V(i-1,n)+PSI(i-1,n);
% $$$ end
% $$$ for j=(n-1):-1:1 %Right edge
% $$$     PHI(m,j) = -V(m,j)*yres+PHI(m,j+1);
% $$$     PSI(m,j) = -U(m,j)*yres+PSI(m,j+1);
% $$$ end
% $$$ for i=(m-1):-1:2 %Bottom edge
% $$$     PHI(i,1) = -U(i,1)*xres+PHI(i+1,1);
% $$$     PSI(i,1) = V(i,1)*xres+PSI(i+1,1);
% $$$ end
% $$$ 
% $$$ %Make RHS matrices:
% $$$ %PSI:
% $$$ RPSI = zeros((m-2)*(n-2),1);
% $$$ cnt = 1;
% $$$ for j=2:(n-1)
% $$$     RPSI(((j-2)*(m-2)+1):((j-1)*(m-2))) = -Curl(2:(m-1),j);
% $$$ end
% $$$ RPSI(1:(m-2)) = RPSI(1:(m-2)) - PSI(2:(m-1),1)/yres/yres;
% $$$ RPSI(((m-2)*(n-3)+1):((m-2)*(n-2))) = RPSI(((m-2)*(n-3)+1):((m-2)*(n-2)))-PSI(2:(m-1),n)/yres/yres;
% $$$ RPSI(1:(m-2):((m-2)*(n-3)+1)) = RPSI(1:(m-2):((m-2)*(n-3)+1))-PSI(1,2:(n-1))'/xres/xres;
% $$$ RPSI((m-2):(m-2):((m-2)*(n-2))) = RPSI((m-2):(m-2):((m-2)*(n-2)))-PSI(m,2:(n-1))'/xres/xres;
% $$$ %PHI:
% $$$ RPHI = zeros((m-2)*(n-2),1);
% $$$ cnt = 1;
% $$$ for j=2:(n-1)
% $$$     RPHI(((j-2)*(m-2)+1):((j-1)*(m-2))) = Div(2:(m-1),j);
% $$$ end
% $$$ RPHI(1:(m-2)) = RPHI(1:(m-2)) - PHI(2:(m-1),1)/yres/yres;
% $$$ RPHI(((m-2)*(n-3)+1):((m-2)*(n-2))) = RPHI(((m-2)*(n-3)+1):((m-2)*(n-2)))-PHI(2:(m-1),n)/yres/yres;
% $$$ RPHI(1:(m-2):((m-2)*(n-3)+1)) = RPHI(1:(m-2):((m-2)*(n-3)+1))-PHI(1,2:(n-1))'/xres/xres;
% $$$ RPHI((m-2):(m-2):((m-2)*(n-2))) = RPHI((m-2):(m-2):((m-2)*(n- 2)))-PHI(m,2:(n-1))'/xres/xres;
% $$$ 
% $$$ %Solve:
% $$$ PSIlist = LAP\RPSI;
% $$$ PHIlist = LAP\RPHI;
% $$$ 
% $$$ %Reorder:
% $$$ for j=2:(n-1)
% $$$     PSI(2:(m-1),j) = PSIlist(((j-2)*(m-2)+1):((j-1)*(m-2)));
% $$$     PHI(2:(m-1),j) = PHIlist(((j-2)*(m-2)+1):((j-1)*(m-2)));
% $$$ end
% $$$ 
% $$$ %Interpolate back to orginal grid:
% $$$ psi = interp2(X',Y',PSI',x,y,'spline');
% $$$ phi = interp2(X',Y',PHI',x,y,'spline');

end


