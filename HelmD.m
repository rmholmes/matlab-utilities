
%This function calculates the streamfunction and velocity potential
%by converting to a regular staggered grid and outputing on this
%grid.

function [PSI,PHI,LON_p,LAT_p,LON_r,LAT_r,dx,dy] = HelmD(lon,lat,u,v)
% $$$           ,Udiv,Vdiv,Urot,Vrot, ...
% $$$           LON_u,LAT_u,LON_v,LAT_v

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
%Get LON and LAT coordinates:
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
F = TriScatteredInterp(lon(:),lat(:),u(:));
U = F(LON_u,LAT_u);
F = TriScatteredInterp(lon(:),lat(:),v(:));
V = F(LON_v,LAT_v);

%Calculate divergence:
DIV = Dsimp2(dx,U,'x')+Dsimp2(dy,V,'y');

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
% $$$ Udiv = Udiv(:,2:(end-1));
% $$$ Vdiv = Vdiv(2:(end-1),:);
Urot = Urot(:,2:(end-1));
Vrot = Vrot(2:(end-1),:);

%Calculate curl:
tmp = Dsimp2(dx,Vrot,'x');
tmp2 = Dsimp2(dy,Urot,'y');
Currot = tmp(:,2:(end-1))-tmp2(2:(end-1),:);

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

% $$$ LON_u = LON_u(3:(end-2),3:(end-2));
% $$$ LAT_u = LAT_u(3:(end-2),3:(end-2));
% $$$ LON_v = LON_v(3:(end-2),3:(end-2));
% $$$ LAT_v = LAT_v(3:(end-2),3:(end-2));

%This function takes a central derivative on a staggered grid of a
%field F defined on a regular grid.
function dF = Dsimp2(dx,F,type)

    if (strcmp(type,'x'))
        dF = (F(2:end,:)-F(1:(end-1),:))/dx;
    else
        dF = (F(:,2:end)-F(:,1:(end-1)))/dx;
    end
end

end

