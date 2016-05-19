function PlotVectorField(lon,lat,xvec,yvec,scale,Escale,pts,varargin);
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%
% This function plots a horizontal vector field over the current
% axis 
%
% INPUTS
%
% (lon,lat,xvec,yvec) = longitude,latitude,vector component in x
% direction, vector component in y direction.
%
% scale = the multiplicative factor to convert the vector length to
% meters
% Escale = size of example vector plotted in top left (units of xvec/yvec)
% pts = [xpts,ypts] number of vectors to draw in x/y directions
%
% varargin contains optional input (pairs) for marker/line/colors etc. as
% follows:
%
% '-vcol','b' - set the vector line color to blue
% '-mcol','b' - set the vector line color to blue
% '-mksz',15  - set the marker size to 15
% '-lwid',2   - set the vector line width to 2. 
%
%---------------------------------------------------------------------
%
% Dependencies; none. 
%
% Ryans ROMS Matlab and netcdf Utilities 17/7/13
%
%---------------------------------------------------------------------
%---------------------------------------------------------------------
    
mksz = 12;
lwid = 1;
vcol = 'k';
mcol = 'k';
lg = length(varargin);
if (lg>1)
    if (mod(lg,2)==0)
        for i=1:(lg/2)
            if (strcmp(varargin{2*i-1},'-vcol'))
                vcol = varargin{2*i};
            elseif (strcmp(varargin{2*i-1},'-mcol'))
                mcol = varargin{2*i};
            elseif (strcmp(varargin{2*i-1},'-lwid'))
                lwid = varargin{2*i};
            elseif (strcmp(varargin{2*i-1},'-mksz'))
                mksz = varargin{2*i};
            else
                warning('Bad input to varargin');
            end
        end
    else
        warning('Bad input to varargin');
    end
end

    
%Restrict to current axis limits:
lons = get(gca,'xlim');
lonmin = lons(1);
lonmax = lons(2);
lats = get(gca,'ylim');
latmin = lats(1);
latmax = lats(2);
[tmp latmini] = min(abs(lat(1,:)-latmin));
[tmp latmaxi] = min(abs(lat(1,:)-latmax));
[tmp lonmini] = min(abs(lon(:,round((latmaxi+latmini)/2))-lonmin));
[tmp lonmaxi] = min(abs(lon(:,round((latmaxi+latmini)/2))-lonmax));
lon = lon(lonmini:lonmaxi,latmini:latmaxi);
lat = lat(lonmini:lonmaxi,latmini:latmaxi);
xvec = xvec(lonmini:lonmaxi,latmini:latmaxi);
yvec = yvec(lonmini:lonmaxi,latmini:latmaxi);



ranX = max(max(lon))-min(min(lon));
dx = ranX/pts(1);
ranY = max(max(lat))-min(min(lat));
dy = ranY/pts(2);
[X,Y] = meshgrid((min(min(lon))+dx):dx:(max(max(lon))-dx), ...
                 (min(min(lat))+dy):dy:(max(max(lat))-dy));
uint = interp2(lon',lat',xvec',X,Y,'linear');
vint = interp2(lon',lat',yvec',X,Y,'linear');
hold on;
plot(X(~isnan(uint)),Y(~isnan(uint)),'.','Color',mcol,'MarkerSize',mksz);
for i=1:length(X(:,1))
    for j=1:length(Y(1,:))
        if (~isnan(uint(i,j)))
        plot([X(i,j) X(i,j)+uint(i,j)*scale],[Y(i,j) Y(i,j)+vint(i, ...
                                                          j)* ...
                            scale],'-','Color',vcol,'LineWidth',lwid);
        end
    end
end
plot([lonmin+0.1 lonmin+0.1+scale*Escale],[latmax-0.1 latmax-0.1],'-','Color',vcol,'LineWidth',lwid);
plot([lonmin+0.1],[latmax-0.1],'.','Color',mcol,'MarkerSize',mksz);

end