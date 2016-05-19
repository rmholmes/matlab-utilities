
function [flt_timeseries] = FltInterp(flt_pos,hname,dname,VarOp, ...
                                      h_filter,v_filter,t_filter)
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%
% This function generates a time series of the chosen variable along
% the float tracks given in flt_pos.
%
% OUTPUTS:
%
% flt_timeseries = time series of interpolated values for each float
% (#flt x #time)
%
% INPUTS:
% hname = filename of history file (variables without underscores)
%
% dname = filename of diagnostics file (variables with underscores)
%
% flt_pos = time series of float positions
% ((flt)x(time)x(Xgrid,Ygrid,depth)) (positions must be at history
% field times)
%
% VarOp = GetVar field operation cell.
%
% slice = GetVar slice to extract spatial ({[xi xf],[yi yf],[zi zf]})
%
% h_filter = radius of circular horizontal box filter window. (2
% elements for two spatial filters)
%
% v_filter = # pts of vertical box filter window.
%
% t_filter = # pts of time box filter window.
%
%---------------------------------------------------------------------
%
% Dependencies; GetVar.
%
% Ryans ROMS Matlab and netcdf Utilities 22/7/13
%
%---------------------------------------------------------------------
%---------------------------------------------------------------------

%----options-----

%buffer outside max/min lon/lat/depth
xbuf = 5;
ybuf = 5;
zbuf = 2;

%Remove negative initial z-values data:
tmp = flt_pos(:,:,4);
tmp(tmp<0) = NaN;
flt_pos(:,:,4) = tmp;

%---setup-------

%Get file limits:
ncid = netcdf.open(hname,'NC_NOWRITE');
[tmp xLh] = netcdf.inqDim(ncid,netcdf.inqDimID(ncid,'xi_rho'));
[tmp yLh] = netcdf.inqDim(ncid,netcdf.inqDimID(ncid,'eta_rho'));
[tmp zLh] = netcdf.inqDim(ncid,netcdf.inqDimID(ncid,'s_rho'));

%----extract slice-------------
slice = cell(4);

timeF = flt_pos(1,:,1)';
tLF = length(timeF); %number of float times
fL = length(flt_pos(:,1,1));
timeH = ncread(hname,'ocean_time')/86400;
[tmp ti] = min(abs(timeH-flt_pos(1,1,1)));
[tmp tf] = min(abs(timeH-flt_pos(1,end,1)));
tL = tf-ti+1; 

%Find common times:
tvecF = zeros(0,0);
tvecH = zeros(0,0);
cnt = 1;
for t=ti:tf
    ind = find(timeF == timeH(t));
    if (length(ind)==1)
       tvecF(cnt) = ind;
       tvecH(cnt) = t;
       cnt = cnt+1;
    end
end
ti = min(tvecH);
tf = max(tvecH);
tL = tf-ti+1; %number of history times
tLa = length(tvecH); %number of common times
['Interpolating for ' num2str(tLa) ' common times, where there were ' ...
 num2str(tL) ' history times and ' num2str(tLF) 'float times']

xi = max(1,round(min(min(flt_pos(:,:,2))))-xbuf);
xf = min(xLh,round(max(max(flt_pos(:,:,2))))+xbuf);
xL = xf-xi+1;

yi = max(1,round(min(min(flt_pos(:,:,3))))-ybuf);
yf = min(yLh,round(max(max(flt_pos(:,:,3))))+ybuf);
yL = yf-yi+1;

zi = max(1,round(min(min(flt_pos(:,:,4))))-zbuf);
zf = min(zLh,round(max(max(flt_pos(:,:,4))))+zbuf);
zL = zf-zi+1;

slice = {[xi xf],[yi yf],[zi zf],[ti tf]};
% $$$ Var = GetVar(hname,dname,VarOp,slice);
% $$$ vard = length(size(Var));
% $$$ Var(isnan(Var)) = 0;
% $$$ 
% $$$ %Filter in vertical and horizontal:
% $$$ status = 'filtering variable...'
% $$$ if (vard == 4)
% $$$ for t = 1:tL
% $$$         Var(:,:,:,t) = filter_field(filter_field(Var(:,:,:,t),h_filter,['-' ...
% $$$                             's']),v_filter,'-v');
% $$$ end
% $$$ elseif (vard == 3)
% $$$ for t = 1:tL
% $$$         Var(:,:,t) = filter_field(Var(:,:,t),h_filter,'-s');
% $$$ end
% $$$ end
status = 'interpolating variable...'
 
%Deal with NaNs:
NaNs = isnan(flt_pos(:,:,9));
flt_pos(:,:,2) = Repl(flt_pos(:,:,2),NaN,xi);
flt_pos(:,:,3) = Repl(flt_pos(:,:,3),NaN,yi);
flt_pos(:,:,4) = Repl(flt_pos(:,:,4),NaN,zi);

%Interpolate onto float tracks:
flt_timeseries = NaN*zeros(fL,tLF);

% $$$ if (vard == 4)
 [X,Y,Z] = meshgrid(1:xL,1:yL,1:zL);
for tii = 1:tLa
    ['Doing ' num2str(tii) ' of ' num2str(tLa) '...']
    slice{4} = [tvecH(tii) tvecH(tii)];
    Var = GetVar(hname,dname,VarOp,slice);
    Var(isnan(Var)) = 0;
    flt_timeseries(:,tvecF(tii)) = interp3(X,Y,Z,permute(Var,[2 1 3 4]),flt_pos(:,tvecF(tii),2)-xi+1, ...
                                  flt_pos(:,tvecF(tii),3)-yi+1, flt_pos(:,tvecF(tii), ...
                                                      4)-zi+1,'spline');
end
% $$$ elseif (vard == 3)
% $$$ [X,Y] = meshgrid(1:xL,1:yL);
% $$$ for tii = 1:tLa
% $$$     flt_timeseries(:,tvecF(tii)) = interp2(X,Y,permute(Var(:,:,tvecH(tii)),[2 1 3]),flt_pos(:,tvecF(tii),2)-xi+1, ...
% $$$                                   flt_pos(:,tvecF(tii),3)-yi+1,'spline');
% $$$ end
% $$$ end

%put back in NaNs:
flt_timeseries(NaNs) = NaN;

%Filter in time:
flt_timeseries = filter_field(flt_timeseries,t_filter,'-t');
status = 'done!'

end
