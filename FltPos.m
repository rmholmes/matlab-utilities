
function [flt_pos] = FltPos(fname,flts,tlims,Rst)
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%
% This function generates the standard flt_pos matrix for a set of
% floats (flts) from the file fname for input into FltInterp
%
% OUTPUTS:
%
% flt_pos = time series of float positions
% ((flt)x(time)x(Xgrid,Ygrid,depth)) (positions must be at history
% field times)
%
% INPUTS:
%
% flts = vector of float numbers to consider
%
% fname = float output .nc file.
%
% tlims = [tst tend] -> time limits in days to output flt_pos
% on. If 0 then all.
%
% Rst = 0 -> don't restrict the time dimension to only non NaNs.
% Rst = anything else -> restrict the time dimension to only times
% when none of the floats have NaNs.
%
%---------------------------------------------------------------------
%
% Dependencies; none
%
% Ryans ROMS Matlab and netcdf Utilities 5/8/13
%
%---------------------------------------------------------------------
%---------------------------------------------------------------------

%-------------file
%information--------------------------------------
ncid = netcdf.open(fname,'NC_NOWRITE');
ftime = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'ocean_time'),'double')/86400;
XgridID = netcdf.inqVarID(ncid,'Xgrid');
YgridID = netcdf.inqVarID(ncid,'Ygrid');
ZgridID = netcdf.inqVarID(ncid,'Zgrid');
depthID = netcdf.inqVarID(ncid,'depth');
lonID = netcdf.inqVarID(ncid,'lon');
latID = netcdf.inqVarID(ncid,'lat');
tempID = netcdf.inqVarID(ncid,'temp');
saltID = netcdf.inqVarID(ncid,'salt');
drID = netcdf.inqDimID(ncid,'drifter');
try
tID = netcdf.inqDimID(ncid,'ftime');
catch exception
tID = netcdf.inqDimID(ncid,'ocean_time');
end
[tmp nflt] = netcdf.inqDim(ncid,drID);
% $$$ [tmp tL] = netcdf.inqDim(ncid,tID);


%-------------Time limits-----------------------------------------
if (tlims == 0)
    mint = 1;
    maxt = length(ftime);
    tL = maxt;
else
    [tmp mint] = min(abs(ftime-tlims(1)));
    [tmp maxt] = min(abs(ftime-tlims(2)));
    tL = maxt-mint+1;
end

if (mint > maxt)
    temp = mint;
    mint = maxt;
    maxt = temp;
    tL = maxt-mint+1;
end

%---------------Get flt_pos -------------------------------------

if (flts == 0)
    flts = 1:nflt;
end
% $$$ if (length(flts) == nflt)
flt_pos = zeros(nflt,tL,10);
flt_pos(:,:,1) = repmat(ftime(mint:maxt)',[nflt 1]);
flt_pos(:,:,2) = netcdf.getVar(ncid,XgridID,[0 mint-1],[nflt ...
                    tL],'double')+1;
flt_pos(:,:,3) = netcdf.getVar(ncid,YgridID,[0 mint-1],[nflt ...
                    tL],'double')+1;
flt_pos(:,:,4) = netcdf.getVar(ncid,ZgridID,[0 mint-1],[nflt ...
                    tL],'double')+0.5;
flt_pos(:,:,5) = netcdf.getVar(ncid,depthID,[0 mint-1],[nflt ...
                    tL],'double');
flt_pos(:,:,6) = netcdf.getVar(ncid,lonID,[0 mint-1],[nflt ...
                    tL],'double');
flt_pos(:,:,7) = netcdf.getVar(ncid,latID,[0 mint-1],[nflt ...
                    tL],'double');
flt_pos(:,:,9) = netcdf.getVar(ncid,tempID,[0 mint-1],[nflt ...
                    tL],'double');
flt_pos(:,:,10) = netcdf.getVar(ncid,saltID,[0 mint-1],[nflt ...
                    tL],'double');
flt_pos = flt_pos(flts,:,:);

% $$$ for f = 1:length(flts)
% $$$     fltsn = flts(f);
% $$$     flt_pos(f,:,1) = ftime(mint:maxt)';
% $$$     flt_pos(f,:,2) = netcdf.getVar(ncid,XgridID,[fltsn-1 mint-1],[1 ...
% $$$                         tL],'double')'+1;
% $$$     flt_pos(f,:,3) = netcdf.getVar(ncid,YgridID,[fltsn-1 mint-1],[1 ...
% $$$                         tL],'double')'+1;
% $$$     flt_pos(f,:,4) = netcdf.getVar(ncid,ZgridID,[fltsn-1 mint-1],[1 ...
% $$$                         tL],'double')'+0.5;
% $$$     flt_pos(f,:,5) = netcdf.getVar(ncid,depthID,[fltsn-1 mint-1],[1 ...
% $$$                         tL],'double')';
% $$$     flt_pos(f,:,6) = netcdf.getVar(ncid,lonID,[fltsn-1 mint-1],[1 ...
% $$$                         tL],'double')';
% $$$     flt_pos(f,:,7) = netcdf.getVar(ncid,latID,[fltsn-1 mint-1],[1 ...
% $$$                         tL],'double')';
% $$$ % $$$     flt_pos(f,:,8) = netcdf.getVar(ncid,rhoID,[fltsn-1 mint-1],[1 ...
% $$$ % $$$                         tL],'double')';
% $$$     flt_pos(f,:,9) = netcdf.getVar(ncid,tempID,[fltsn-1 mint-1],[1 ...
% $$$                         tL],'double')';
% $$$     flt_pos(f,:,10) = netcdf.getVar(ncid,saltID,[fltsn-1 mint-1],[1 ...
% $$$                         tL],'double')';
% $$$ end
netcdf.close(ncid);
%Make bad values NaNs:
flt_pos(abs(flt_pos)>1e10) = NaN;

%Restrict time period if some of these floats have invalid values:
if (Rst ~= 0)
swi = 0;
mintaft = mint;
maxtaft = maxt;
for t = 2:length(flt_pos(1,:,1))
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
end

%Make sure time goes forwards:
if (flt_pos(1,1,1)>flt_pos(1,end,1))
    flt_pos = flipdim(flt_pos,2);
end
end