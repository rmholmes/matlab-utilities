function [SST,lon,lat] = plotTRMM(dvec,region)
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%
% This function plots 3-day averaged TRMM SST for the specified
% date in the specified region on the current axes.
%
% INPUTS:
%
% dvec = date vector for plotting, which is rounded to the nearest
% day (middle day of 3-day period).
%
% region = [lonW lonE latS latN] region to plot in.
%
%---------------------------------------------------------------------
%
% Dependencies; none. 
%
% Ryans ROMS Matlab and netcdf Utilities 22/10/14
%
%---------------------------------------------------------------------
%---------------------------------------------------------------------

%--------------TRMM information---------------------------------------
post = 'v4_d3d';
SSTwebsitebase = 'ftp://ftp.ssmi.com/tmi/bmaps_v04/';
fnSSTbase = '/mnt/Data1/ryan/TRMM/SST/';

%TRMM lon/lat:
lon_lin = 0.25*(1:1440)-0.125-360;
lat_lin = 0.25*(1:320)-40.125;
[lon,lat] = ndgrid(lon_lin,lat_lin);

[tmp indE] = min(abs(lon_lin-region(1)));
[tmp indW] = min(abs(lon_lin-region(2)));
[tmp indS] = min(abs(lat_lin-region(3)));
[tmp indN] = min(abs(lat_lin-region(4)));

lon = lon(indE:indW,indS:indN);
lat = lat(indE:indW,indS:indN);

%--------------Get timing information---------------------------------
year = dvec(1);
month = dvec(2);
day = round(dvec(3)+dvec(4)/24+dvec(5)/24/60+dvec(6)/24/60/60);

%add day for 3-day average:
day = day+1;
if (day>eomday(year,month))
    day = day-eomday(year,month);
    month = month+1;
end
if (month>12)
    year = year+1;
    month = 1;
end

yearstr = num2str(year);
if (month<10)
    monthstr = ['0' num2str(month)];
else
    monthstr = num2str(month);
end
if (day<10)
    daystr = ['0' num2str(day)];
else
    daystr = num2str(day);
end

%--------------Get TRMM data-------------------------------------------
fnTRMM = [fnSSTbase 'tmi_' yearstr monthstr daystr ...
          post];

%Download data:
if (exist(fnTRMM) ~= 2)
    fnEXTbase = [SSTwebsitebase 'y' yearstr '/m' monthstr '/'];
    fnEXT = [fnEXTbase 'tmi_' yearstr monthstr daystr post ...
             '.gz'];
    system(['wget ' fnEXT ' -O ' fnTRMM '.gz'],'-echo');
    system(['gunzip ' fnTRMM '.gz'],'-echo');
end

%--------------Plot----------------------------------------------------
[sst,wind11,wind37,vapor,cloud,rain] = ...
    read_tmi_averaged_v4(fnTRMM);
sst(abs(sst)>100) = NaN;
sst = sst(indE:indW,indS:indN);
pcolPlot(lon,lat,sst);

end

