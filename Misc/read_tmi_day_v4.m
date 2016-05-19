function [time,sst,wind11,wind37,vapor,cloud,rain]=read_tmi_day_v4(data_file)
% [time,sst,wind11,wind37,vapor,cloud,rain]=read_tmi_day_v4(data_file);
%
% this subroutine reads compressed or uncompressed RSS TMI daily byte maps
% (version-4 data released Sepetember 2006)
%
% input arguments:
% data_file = the full path and name of the uncompressed data file
%
% the function returns these products:
% time  = time of observation in fractional hours GMT
% sst   = sea surface temperature in deg C
% wind11= wind speed derived using 11 GHz channel in meters/second
% wind37= wind speed derived using 37 GHz channel in meters/second
% vapor = atmospheric water vapor in millimeters
% cloud = liquid cloud water in millimeters
% rain  = rain rate in millimeters/hour
%
% longitude is 0.25*xdim-0.125
% latitude is 0.25*ydim-40.125
%
% For detailed data description, see 
% http://www.remss.com/tmi/tmi_description.html
% Remote Sensing Systems
% support@remss.com


scale=[0.1,.15,.2,.2,.3,.01,.1];
offset=[0.,-3.,0.,0.,0.,0.,0.];
xdim=1440;ydim=320;tdim=2;numvar=7;
mapsiz=xdim*ydim*tdim;

if ~exist(data_file,'file')
   disp(['file not found: ' data_file]);
   time=[];sst=[];wind11=[];wind37=[];vapor=[];cloud=[];rain=[];
   return;
end;

if ~isempty(regexp(data_file,'.gz', 'once'))
    data_file=char(gunzip(data_file));
end

fid=fopen(data_file,'rb');
data=fread(fid,mapsiz*numvar,'uint8');
fclose(fid);
disp(data_file);
map=reshape(data,[xdim ydim numvar tdim]);

for i=1:numvar
    tmp=map(:,:,i,:);
    ia=find(tmp<=250);tmp(ia)=tmp(ia)*scale(i)+offset(i);
    map(:,:,i,:)=tmp;
end;

time   = squeeze(map(:,:,1,:));
sst    = squeeze(map(:,:,2,:));
wind11 = squeeze(map(:,:,3,:));
wind37 = squeeze(map(:,:,4,:));
vapor  = squeeze(map(:,:,5,:));
cloud  = squeeze(map(:,:,6,:));
rain   = squeeze(map(:,:,7,:));

return;
end
