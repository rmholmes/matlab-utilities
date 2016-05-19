function [lon,lat,sst] = plot_modis(filename,lonlim,latlim,plotting)

PI=hdfinfo(filename);

% And write it into a structure
pin=[];
lims = length(PI.Attributes)-3;
for k=1:lims,
    nm=PI.Attributes(k).Name;nm(nm==' ')='_';
    if isstr(PI.Attributes(k).Value),
        pin=setfield(pin,nm,PI.Attributes(k).Value);
    else
        pin=setfield(pin,nm,double(PI.Attributes(k).Value));
    end
end;

% lon/lat of grid corners
lon = [pin.Westernmost_Longitude:pin.Longitude_Step:pin.Easternmost_Longitude];
lat = [pin.Northernmost_Latitude:-pin.Latitude_Step:pin.Southernmost_Latitude];

% $$$ % Get the indices needed for the area of interest
% $$$ [mn,ilt]=min(abs(lat-max(latlim)));
% $$$ [mn,ilg]=min(abs(lon-min(lonlim)));
% $$$ ltlm=fix(diff(latlim)/pin.Latitude_Step);
% $$$ lglm=fix(diff(lonlim)/pin.Longitude_Step);

% load the subset of data needed for the map limits given
sst=hdfread(filename,'l3m_data');%,'Index',{[ilt ilg],[],[ltlm lglm]});
sst = pin.Slope*double(sst)+pin.Intercept;
sst(sst>35) = NaN;
sst(sst<5) = NaN;

%flip around -180:
sst2 = sst;
inds = lon>0;
sst2(:,1:length(find(inds))) = sst(:,inds);
sst2(:,(length(find(inds))+1):end) = sst(:,1:length(find(inds)));
sst = sst2;
lon2 = lon;
lon2(1:length(find(inds))) = lon(inds);
lon2((length(find(inds))+1):end) = lon(1:length(find(inds)));
lon = lon2;
lon(lon>0) = lon(lon>0)-360;

[tmp ind1] = min(abs(lon-min(lonlim)));
[tmp ind2] = min(abs(lon-max(lonlim)));
[tmp ind3] = min(abs(lat-min(latlim)));
[tmp ind4] = min(abs(lat-max(latlim)));
lon = lon(ind1:ind2);
lat = lat(ind4:ind3);
sst = sst(ind4:ind3,ind1:ind2)';
[lon,lat] = ndgrid(lon,lat);

if (plotting)
pcolPlot(lon,lat,sst);
end
end
