function [h,lon,lat]=plt_front(frontfile,color,marker,size)
%
% [h,lon,lat]=plt_front(frontfile,color,marker,size)
%
% use:
% [h,lon,lat]=plt_front('../EXPS_2008/245/gsfm.2008245.dig','k','.',5)
% lon,lat-> front position
% h -> graphic handle

%
% plt_front.m
% purpose: plot a digitized fron in a SST image
% author: Filipe Fernandes
% date: 12-Oct-08
%
% obs: function version
%

%% load front file
 load(frontfile);

% plot points
h=plot(-gsfm(:,3),gsfm(:,2),'Color',color,'Marker',marker,'LineStyle','none','MarkerSize',size);
lon=-gsfm(:,3); lat=gsfm(:,2);