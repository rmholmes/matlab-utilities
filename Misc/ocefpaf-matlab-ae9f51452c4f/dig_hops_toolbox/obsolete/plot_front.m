%
% plot_front.m
% purpose: plot a digitized fron in a SST image
% author: Filipe Fernandes
% date: 26-May-08
%
% obs: plot front only over a already plotted image
%      need a fermi image name loaded in the work space
%

%% date vector
 year=['20',arqi(4:5)]; day=datenum([arqi(9:10),'-',arqi(6:8),'-',year])-datenum(['00-00-',year]);

%% load front file
 [arq,paf] = uigetfile(['../EXPS_',num2str(year),'/',num2str(day),'/gsfm.*'],'Choose the front file'); load([paf,arq]);
%   eval(['load gsfm.',year,num2str(day)])

% plot points
plot(-gsfm(:,3),gsfm(:,2),'.','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')