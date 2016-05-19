%
% plot_rings.m
% purpose: plot digitized ring file on SST image
% author: Filipe Fernandes
% date: 16-Sep-08
%
% obs: image has to be already load with filename still in the workspace

%% date vector
 year=['20',arqi(4:5)]; day=datenum([arqi(9:10),'-',arqi(6:8),'-',year])-datenum(['00-00-',year]);

%% load front and ring files
%   name=[year,num2str(day)];
 [arq,paf] = uigetfile(['../EXPS_',num2str(year),'/',num2str(day),'/ring*.dig'],'Choose the ring file');
 fid = fopen([paf,arq], 'r');
 tline = fgetl(fid); x=str2num(tline);

 for i=1:x
  tline = fgetl(fid);
   if ~ischar(tline), break, end
  ring(i,:)=sscanf(tline,'%*s %d %d %g %g %d %d %d %d',[1,inf]);
 end

[a,b]=size(ring);
for i=1:a; h(i)=my_range_ring(-ring(i,4),ring(i,3),ring(i,6),[-85 -45],[20 55]); set(h(i),'Color','k'); end
plot(-ring(:,4),ring(:,3),'kx');