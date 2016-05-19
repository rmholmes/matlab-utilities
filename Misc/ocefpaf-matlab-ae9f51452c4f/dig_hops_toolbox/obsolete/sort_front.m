function sort_front (front_file)
%
% [front_file] = sort_front (front_file)
%
% input: unsorted front file after digitalization
%        ASCII -> 'gsfm.yyyyddd.dig'
%
% output: sorted (by longitude) front file
%        ASCII -> 'gsfm.yyyyddd.dig'

%
% sort_front.m
% purpose: Sort the digitized front for FM
% author: Filipe Fernandes
% date: 18-Aug-08
%
% obs:
%

eval(['load ' front_file])

 %% lon, lat from file
 lon = -gsfm(:,3); lat = gsfm(:,2);

 %% sort by longitude
 [xdata,i] = sort(lon,'ascend');
 ydata=lat(i);

 %% rewrite file
 day = gsfm(:,1);

 %% Save gsfm.yyyyddd file, day in Julians
 aux=[day ydata -xdata];
 [m,n]=size(aux);
 str2='%8.4f '; str1='%6.0f ';
 fmt=[repmat(str1,1,1),repmat(str2,1,2)];
 fid1=fopen(front_file,'w');
 fprintf(fid1,[fmt, '\n'],aux');
 fclose(fid1); clear fid1

 %% show the first 10 lines of the file
 disp(front_file)
 eval(['!head ',front_file])