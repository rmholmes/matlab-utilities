function [gtime] = oasp_greg(jhours)
%
% oasp_greg.m - Converts Julian hours from OASP to Gregorian calendar.
%
% use:  [gtime] = oasp_greg(jhours);
% input:
%    jhours from OASP
%
% output:
%    calendar date [yyyy mo da hr mi sec]. 
%
% example:
%    [1998,01,11]=oasp_greg(859296)
%

% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     05-Oct-2009
% modified: 05-Oct-2009
%
% obs: based on RPStuff gregorian.m
%

      jd=jhours/24;

% Add 0.2 milliseconds before Gregorian calculation to prevent
% roundoff error resulting from math operations on time 
% from occasionally representing midnight as 
% (for example) [1990 11 30 23 59 59.99...] instead of [1990 12 1 0 0 0]);
% If adding a 0.2 ms to time (each time you go back and forth between 
% Julian and Gregorian) bothers you more than the inconvenient representation
% of Gregorian time at midnight you can comment this line out...

      jd=jd+2.e-9;    

%      if you want Julian Days to start at noon...
%      h=rem(jd,1)*24+12;
%      i=(h >= 24);
%      jd(i)=jd(i)+1;
%      h(i)=h(i)-24;

      secs=rem(jd,1)*24*3600;

      j = floor(jd) + 693902; %-1721119;
      in = 4*j -1;
      y = floor(in/146097);
      j = in - 146097*y;
      in = floor(j/4);
      in = 4*in +3;
      j = floor(in/1461);
      d = floor(((in - 1461*j) +4)/4);
      in = 5*d -3;
      m = floor(in/153);
      d = floor(((in - 153*m) +5)/5);
      y = y*100 +j;
      mo=m-9;
      yr=y+1;
      i=(m<10);
      mo(i)=m(i)+3;
      yr(i)=y(i);
      [hour,min,sec]=s2hms(secs);
      gtime=[yr(:) mo(:) d(:) hour(:) min(:) sec(:)];