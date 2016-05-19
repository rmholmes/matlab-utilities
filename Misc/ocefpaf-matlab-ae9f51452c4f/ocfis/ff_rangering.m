function [h] = ff_rangering(long,lat,range,lonlim,latlim,varargin)
% ff_rangering.m -> plot circles using a center and a range
%
%   call:  [output] = ff_rangering(input);
%
%  input:
%          long, lat -> latitude/longitude center
%              range -> radius
%     lonlim, latlim -> map box (needed to call projection)
%
% output:  h -> figure handle
%
% example: [h] = ff_rangering(long,lat,range,lonlim,latlim,varargin)
%
% need: m_map packages
%

%
% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     09-Jan-2010
% modified: 09-Jan-2010
%
% obs: basically it is a copy of m_range_ring.m from m_map package
%

m_proj('mercator','long',[lonlim(1) lonlim(2)],'lat',[latlim(1) latlim(2)],'on')

global MAP_VAR_LIST

pi180=pi/180; earth_radius=6378.137; n=72;

if length(varargin)>0 & ~ischar(varargin{1}),
 n=varargin{1};varargin(1)=[];
end;

c=range(:)'/earth_radius;

h=[];
for k=1:length(long),
  rlat=lat(k)*pi180;
  rlong=long(k)*pi180;
  if long(k)<MAP_VAR_LIST.longs(1), rlong=rlong+2*pi; end
  if long(k)>MAP_VAR_LIST.longs(2), rlong=rlong-2*pi; end

  x=sin([0:n-1]'/(n-1)*2*pi)*c;
  y=cos([0:n-1]'/(n-1)*2*pi)*c;
  on=ones(n,1);

  Y=(asin(on*cos(c)*sin(rlat) + (on*cos(rlat)*(sin(c)./c)).*y))/pi180;
  switch lat(k),
    case 90,
      X=(rlong+atan2(x,-y))/pi180;
    case -90,
      X=(rlong+atan2(x,y))/pi180;
    otherwise
      X=(rlong+atan2(x.*(on*sin(c)),on*(cos(rlat)*cos(c).*c) - (on*sin(rlat)*sin(c)).*y ) )/pi180;
  end

  nz=zeros(1,length(range(:)));
  X=X+cumsum([nz;diff(X)<-300]-[nz;diff(X)>300])*360;

  kk=find(X(1,:)~=X(end,:));
  X2=X(:,kk)+360;X2(X2>MAP_VAR_LIST.longs(2))=NaN;
  X3=X(:,kk)-360;X3(X3<MAP_VAR_LIST.longs(1))=NaN;

  [XX,YY]=m_ll2xy([X,X2,X3],[Y,Y(:,kk),Y(:,kk)],'clip','on');

%% get rid of 2-point lines (these are probably clipped lines spanning the window)
  fk=finite(XX(:));
  st=find(diff(fk)==1)+1;
  ed=find(diff(fk)==-1);
  if length(st)<length(ed), st=[1;st]; end
  if length(ed)<length(st), ed=[ed;length(fk)]; end
  k=find((ed-st)==1);
  XX(st(k))=NaN;
  [XX,YY]=m_xy2ll(XX,YY);
  h=[h;line(XX,YY,varargin{:},'tag','m_range_ring')];
end

if nargout==0,
 clear h
end