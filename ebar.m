
function ebar(X,Y,SD,W,type,varargin)
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%
% This function plots error bars on the current plot.
%
% INPUTS:
% (X,Y) = points to center Error bars on.
%
% SD = half width of error bar (vector)
%
% W = length of side bars on error bars (scalar)
%
% type = error bars in x ('-x') or in y ('-y')
%
% varargin = any plotting/line specs.
%
%---------------------------------------------------------------------
%
% Dependencies; none. 
%
% Ryans ROMS Matlab and netcdf Utilities 17/7/13
%
%---------------------------------------------------------------------
%---------------------------------------------------------------------

switch type

  case '-y'
    Xe = [];Ye=[]
    for i=1:length(X)
        Xe = [Xe X(i) X(i) NaN X(i)-W/2 X(i)+W/2 NaN X(i)-W/2 X(i)+W/2 ...
             NaN];
        Ye = [Ye Y(i)-SD(i) Y(i)+SD(i) NaN Y(i)+SD(i) Y(i)+SD(i) NaN ...
             Y(i)-SD(i) Y(i)-SD(i) NaN];
% $$$         plot([X(i)-W/2 X(i)+W/2],[Y(i) Y(i)]+SD(i),'-',varargin{1: end});
% $$$         plot([X(i)-W/2 X(i)+W/2],[Y(i) Y(i)]-SD(i),'-',varargin{1: ...
% $$$                             end});
% $$$         plot([X(i) X(i)],[Y(i)-SD(i) Y(i)+SD(i)],'-',varargin{1:end});
    end
    plot(Xe,Ye,'-',varargin{1:end});
    
  case '-x'
    Xe = [];Ye=[];
    for i=1:length(X)
        Xe = [Xe X(i)-SD(i) X(i)+SD(i) NaN X(i)+SD(i) X(i)+SD(i) NaN ...
              X(i)-SD(i) X(i)-SD(i) NaN];
        Ye = [Ye Y(i) Y(i) NaN Y(i)-W/2 Y(i)+W/2 NaN Y(i)-W/2 Y(i)+W/2 ...
              NaN];
% $$$         plot([X(i)-SD(i) X(i)+SD(i)],[Y(i) Y(i)],'-',varargin{1:end});
% $$$         plot([X(i) X(i)]+SD(i),[Y(i)-W/2 Y(i)+W/2],'-',varargin{1: end});
% $$$         plot([X(i) X(i)]-SD(i),[Y(i)-W/2 Y(i)+W/2],'-',varargin{1: end});
    end
    plot(Xe,Ye,'-',varargin{1:end});
end


end

    