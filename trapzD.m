function Int = trapzD(x,F,lower,upper,dim)
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%
% This function performs a trapezoidal integration but with the
% definite integration limits lower and upper.
%
% INPUTS:
%
% x  -> as for trapz.
% F -> as for trapz.
%
% lower -> lower limit of integration.
%
% upper -> upper limit of integration.
%
% dim -> dimension to integrate along.
%---------------------------------------------------------------------
%
% Dependencies; none. 
%
% Ryans ROMS Matlab and netcdf Utilities 11/11/13
%
%---------------------------------------------------------------------
%---------------------------------------------------------------------

%Check limits:
if (lower < min(x) | upper > max(x))
% $$$     warning('Bad limits! x doesn''t cover limits');
    if (upper > x(end))
        upper = x(end);
    end
    if (lower < x(1))
        lower = x(1);
    end
end
%Find the location of limits:
[tmp lowind] = min(abs(x-lower));
if (x(lowind)>lower)
    lowind = lowind-1;
end
if (x(lowind)>=upper)
    lowind = length(x)-1;
end

[tmp upind] = min(abs(x-upper));
if (x(upind)>=upper)
    upind = upind-1;
end
%Integrate between rounded grid points:
if (dim == 1)
    Int = trapz(x((lowind+1):upind),F((lowind+1):upind,:,:,:),dim);
elseif (dim == 2)
    Int = trapz(x((lowind+1):upind),F(:,(lowind+1):upind,:,:),dim);
elseif (dim == 3)
    Int = trapz(x((lowind+1):upind),F(:,:,(lowind+1):upind,:),dim);
elseif (dim == 4)
    Int = trapz(x((lowind+1):upind),F(:,:,:,(lowind+1):upind),dim);
end

%Add interpolated portions:
ldxP = x(lowind+1)-lower;
ldxM = lower-x(lowind);
ldx  = x(lowind+1)-x(lowind);
udxM = upper-x(upind);
udxP = x(upind+1)-upper;
udx  = x(upind+1)-x(upind);

if (dim == 1)
    
    Lint = F(lowind,:,:,:)+ldxM/ldx*(F(lowind+1,:,:,:)-F(lowind,:,: ...
                                                      ,:));
    Int = Int + (F(lowind+1,:,:,:)+Lint)/2*ldxP;
    
    Uint = F(upind,:,:,:)+udxM/udx*(F(upind+1,:,:,:)-F(upind,:,: ...
                                                      ,:));
    Int = Int + (F(upind,:,:,:)+Uint)/2*udxM;    
    
elseif (dim == 2)

    Lint = F(:,lowind,:,:)+ldxM/ldx*(F(:,lowind+1,:,:)-F(:,lowind,: ,:));
    Int = Int + (F(:,lowind+1,:,:)+Lint)/2*ldxP;
    
    Uint = F(:,upind,:,:)+udxM/udx*(F(:,upind+1,:,:)-F(:,upind,: ,:));
    Int = Int + (F(:,upind,:,:)+Uint)/2*udxM;    

elseif (dim == 3)

    Lint = F(:,:,lowind,:)+ldxM/ldx*(F(:,:,lowind+1,:)-F(:,:,lowind,:));
    Int = Int + (F(:,:,lowind+1,:)+Lint)/2*ldxP;
    
    Uint = F(:,:,upind,:)+udxM/udx*(F(:,:,upind+1,:)-F(:,:,upind,:));
    Int = Int + (F(:,:,upind,:)+Uint)/2*udxM;    

elseif (dim == 4)

    Lint = F(:,:,:,lowind)+ldxM/ldx*(F(:,:,:,lowind+1)-F(:,:,:,lowind));
    Int = Int + (F(:,:,:,lowind+1)+Lint)/2*ldxP;
    
    Uint = F(:,:,:,upind)+udxM/udx*(F(:,:,:,upind+1)-F(:,:,:,upind));
    Int = Int + (F(:,:,:,upind)+Uint)/2*udxM;    

end

end
