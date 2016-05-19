function int = gtrapz(x,y,x1,x2)
% This is a generalized trapz integration function to integrate a
% function over an interval (x1,x2) whoose limits don't
% neccessarily lie on an x element.

    [x,I] = sort(x);
    y = y(I);
    if (x1<min(x))
        x1 = min(x);
    end
    if (x2>max(x))
        x2 = max(x);
    end
   
    if (x2<x1)
        tmp = x1;
        x1 = x2;
        x2 = tmp;        
        mult = -1;
    else
        mult = 1;
    end
    
%find first and last elements inside limits:
inds = find(x>=x1,1,'first');
indl = find(x<=x2,1,'last');

int = 0;
if (indl>inds)
%Do central part of integral using trapz:
int = int+trapz(x(inds:indl),y(inds:indl));
end

if (indl>=inds)
%Do starting end:
ys = interp1(x,y,x1,'linear');
int = int + (x(inds)-x1)*(y(inds)+ys)/2;

%Do finishing end:
yl = interp1(x,y,x2,'linear');
int = int + (x2-x(indl))*(y(indl)+yl)/2;

else
ys = interp1(x,y,x1,'linear');
yl = interp1(x,y,x2,'linear');
int = (x2-x1)*(ys+yl)/2;
end
int = int*mult;
end