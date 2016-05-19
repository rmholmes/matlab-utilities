function r = mrdivide(x,y)
% Division x/y for AD objects
if isnumeric(y)
    tmp=y;
    y=[];
    y.tc=[tmp 0 0];
end
if isnumeric(x)
    tmp=x;
    x=[];
    x.tc=[tmp 0 0];
end
    xt1=x.tc(1);
    yt1=y.tc(1);
    r.tc(3) = (x.tc(3)*yt1^2 - 2*x.tc(2)*y.tc(2)*yt1 ...
                 + 2*xt1*y.tc(2)^2 - xt1*y.tc(3)*yt1)/yt1^3;
    r.tc(2) = (x.tc(2)*yt1 - xt1*y.tc(2))/yt1^2;
    r.tc(1) = xt1/yt1;
    r = class(r, 'AD');




