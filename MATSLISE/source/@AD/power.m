function r = power(x,y)
% x^y for AD objects.
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
L=log(x.tc(1));
yt1=y.tc(1);
yt2=y.tc(2);
xt1=x.tc(1);
r.tc(3) = (xt1^(yt1+2)*yt2^2*L^2 + ...
                  2*xt1^(yt1+1)*yt2*L*x.tc(2)*yt1 + ...
                  xt1^yt1*x.tc(2)^2*yt1^2 + ...
                  xt1^(yt1+2)*y.tc(3)*L + ...
                  2*xt1^(yt1+1)*x.tc(2)*yt2 + ...
                  xt1^(yt1+1)*x.tc(3)*yt1 - ...
                  xt1^yt1*yt1*x.tc(2)^2)/xt1^2;
    r.tc(2) = xt1^(yt1-1)*(yt2*L*xt1 ...
                 + x.tc(2)*yt1);
    r.tc(1) = xt1^yt1;
    r=class(r,'AD');
