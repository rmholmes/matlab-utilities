function schrodObject = liouvilletransformation(slpObject)
%implementation of the Liouville transformation
%similar procedure as in Fortran code SLCPM12

p=slpObject.p;
w=slpObject.w;
equalpw = strcmp(func2str(w),func2str(p));
q=slpObject.q;

function x=fcn(r) 
    if equalpw
        x = ones(1,length(r));
    else
        x = sqrt(arrayfun(w,r)./arrayfun(p,r));
    end
end

function r = x2r(x)
%conversion from Schrod variable to SL variable

%preallocation:
r=zeros(1,length(x));

for j = 1:length(x)
  %determination of the interval
  tmp=find(xG>x(j));
  if ~isempty(tmp)
    i=tmp(1)-1;
  else
    i=length(xG)-1;  
  end
  
  r1 = rG(i);
  r2 = rG(i+1);
  g1 = xG(i)-x(j);
  gcc = fcn([r1 r2]);
  gc1 = gcc(1);
  g2 = xG(i+1)-x(j);
  gc2 = gcc(2);
  loop1 = true;
  its=0;
  while loop1
    its=its+1;
    gleft = g1;
    n = 0;
    while true
     n = n+1;
     rn = 0.1*n;
     rn2 = rn*rn;
     a1 = rn2*(3-2*rn);
     a0 = 1-a1;
     b0 = rn-rn2*(2-rn);
     b1 = rn2*(rn-1);
     gright = a0*g1+a1*g2+(b0*gc1+b1*gc2)*(r2-r1);      
     if (gleft*gright <= 0)
        break;
     else
        gleft=gright;
     end
    end
    %startpoint for Newton procedure
    rs = r1+(rn-0.1*gright/(gright-gleft))*(r2-r1);
    tol = 10*eps*max(abs(rs),1); 
    if (abs(r2-r1) < tol)
      r(j) = rs;
      break; 
    end 
    it=0;
    while true
    %Newton procedure
      it = it+1;
      if it >= 100
        r(j) = rs;
        loop1 = false;
        break;
      end
      qq = gauss(@fcn,rG(i),rs);
      gs = qq - x(j) + xG(i);
      gsp = fcn(rs);    
      drnewton = -gs/gsp;
      rnewton = rs+drnewton;
      if ((rnewton-r1)*(rnewton-r2)<=0)
        rs = rnewton;
        tol = 10*eps*max(abs(rs),1); 
        if (abs(drnewton) <= tol) 
            r(j) = rs;
            loop1 = false;
            break;
        end
      else
        if(gs*g1 < 0) 
          g1 = gs;
          r1 = rs;
          gc1 = gsp;
        else
          g2 = gs;
          r2 = rs;
          gc2 = gsp;
        end
      end
    end %while
 end %while
end %for
end


function pot=get_schrod_pot(x)
   r=x2r(x);
   po=zeros(1,length(r));
   dpo=zeros(1,length(r));
   ddpo=zeros(1,length(r));
   wo=zeros(1,length(r));
   dwo=zeros(1,length(r));
   ddwo=zeros(1,length(r));
   for i=1:length(r)
     rad=AD(r(i));  %automatic differentation
     tmpp=p(rad);
     po(i)=tmpp(1);
     tmpw=w(rad);
     wo(i)=tmpw(1);
     if isa(tmpp,'AD')
       dpo(i)=tmpp(2);  %als p='1' dan tmpp=1
       ddpo(i)=tmpp(3);
     end
     if isa(tmpw,'AD')
       dwo(i)=tmpw(2);   
       ddwo(i)=tmpw(3);
     end
   end   
   f1=dpo.*wo+po.*dwo;
   f2=ddpo.*wo+2*dpo.*dwo+po.*ddwo;
   f3=f1.*(dpo.*wo-po.*dwo);
   pot=((-3*f1.*f1+2*f3)./(po.*wo)+4*f2)./(16.*wo.*wo);
   pot=pot+q(r)./wo;
end


function x=m(r)
        x=(p(r).*w(r))^(-1/4);
end
function x=dm(r)
       rad=AD(r); tmpp=p(rad); tmpw=w(rad);
       if isa(tmpp,'AD')
         dp=tmpp(2);
       else
         dp=0;
       end
       if isa(tmpw,'AD')
         dw=tmpw(2);
       else
         dw=0;   
       end
       x=-0.25*m(r).^5*(dp.*tmpw(1)+tmpp(1).*dw); 
end

function x = r2x(r)
%conversion from SL variable to Schrod variable

lenr=length(r);
x = zeros(lenr,1);
for j = 1:lenr   
  %determination of the interval:
  if r(j)<rG(1)
      x(j)=0;
  else
      for i = 1:length(rG)-1
       if ((r(j)-rG(i))*(r(j)-rG(i+1)) <= 0) 
         break; 
       else
         x(j) = x(j) + Q(i);   
       end
      end  
      %only the integral from rG(i) up to r is evaluated
      x(j) = x(j) + gauss(@fcn,rG(i),r(j));     
  end
end
end

[mint,maxt,rG,xG,Q] = sl2sch(@fcn,slpObject.xmin,slpObject.xmax);
a0 = slpObject.a0*m(slpObject.xmin)^2+slpObject.b0*p(slpObject.xmin)*dm(slpObject.xmin)*m(slpObject.xmin);
b0 = slpObject.b0;
if isinf(a0) || (a0==0 && b0==0) || isnan(a0)
     a0=1;
     b0=0;
end
a1=slpObject.a1*m(slpObject.xmax)^2+slpObject.b1*p(slpObject.xmax)*dm(slpObject.xmax)*m(slpObject.xmax);
b1 = slpObject.b1;
if isinf(a1) || (a1==0 && b1==0) || isnan(a1)
     a1=1;
     b1=0;
end
schrodObject = slp('1',@get_schrod_pot,'1',mint,maxt,a0,b0,a1,b1);
schrodObject = addLiouvilleTransInfo(schrodObject,@x2r,@r2x,p,w);
end

function [xmin,xmax,rG,xG,q] = sl2sch(fcn,rmin,rmax)
% For given Sturm-Liouville problem, this returns
% the suitable x_min and x_max for the corresponding
% Schrodinger- problem. Also the rG - xG -points are
% returned.

i = 1;
xmin = 0;
rG(1) = rmin;
xG(1) = 0;
h = (rmax-rmin)/10;
hmin = 10^4*eps*abs(rmax-rmin); 
while true
  rG(i+1) = rG(i) + h;
  if (rG(i+1) > (rmax-hmin))
    rG(i+1) = rmax;
    h = rG(i+1) - rG(i);
  end
  hi = h;
  accurate = false;
  while ~(accurate) 
    mid = (rG(i)+rG(i+1))/2;
    q0 = gauss(fcn,rG(i),rG(i+1));
    q1 = gauss(fcn,rG(i),mid); 
    q2 = gauss(fcn,mid,rG(i+1));
    aq = abs(q1+q2);
    diffaq = abs(q1+q2-q0);
    tol=eps*max(aq,1); 
    if(diffaq < tol) || hi< hmin
       accurate = true;
    else
       hi = hi/2;
       rG(i+1) = rG(i) + hi;
    end
  end %while ~(accurate)
  q(i) = q1 + q2;
  xG(i+1) = xG(i) + q(i); 
  if(rG(i+1) >= rmax) 
    xmax = xG(i+1);
    return;
  end
  i = i+1;
end %while
end



function res = gauss(func,a,b)
%SLP\GAUSS res = gauss(func,a,b)
% calculates the integral of func by the 
% twelve-point Gauss quadrature formula.
x = [0.125233408511468915472441369464,...
   0.36783149899818019375269153664,...
   0.5873179542866174472967024189,...
   0.7699026741943046870368938332,...
   0.904117256370474856678465866,...
   0.9815606342467192506905490901];
w = [0.24914704581340278500056243604,...
   0.233492536538354808760849899,...
   0.2031674267230659217490645,...
   0.160078328543346226334653,...
   0.106939325995318430960255,...
   0.047175336386511827194616];

m = (b+a)/2;
d = (b-a)/2;
dx = d*x;
ff=func([m+dx m-dx]);
fps=ff(1:6);
fms=ff(7:end);
res = d*(w*(fps+fms)');
end


