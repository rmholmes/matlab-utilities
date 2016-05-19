function [e_loc,qq,ww,pp] = getErrorSLP(qfun,pfun,wfun,h,x)
%used in the automatic stepsize selection procedure, for problems in the
%full SLP form (not Liouville normal form)
[qq0,qq1,qq2]=getLegendrePolyn(qfun,x,h);
[ww0,ww1,ww2]=getLegendrePolyn(wfun,x,h);
[pp0,pp1,pp2]=getLegendrePolyn(@(x)1/pfun(x),x,h);
qq=[qq0 qq1 qq2];
ww=[ww0 ww1 ww2];
pp=[pp0 pp1 pp2];
qf=qfun(x+h); wf=wfun(x+h); pf=pfun(x+h);
h2=h^2;
if isnan(wf/pf)
    qe=abs(qf-(qq0+qq1/h2+qq2/h2));%/max(1,abs(qf)); 
    we=abs(wf-(ww0+ww1/h2+ww2/h2));%/max(1,abs(wf));
    pe=abs(1/pf-(pp0+pp1/h2+pp2/h2));%/max(1,abs(1/pf));
    e_loc=max([qe we pe]);
else
    [w0,w1,w2]=getLegendrePolyn(@(x)wfun(x)/pfun(x),x,h); 
    if qq0==0
      qe=abs(1/pf-(pp0+pp1/h2+pp2/h2)); 
    else
      [q0,q1,q2]=getLegendrePolyn(@(x)qfun(x)/pfun(x),x,h);
      qe=abs(qf/pf-(q0+q1/h2+q2/h2));%/max(1,abs(qf));
    end
    we=abs(wf/pf-(w0+w1/h2+w2/h2));%/max(1,abs(wf));
    e_loc=max([qe we]);
end
% if x==0
%     e_loc=e_loc*1000;
% end





function [V0,V1,V2]=getLegendrePolyn(fun,x,h)
%Gauss Legendre (3 Legendre points)
%realizes approximation by a series over shifted Legendre polynomials
% the F_s coefficients are returned
%cc=sqrt(15)/5;
%c=h/2*[-cc+1 1 cc+1];
%f=arrayfun(fun,x+c);
%V0=1/18*(8*f(2)+5*f(1)+5*f(3));
%V1=(3*(-1/18*15^(1/2)*(f(1)-f(3))))*h^2;
%V2=5/9*(-2*f(2)+f(1)+f(3))*h^2;
%4 nodes:
c1=sqrt((3-2*sqrt(6/5))/7);
c2=sqrt((3+2*sqrt(6/5))/7);
c=h/2*[-c2+1 -c1+1 c1+1 c2+1];
f=arrayfun(fun,x+c);
sq30=30^(1/2);
V0=1/360*sq30*((-5+3*sq30)*f(1)+(5+3*sq30)*f(2)+(5+3*sq30)*f(3)+(-5+3*sq30)*f(4));
V1=3/360*sq30*((-10*c(1)+5*h+6*c(1)*sq30-3*sq30*h)*f(1)+(10*c(2)-5*h+6*c(2)*sq30-3*sq30*h)*f(2)+(10*c(3)-5*h+6*c(3)*sq30-3*sq30*h)*f(3)+(-10*c(4)+5*h+6*c(4)*sq30-3*sq30*h)*f(4))*h;
V2=5*(-7*sq30*f(2)-7*sq30*f(3)+7*sq30*f(1)+7*f(4)*sq30)*h^2/360;
%V3=7*(1/12600*30^(1/2)*(15*(525-70*30^(1/2))^(1/2)+2*(525-70*30^(1/2))^(1/2)*30^(1/2))*f(2)+1/12600*30^(1/2)*(-15*(525-70*30^(1/2))^(1/2)-2*(525-70*30^(1/2))^(1/2)*30^(1/2))*f(3)+1/12600*30^(1/2)*(-15*(525+70*30^(1/2))^(1/2)+2*(525+70*30^(1/2))^(1/2)*30^(1/2))*f(1)+1/12600*30^(1/2)*(15*(525+70*30^(1/2))^(1/2)-2*(525+70*30^(1/2))^(1/2)*30^(1/2))*f(4));
