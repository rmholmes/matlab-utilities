function meshData = computeMesh(slpObject,tol,varargin)
%constructs an appropriate mesh for the SLP-problem corresponding to
%slpObject. The size of the mesh intervals depends on a user input
%tolerance tol.
%For infinite problems, the mesh is constructed on a truncated interval,
%the size of this interval depends on the third input argument which is the
%index of the highest eigenvalue to be computed.

if nargin>2
    nmax=varargin{1};
end


%can half-range reduction be used?
if nargin>2 && testHalfRangeReduction(slpObject)
    %mesh needs to constructed only for the half-range
    %construct half-range problem:
    slp2=slp(slpObject.p,slpObject.q,slpObject.w,0,slpObject.xmax,0,1,slpObject.a1,slpObject.b1,slpObject.jumps(slpObject.jumps>0));
    %construct mesh for half-range problem:
    meshData=computeMesh(slp2,tol,nmax);
    meshData.halfRangeReduction=true;
    return;
else
    meshData.halfRangeReduction=false;
end

%initialisations
LNF=liouvilleNormalForm(slpObject); %is the problem in Liouville normal form?
meshData.LNF=LNF;
pot = slpObject.q; %the potential
wfunc = slpObject.w;
xmin = slpObject.xmin;   %mesh interval
xmax = slpObject.xmax;

hmin=1e-3;              %initial value for stepsize h of first interval around singularity
if xmax-xmin<1e-2
    hmin=(xmax-xmin)/10;
end

meshData.Infa=xmin;
meshData.Infb=xmax;

c=slpObject.classification;

%infinite endpoints
if isinf(xmin)||isinf(xmax)
    %find position of minimum in the potential
    xs=linspace(max(-20,xmin),min(20,xmax),400);
    [vmin,ix]=min(arrayfun(pot,xs));
    if pot(xs(end))==vmin && pot(xs(floor(length(xs)/2)))==vmin
        xvmin=xs(floor(length(xs)/2));
    else
        xvmin=xs(ix);
    end
end
if isinf(xmin)
    %compute an initial truncation point
    xmin=truncateInfiniteEndpoint(pot,wfunc,c,xmin,xmax,vmin,xvmin,max(nmax,7),-1);
    meshData.Infa=xmin;
end
if isinf(xmax)
    xmax=truncateInfiniteEndpoint(pot,wfunc,c,xmin,xmax,vmin,xvmin,max(nmax,7),1);
    meshData.Infb=xmax;
end


meshData.radial=false;
    
if LNF %problem in Liouville normal form
    %construct mesh for CPM{16,14} method
    %
    % this furnishes the vector of stepsizes (h) and the elements of the one
    % step propagator (matrices cu, cup, cv, cvp) + Vbar
    % and imatch (the interval at the rhs of which is the matching point)
    x = xmin;
    i=0;
    h = (xmax-xmin)/50;
    %Gauss-Legendre nodes and weights:
    nodes = [0.095012509837637440185, 0.281603550779258913230,...
       0.458016777657227386342, 0.617876244402643748447,...
       0.755404408355003033895, 0.865631202387831743880,...
       0.944575023073232576078, 0.989400934991649932596];
    weigths = [0.189450610455068496285, 0.182603415044923588867,...
       0.169156519395002538189, 0.149595988816576732081,...
       0.124628971255533872052, 0.095158511682492784810,...
       0.062253523938647892863, 0.027152459411754094852];
    P = 16;
    l1= legendre((1-nodes)/2,P);
    l2= legendre((1+nodes)/2,P);

    %trunca and truncb >0 for singular problems
    meshData.trunca=0; %number of intervals at each endpoint which need to be refined to check accuracy 
    meshData.truncb=0;

    if ~c.regular(1) && ~isinf(slpObject.xmin) % singular finite left endpoint
        if slpObject.xmin==0 && (c.regular(2) || isinf(slpObject.xmax)) 
            [x,l,Vs]=radialSchrodinger(pot,slpObject.xmax);
            %determines a narrow subinterval I1 around the origin, and
            %compute data for the algorithm consistent with the singular
            %nature around the origin
            if x>1e-12 && abs(imag(l))<1e-12 && abs(Vs(end))<1e4 %pryce8
                % potential does not fit the form l*(l+1)/x^2 +S(x)/x+R(x)
                meshData.r0=x;
                meshData.radial=true;
                meshData.l=l;
                meshData.Vs=Vs;
                meshData.Infa=x;
            else
                x=xmin;
            end
        end
        if ~meshData.radial
             %take as first meshinterval an interval of length hmin
             h=hmin;
             while pot(x+h)>1e10
                 h=h*2;
             end
             i=1;
             m = (x+x+h)/2; d = h/2;
             ff=arrayfun(pot,[m+d*nodes m-d*nodes]); %arrayfun: gives no problems when no vectorization is used in function definitions
             fpm=ff(1:8); fmm=ff(9:end);
             qpot = d*((weigths.*fmm)*l1+(weigths.*fpm)*l2);
             V0 = qpot(1)/h; %constant approximation of potential
             V = h*qpot(2:P+1).*(2*(1:P)+1);
             meshData.h(i)=h; meshData.V0(i)=V0;
             [cu(:,i),cup(:,i),cv(:,i),cvp(:,i)] = computeCoeffsBP(V); %CP coefficients
             [rcu(:,i),rcup(:,i),rcv(:,i),rcvp(:,i)] =computeCoeffsRP(V); 
             x = x + h;
             meshData.trunca=1;
        end
    end


    xmax2=xmax;
    if ~c.regular(2) && ~isinf(slpObject.xmax)
        xmax2=xmax-hmin;
    end

   %construct rest of the mesh (regular part)
    jumps=slpObject.jumps;
    if isempty(jumps)
        hmax=(xmax2-x)/2;
        while x < xmax2
            i=i+1;
            [h,V0,V]=getStep(x,h,min(hmax,xmax2-x),pot,l1,l2,tol);
            meshData.h(i)=h; meshData.V0(i)=V0;
            [cu(:,i),cup(:,i),cv(:,i),cvp(:,i)] = computeCoeffsBP(V); 
             %reference partition = basic partition but with a higher order
             %method
            [rcu(:,i),rcup(:,i),rcv(:,i),rcvp(:,i)] =computeCoeffsRP(V); 
            x = x + h;
        end
    else %problem with discontinuities
        %jumps need to be among the meshpoints
       for jump=[jumps xmax2]
           hmax=(jump-x);
           while x< jump
               i=i+1;
               [h,V0,V]=getStep(x,h,min(hmax,jump-x),pot,l1,l2,tol);
               meshData.h(i)=h; meshData.V0(i)=V0;
               [cu(:,i),cup(:,i),cv(:,i),cvp(:,i)] = computeCoeffsBP(V); 
               [rcu(:,i),rcup(:,i),rcv(:,i),rcvp(:,i)] =computeCoeffsRP(V); 
               x = x + h;
           end
       end
    end
    if ~c.regular(2) && ~isinf(slpObject.xmax)  %from xmax2 to xmax
        % last interval (for singular endpoint b)
         h=hmin;
         i=i+1;
         m = (x+x+h)/2; d = h/2;
         ff=arrayfun(pot,[m+d*nodes m-d*nodes]); %arrayfun: no problems due to forgotten vectorization in function definitions
         fpm=ff(1:8); fmm=ff(9:end);
         qpot = d*((weigths.*fmm)*l1+(weigths.*fpm)*l2);
         V0 = qpot(1)/h;
         V = h*qpot(2:P+1).*(2*(1:P)+1);
         meshData.h(i)=h; meshData.V0(i)=V0;
         [cu(:,i),cup(:,i),cv(:,i),cvp(:,i)] = computeCoeffsBP(V); 
         [rcu(:,i),rcup(:,i),rcv(:,i),rcvp(:,i)] =computeCoeffsRP(V); 
         meshData.truncb=1;
    end

    meshData.cu=cu';  %omdat kolommen toevoegen sneller gaat dan rijen toevoegen
    meshData.cup=cup'; meshData.cv=cv'; meshData.cvp=cvp';
    meshData.rcu=rcu';  meshData.rcup=rcup'; meshData.rcv=rcv'; meshData.rcvp=rcvp';


    %reference partition (used to form an estimation of the error in the
    %calculated eigenvalues): has equal stepsizes and meshpoints as the basic 
    %partition but a higher order method is used to advance the solution.
    
    %select an appropriate matching point
    [~,meshData.imatch] = min(meshData.V0(meshData.trunca+1:end-meshData.truncb));
    if meshData.trunca && meshData.truncb
        %avoid matchpoint near singular endpoints
        xs=xmin:(xmax-xmin)/100:xmax;
        m= median(arrayfun(pot,xs));
        ind=find(meshData.V0(meshData.trunca+1:end-meshData.truncb)>m);
        if ~isempty(ind)
          meshData.imatch=ind(1);
        end
    end

else %not in LNF
    % the sixth order CP algorithm for general Sturm-Liouville problems is used
    x=xmin;
    i=0;
    h=(xmax-xmin)/50;
    meshData.trunca=0;
    meshData.truncb=0;
   
    if ~c.regular(1) && ~isinf(slpObject.xmin)
         h=hmin;
         i=1;
         %[~,qq0,qq1,qq2,ww0,ww1,ww2,pp0,pp1,pp2] = getErrorSLP(slpObject.q,slpObject.p,slpObject.w,h,x);
         [~,qq,ww,pp] = getErrorSLP(slpObject.q,slpObject.p,slpObject.w,h,x);
         meshData.h(i)=h;
         meshData.Q0(i)=qq(1);meshData.Q1(i)=qq(2);meshData.Q2(i)=qq(3); %store data related to this mesh for later use
         meshData.W0(i)=ww(1);meshData.W1(i)=ww(2);meshData.W2(i)=ww(3);
         meshData.P0(i)=pp(1);meshData.P1(i)=pp(2);meshData.P2(i)=pp(3);
         x = x + h;
         meshData.trunca=1;
    end

    xmax2=xmax;
    if ~c.regular(2) && ~isinf(slpObject.xmax)
        xmax2=xmax-hmin;
    end
    
    jumps=slpObject.jumps;
    if isempty(jumps) 
        hmax=(xmax2-x)/10;
        while x<xmax2
               hmax=min(hmax,xmax2-x);
               hnew = 2*h;
               i=i+1;   
               it = 0;
               while abs(hnew/h-1) > 0.1 && it < 20  
                 h=max(min([hnew,hmax]),eps*10);      
                 it = it +1;         
                 [err,qq,ww,pp] = getErrorSLP(slpObject.q,slpObject.p,slpObject.w,h,x); %get_error geeft in feite V3/h^2
                 err=h^2*err;        
                 hnew = 5*h;
                 if(err > eps) 
                    hnew = h * (tol / err)^(1/6);
                 end
                 if abs(h -hmax)<eps && hnew>h, break; end
               end
               x=min(x+h,xmax2);
               if h>1e-14 || abs(x-xmax2)>eps
                   meshData.h(i)=h;
                   meshData.Q0(i)=qq(1);meshData.Q1(i)=qq(2);meshData.Q2(i)=qq(3); %store data related to this mesh for later use
                   meshData.W0(i)=ww(1);meshData.W1(i)=ww(2);meshData.W2(i)=ww(3);
                   meshData.P0(i)=pp(1);meshData.P1(i)=pp(2);meshData.P2(i)=pp(3);
               end
        end
    else
        for jump=[jumps xmax2];
             while x<jump
               hnew = 2*h;
               i=i+1;   
               it = 0;
               while abs(hnew/h-1) > 0.1 && it < 20
                 h=max(min(hnew,jump-x),eps*10);      
                 it = it +1;         
                 [err,qq,ww,pp] = getErrorSLP(slpObject.q,slpObject.p,slpObject.w,h,x); %get_error geeft in feite V3/h^2
                 err=h^2*err;        
                 hnew = 5*h;
                 if(err > eps) 
                    hnew = h * (tol / err)^(1/6);
                 end
                 if abs(h -(jump-x))<eps && hnew>h, break; end
               end
               x=x+h;
               meshData.h(i)=h;
               meshData.Q0(i)=qq(1);meshData.Q1(i)=qq(2);meshData.Q2(i)=qq(3); %store data related to this mesh for later use
               meshData.W0(i)=ww(1);meshData.W1(i)=ww(2);meshData.W2(i)=ww(3);
               meshData.P0(i)=pp(1);meshData.P1(i)=pp(2);meshData.P2(i)=pp(3);
            end
        end
    end
   
    
    if ~c.regular(2) && ~isinf(slpObject.xmax)
         h=hmin;
         i=i+1;
         [~,qq,ww,pp] = getErrorSLP(slpObject.q,slpObject.p,slpObject.w,h,x);
         meshData.h(i)=h; 
         meshData.Q0(i)=qq(1);meshData.Q1(i)=qq(2);meshData.Q2(i)=qq(3); %store data related to this mesh for later use
         meshData.W0(i)=ww(1);meshData.W1(i)=ww(2);meshData.W2(i)=ww(3);
         meshData.P0(i)=pp(1);meshData.P1(i)=pp(2);meshData.P2(i)=pp(3);
         meshData.truncb=1;
    end
end

function [r0,l,Vs]=radialSchrodinger(q,xmax)
 %find a good value for r0:
 %r0 should be small enough to ensure
 %both a strong domination of the reference potential (proportional to
 %r^(-2))
 %and a well fitting by parabolae of the
 %other potential terms
  if isinf(xmax)
      xmax=100;
  end
  r0=min(0.1,xmax/50);
  t1=2;t2=1;
  g=@(x)arrayfun(q,x).*x.^2;
  while abs(t1-t2)/abs(t1)>1e-6
      m = r0/2;
      x=[(1/3)*sqrt(5-2*sqrt(10/7)) (1/3)*sqrt(5+2*sqrt(10/7))]; %5 Legendre points in [0,r0]
      x=sort([m+m*x m m-m*x]);
      warning('off','all');
      p=polyfit(x,g(x),4); %gives coefficients of third degree polynomial through Legendre points
      warning('on','all');
      r1=r0/2;
      t1=q(r1); t2=polyval(p,r1)./r1.^2;
      r0=r1;
  end
r0=r0*2;
l=-1/2+sqrt(1+4*p(5))/2;
Vs=fliplr(p(1:4));


  

function [h,V0,V]=getStep(x,h,hmax,pot,l1,l2,tol)
%automatic stepsize selection algorithm.
%Gauss-Legendre nodes and weights:
nodes = [0.095012509837637440185, 0.281603550779258913230,...
   0.458016777657227386342, 0.617876244402643748447,...
   0.755404408355003033895, 0.865631202387831743880,...
   0.944575023073232576078, 0.989400934991649932596];
weigths = [0.189450610455068496285, 0.182603415044923588867,...
   0.169156519395002538189, 0.149595988816576732081,...
   0.124628971255533872052, 0.095158511682492784810,...
   0.062253523938647892863, 0.027152459411754094852];
it=0;
hnew=2*h;
P = 16;
while abs(hnew/h-1) > 0.1 && it < 20
    h=max(min([hnew,hmax]),eps*10);

    %[V0,V] = calculateV(pot,x,h);
    m = (x+x+h)/2; d = h/2;
    ff=arrayfun(pot,[m+d*nodes m-d*nodes]); %arrayfun: no problems due to forgotten vectorization in function definitions
    fpm=ff(1:8); fmm=ff(9:end);

    qpot = d*((weigths.*fmm)*l1+(weigths.*fpm)*l2);
    V0 = qpot(1)/h;
    V = h*qpot(2:P+1).*(2*(1:P)+1);

    err = estimateError(V);
    hnew = 5*h;
    if err > eps
         hnew = h * (tol / err)^(1/15);
    end
    if (abs(h -hmax)<eps) && hnew>h, break; end
    it=it+1;
end
if it==20 && hnew<h
    h=hnew;
    m = (x+x+h)/2; d = h/2;
    ff=arrayfun(pot,[m+d*nodes m-d*nodes]); %arrayfun: no problems due to forgotten vectorization in function definitions
    fpm=ff(1:8); fmm=ff(9:end);
    qpot = d*((weigths.*fmm)*l1+(weigths.*fpm)*l2);
    V0 = qpot(1)/h;
    V = h*qpot(2:P+1).*(2*(1:P)+1);
end

function hr=testHalfRangeReduction(slpObject)
hr=false;    
jumps=slpObject.jumps;
jumpsneg=jumps(jumps<0); jumpspos=jumps(jumps>0);
if any(abs(jumpsneg)~=jumpspos)
    return;
end
xmax=slpObject.xmax;
if isinf(slpObject.xmax)
    xmax=100;
end
xmin=slpObject.xmin;
if isinf(slpObject.xmin)
    xmin=-100;
end
if xmin==-xmax && slpObject.a0==slpObject.a1 && slpObject.b0==-slpObject.b1 && slpObject.q(xmin)==slpObject.q(xmax) && ...
       slpObject.p(xmin)==slpObject.p(xmax)  && slpObject.w(xmin)==slpObject.w(xmax)
    xmin=max(xmin,-10);
    xmax=min(xmax,10);
    xs=0:(xmax-xmin)/100:xmax;
    q1=arrayfun(slpObject.q,xs);
    q2=arrayfun(slpObject.q,-xs);
    if all(q1==q2)
         if ~liouvilleNormalForm(slpObject)
            q1=arrayfun(slpObject.p,xs);
            q2=arrayfun(slpObject.p,-xs);
            w1=arrayfun(slpObject.w,xs);
            w2=arrayfun(slpObject.w,-xs);
            if all(q1==q2) && all(w1==w2)
                hr=true;
            end
         else
             hr=true;
         end
    end
end




function x=truncateInfiniteEndpoint(pot,wfunc,c,xmin,xmax,vmin,xminx,wkb,sign)
%find an appropriate truncation point for the infinite endpoint
if any(c.cont_spectrum)
    x=xminx;
    v=pot(x)/wfunc(x);
    if v < c.cutoff 
        he=1e-5;
        vx=pot(x+sign*he)/wfunc(x+sign*he);
        vx2=pot(x)/wfunc(x);
        while (vx-vx2)/he>1e2 && c.cutoff-vx>10
             he=he*1.2;
             x=x+sign*he;
             vx2=vx;
             vx=pot(x)/wfunc(x);
        end
    end
    if v>c.cutoff && (c.cont_spectrum(2) && sign==1)
        xs=x:sign*x:sign*10*x;
        ind=find(arrayfun(pot,xs)./arrayfun(pot,xs)<c.cutoff);
        if ~isempty(ind)
            x=xs(ind);
        end
    else
        sums = 0;
        vmin =c.cutoff;
        he=0.01;
        vx_prev=vmin;
        vx=v;
        while sums<wkb && abs(vx-c.cutoff)>1e-10
           x=x+sign*he;
           vx=pot(x)/wfunc(x);
           if vx > vmin && vx>vx_prev 
              sums=sums+sqrt(vx-vmin)*he;
           else
             vmin=min(vx,c.cutoff);
             sums=0;
           end
           vx_prev=vx;
           he=min(he*1.1,(abs(x)+1)/10);
        end
    end
else
     if sign>0
         x=max([10,xminx,xmin+2]);
     else
         x=min([-10,xminx,xmax-2]);
     end
     vx=pot(x)/wfunc(x);
     vx2=pot(x/2)/wfunc(x/2);
     while vx>vx2 && sqrt(vx-vx2)*abs(x/2)>wkb*10  %ook voor min inf
      x=x/2;
      vx=pot(x)/wfunc(x);
      vx2=pot(x/2)/wfunc(x/2);
     end
     x=x*1.2;
 end    


