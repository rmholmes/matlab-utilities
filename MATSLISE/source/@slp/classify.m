function [c,slp_struct] = classify(slp_struct)
%CLASSIFY returns classification information for a Sturm-Liouville problem


% SLEDGE classification algorithm
% by C.T. Fulton, S. Pruess and Y. Xie
% classifies the SLP into one of ten spectral categories based on the
% endpoint classification

%  Note:
%     Any computational algorithm must be based on a finite
%     amount of information; hence, there will always be cases that
%     any algorithm misclassifies.  In addition, some problems are
%     inherently ill-conditioned, in that a small change in the 
%     coefficients can produce a totally different classification.

a=slp_struct.xmin;
b=slp_struct.xmax;
Q=slp_struct.q;
W=slp_struct.w;
P=slp_struct.p;
bcs = [slp_struct.a0 slp_struct.b0 slp_struct.a1 slp_struct.b1];
c.msg='';

kmax = 16; %number of points sampled for singular problems.
%  By increasing kmax, the reliability of the classification may increase;
%  however, the computing time may also increase.

if isinf(a)
      z(:,1) = -8.^(1:kmax);    %sample points
else
      z(:,1) =a+(1/8).^(1:kmax);    %sample points
end
if isinf(b)
      z(:,2) = 8.^(1:kmax);    %sample points
else
      z(:,2) =b-(1/8).^(1:kmax);    
end
kused(1) = kmax; kused(2) = kmax;
LNF=true;                       %is the problem in Liouville normal form?
%arrayfun: dan geen problemen als bvb p=1
QZ(:,1)=arrayfun(Q,z(:,1));QZ(:,2)=arrayfun(Q,z(:,2));
PZ(:,1)=arrayfun(P,z(:,1));PZ(:,2)=arrayfun(P,z(:,2));
WZ(:,1)=arrayfun(W,z(:,1));WZ(:,2)=arrayfun(W,z(:,2));
if any(PZ(:,1)<=0) || any(WZ(:,1)<=0) || any(WZ(:,2)<=0) || any(PZ(:,2)<=0)
    error('The functions p(x) or w(x) are not positive');
end
if any(PZ(:,1)~=PZ(1,1)) || any(WZ(:,1)~=WZ(1,1)) || any(PZ(:,2)~=PZ(1,1)) || any(WZ(:,2)~=WZ(1,1)) || any(PZ(:,1)~=WZ(1,1))
    LNF=false;
end
c.LNF=LNF;
for k=1:kmax 
    if abs(log(PZ(k,1)))>300 || isnan(PZ(k,1)) || isnan(QZ(k,1))|| isnan(WZ(k,1)) || ...
            abs(log(1+abs(QZ(k,1))))>300 || abs(log(WZ(k,1)))>300 ||  abs(a-z(k,1))<1e-14 
         kused(1)=k-1;
         break;
     end
end
for k=1:kmax 
    if abs(log(PZ(k,2)))>300 || isnan(PZ(k,2)) || isnan(QZ(k,2))|| isnan(WZ(k,2)) || ...
            abs(log(1+QZ(k,2)))>300 || abs(log(WZ(k,2)))>300 ||  abs(b-z(k,2))<1e-14    
         kused(2)=k-1;
         break;
     end
end
c.bcs = bcs;
[cev1,cont_spectrum1,regular1,limit_circle1,oscillatory1,BC,uncertain1,enda,msg]=classify_endpoint(z(1:kused(1),1),...
    PZ(1:kused(1),1),QZ(1:kused(1),1),WZ(1:kused(1),1),P,Q,W,a,1,LNF);
c.msg=[c.msg msg];
if (cont_spectrum1 || ~oscillatory1 ) && ~regular1
    c.bcs(1)=BC(1);
    c.bcs(2)=-BC(2);
end
[cev2,cont_spectrum2,regular2,limit_circle2,oscillatory2,BC,uncertain2,endb,msg]=classify_endpoint(z(1:kused(2),2),...
    PZ(1:kused(2),2),QZ(1:kused(2),2),WZ(1:kused(2),2),P,Q,W,b,-1,LNF);
c.msg=[c.msg msg];
if (cont_spectrum2 || ~oscillatory2 ) && ~regular2
    c.bcs(3)=BC(1);
    c.bcs(4)=-BC(2);
end
c.cutoff = min(cev1,cev2);
c.end=[enda endb];

if uncertain1 || uncertain2
    c.msg=[c.msg 'Classification is not certain. '];
end

%find the number of eigenvalues below the start of the continuous spectrum
lastev=inf; %lastev=number of eigenvalues below the continuous spectrum
if (~oscillatory1 && ~limit_circle2 && oscillatory2) || ...
    (~oscillatory2 && ~limit_circle1 && oscillatory1)
    lastev=0; %There are no eigenvalues below the start of the continuous spectrum.
elseif (oscillatory1 && ~limit_circle1 && oscillatory2 && limit_circle2) || ...
    (oscillatory2 && ~limit_circle2 && oscillatory1 && limit_circle1)
    lastev=0;
elseif (cont_spectrum1 && oscillatory2) || (cont_spectrum2 && oscillatory1)
    lastev=0;
elseif (cont_spectrum1 && ~oscillatory2) ||(cont_spectrum2 && ~oscillatory1) || ...
        cont_spectrum1 && cont_spectrum2
     end1=enda;
     end2=endb;
      %if singular point, find truncated point:
     if ~cont_spectrum1 && ~isinf(enda) && (isinf(Q(end1))||isnan(Q(end1)) || (Q(end1)-c.cutoff)>1e-10)
       if enda<0  
           tmp=enda:0.001:min(0,end2);
           ind=find((arrayfun(Q,tmp)-c.cutoff)<1e7);
           if ~isempty(ind)
               end1=tmp(ind(1));
           end
       else
           end1=end1+0.001;
       end
     end
     if ~cont_spectrum2 && ~isinf(endb) && (isinf(Q(end2))||isnan(Q(end2)) || (Q(end2)-c.cutoff)>1e-10)
       if endb>0
           tmp=endb:-0.001:max(0,end1);
           ind=find(arrayfun(Q,tmp)-c.cutoff<1e7);
           if ~isempty(ind)
                end2=tmp(ind(1));
           end
       else
           end2=end2-0.001;
       end
     end
     if c.LNF && cont_spectrum2 %our interval must be large enough
       v2=Q(end2);
       if isinf(b) && v2 < c.cutoff 
         while c.cutoff-v2 > 0.1
           end2 = end2*2;
           v2 = Q(end2);
         end
       end
     end
     if ~c.LNF && cont_spectrum2
         end2=end2*2;
     end
     if c.LNF && cont_spectrum1 %our interval must be large enough
       v1=Q(end1);
       if isinf(a) && v1 < c.cutoff 
         while c.cutoff-v1 > 0.1
           end1 = end1*2;
           v1 = Q(end1);
         end
       end
     end
     if ~c.LNF && cont_spectrum1
         end1=end1*2;
     end
     s=slp(P,Q,W,end1,end2,bcs(1),bcs(2),bcs(3),bcs(4));
     meshData = computeMesh(s,1e-5);
     end1s=end1;
     end2s=end2;
     n1n=locateEigenvalue(s,c.cutoff,meshData,false);
     n1=floor(n1n);
     if c.LNF && cont_spectrum2 %our interval must be large enough
       end2=end2*2;
       v2=Q(end2);
       if isinf(b) && v2 < c.cutoff 
         while c.cutoff-v2 > 1e-3
           end2 = end2*2;
           v2 = Q(end2);
         end
       end
     end
     if c.LNF && cont_spectrum1 %our interval must be large enough
       end1=end1*2;
       v1=Q(end1);
       if isinf(a) && v1 < c.cutoff 
         while c.cutoff-v1 > 1e-3
           end1 = end1*2;
           v1 = Q(end1);
         end
       end
     end
     if ~c.LNF && cont_spectrum2
         end2=end2*2;
     end
     if ~c.LNF && cont_spectrum1
         end1=end1*2;
     end
     if end1s~=end1 || end2s~=end2
         s=slp(P,Q,W,end1,end2,bcs(1),bcs(2),bcs(3),bcs(4));
     end
     meshData = computeMesh(s,1e-8);
     if ~regular1 && ~cont_spectrum1 && enda==0 && c.LNF
         [r0,l,Vs]=radialSchrodinger(Q,end1);
         if abs(imag(l))<1e-12
             meshData.l=l; meshData.r0=r0; meshData.Vs=Vs;
             meshData.radial=1;
         end
     end
     n2n=locateEigenvalue(s,c.cutoff,meshData,false);
     n2=floor(n2n);
     lastev=n2;
     end1s=end1;
     end2s=end2;
     if cont_spectrum2 && isinf(b)
           end2 = end2*2;
     end
     if  cont_spectrum1 && isinf(a)
         end1 = end1*2;
     end
     if end1s~=end1 || end2s~=end2
         s=slp(P,Q,W,end1,end2,bcs(1),bcs(2),bcs(3),bcs(4));
         meshData = computeMesh(s,1e-8); 
         if ~regular1 && ~cont_spectrum1 && enda==0 && c.LNF
           [meshData.r0,meshData.l,meshData.Vs]=radialSchrodinger(Q,end1);
           meshData.radial=1;
         end
         n3n=locateEigenvalue(s,c.cutoff,meshData,false);
         n3=floor(n3n);
         if n3>n2
             lastev=inf;
         end
     end
     if n1 >n2 && length(c.msg)<1
           c.msg='The eigenvalue count is uncertain';
     end
end
c.lastev  = lastev;
c.regular=[regular1 regular2];
c.limit_circle=[limit_circle1 limit_circle2];
c.cont_spectrum=[cont_spectrum1 cont_spectrum2];
c.cev=[cev1 cev2];
if cont_spectrum1
    oscillatory1=false;
end
if cont_spectrum2
    oscillatory2=false;
end
c.oscillatory=[oscillatory1 oscillatory2];
%classify the problem in one of the 10 categories
if ~oscillatory1 && ~oscillatory2 && ~cont_spectrum2 && ~cont_spectrum1
    category=1;
elseif (~oscillatory1 && ~cont_spectrum1 && cont_spectrum2) || (~oscillatory2 && ~cont_spectrum2 && cont_spectrum1)
    category=2;
elseif (~oscillatory1 && ~cont_spectrum1 && oscillatory2 && limit_circle2 && ~cont_spectrum2) ||...
        (~oscillatory2 && ~cont_spectrum2 && oscillatory1 && limit_circle1 && ~cont_spectrum1)
    category=3;
elseif (~oscillatory1 && ~cont_spectrum1 && ~limit_circle2 && oscillatory2 && ~cont_spectrum2) ||...
        (~oscillatory2 && ~cont_spectrum2 && ~limit_circle1 && oscillatory1 && ~cont_spectrum1)
    category=4;
elseif limit_circle1 && oscillatory1 && limit_circle2 && oscillatory2
    category=5;
elseif oscillatory1 && ~limit_circle1 && ~cont_spectrum1 && oscillatory2 && ~limit_circle2 && ~cont_spectrum2
    category=6;
elseif oscillatory1 && oscillatory2 && ~(cont_spectrum1 || cont_spectrum2) && ((limit_circle1 && ~limit_circle2) ||(limit_circle2 && ~limit_circle1))
    category=7;
elseif (limit_circle1 && oscillatory1 && cont_spectrum2)||(limit_circle2 && oscillatory2 && cont_spectrum1)
    category=8;
elseif (~limit_circle1 && oscillatory1 && ~cont_spectrum1 && cont_spectrum2) || (~limit_circle2 && oscillatory2 && ~cont_spectrum2 && cont_spectrum1)
    category=9;
elseif cont_spectrum1 && cont_spectrum2
    category=10;
end
c.category=category;


function [r0,l,Vs]=radialSchrodinger(q,r0)
 g=@(x)arrayfun(q,x).*x.^2;
 m = r0/2;
 x=[(1/3)*sqrt(5-2*sqrt(10/7)) (1/3)*sqrt(5+2*sqrt(10/7))]; %5 Legendre points in [0,r0]
 x=sort([m+m*x m m-m*x]);
 warning('off','all');
 p=polyfit(x,g(x),4); %gives coefficients of third degree polynomial through Legendre points
 warning('on','all');
 l=-1/2+sqrt(1+4*p(5))/2;
 Vs=fliplr(p(1:4));

%-----------------------------------------------------------------
function [cev,cont_spectrum,regular,limit_circle,oscillatory,BC,uncertain,ende,msg]=classify_endpoint(z,pz,qz,wz,P,Q,W,endpoint,direction,LNF)
%information about the problem at the singular endpoint (called "endpoint")
%is passed through the variables cev,cont_spectrum,regular,...

%initializations:
U=realmin*10; %small value
cont_spectrum=false;
irreg=false;
oscillatory=false;
limit_circle=false;
KCLASS = 0;
cev=1/U;
BC=[0;0];
ende=endpoint;
msg='';

%seek monomial approximation 
[eq,cq,qosc,~,uncertain]=powerf(z,qz,endpoint,direction);
if abs(cq) <= 1e-14
    cq=0;
    eq=0;
end
if LNF 
   ep = 0;
   cp = pz(1);
   ew = 0;
   cw = wz(1);
   posc = false;
   wosc = false;
else
    [ep,cp,posc,~,uncertain]=powerf(z,pz,endpoint,direction);
    if abs(cp) <= 1e-14 
        ep = 0;
    end
    [ew,cw,wosc,~,uncertain]=powerf(z,wz,endpoint,direction);
    if abs(cw) <= 1e-14 
        ew = 0;
    end
end

if posc || wosc
    if isinf(endpoint)
         msg='p(x) or w(x) is not well approximated by a power potential, classification is uncertain!. ';
         KCLASS = 1;
    end
    limit_circle= true;
    oscillatory = false;
end


if qosc
  msg=' q(x) is not well approximated by a power potential, classification is uncertain!' ;
  uncertain=true;
  if isinf(endpoint)
      %oscillatory coefficient function
      KCLASS=1;
      limit_circle=false;
      cont_spectrum=true;
      oscillatory = false;
      BC(1)=1;
      BC(2)=0;
      kmax = length(qz);
      cev = qz(kmax-1);
      delta = (z(kmax)-z(kmax-1))/(41);
      for i = 0:40
         zz = z(kmax) - i*delta;
         cp=P(zz);
         cq=Q(zz);
         cw=W(zz);
         cev = min(cev,cq);
      end
      if abs(cev) < 1e-14 
          cev = 0;
      end
      if U*abs(cev) >= 1 
           cont_spectrum = false;
      end
  else
      limit_circle=true;
      oscillatory = false;
  end
  eqlnf=0;
end

%analyze this endpoint
if ep<1 && eq > -1 && ew>-1 && ~isinf(endpoint) 
    regular=true;
    limit_circle=true;
    oscillatory=false;
    cont_spectrum=false;
    if ep>0  || eq < 0 || ew < 0
        KCLASS = 2;
    end
   return
end

regular=false;
eta(1) = (ew-ep+2)/2;
if abs(eta(1))<=1e-14 
    eta(1) = 0;
end
if eta(1)~=0 
  eqlnf = (eq-ew)/eta(1);
  if eqlnf>=0
    iqlnf=eqlnf+0.5*eqlnf;
  else
    iqlnf=eqlnf-0.5*eqlnf;
  end
  if abs(iqlnf-eqlnf)< 1e-13
      eqlnf = iqlnf;
  end
  c1 = (cq/cw)*(abs(eta(1))*sqrt(cp/cw))^eqlnf;
  if c1 == 0 
      eqlnf = 0;
  end
  c2 = (ep+ew)/(4*eta(1));
  c2 = c2* (c2-1);
else
  c3 = (cp/cw)* ((ep+ew)/4)^2;
end

if isinf(endpoint)
   sgn=-1; 
   %make an initial estimate for "infinity" 
    if eq>eq || cq==0 
        if ew~=ep
           gamma = ew - ep;
           delta = cw/cp;
        else
           gamma = eq - ep;
           delta = abs(cq)/cp;
           if delta ==0 
               delta = 1;
           end
        end
     else
        if eq>ew
            gamma = eq - ep;
            delta = abs(cq)/cp;
        else
            if ew>ep
                 gamma = ew - ep;
                 delta = cw/cp;
            else
                 gamma=0;
                 if ew== ep
                     delta = abs(cw-cq)/cp;
                 else
                     delta = abs(cq)/cp;
                 end
            end
        end
    end

    if gamma >0.5
       if gamma < 2 %potentiaal stijgt geleidelijk
                ende=80;
       else
         ende=min(max(64/((2*gamma-3)*delta^(1/(gamma+2))),1),80);
         if gamma > 24  %sterke stijging in potentiaal
             ende = 12;
         end
       end
    else
       if gamma <-0.5
          ende=max(min(600*delta^(1/gamma)*5^gamma,120),2);
       else
          ende = 12;
       end 
       if gamma == 0 && cq ~=0
          ende =40;
       end
    end
   if endpoint<0 && ende>0
      ende = -ende;
   end
else
    %test for finite irregular singular points
    sgn = 1;
    if ew-ep>=0
        g = ew - ep + 0.5;
    else
        g = ew -ep - 0.5;
    end
    if eq-ep>=0
        k=eq-ep+0.5;
    else
        k=eq-ep-0.5;
    end
    irreg = true;
    if cq == 0
         if g >= -2 && abs(ew-ep-g) <=1e-14
             irreg = false;
         end
    else
        if ew<=eq && g >= -2 && abs(ew-ep-g)<= 1e-14
            irreg = false;
        end
        if ew>eq && k>=-2 && abs(eq-ep-k) <= 1e-14
            irreg = false;
        end
    end
end

if sgn*eta(1) > 0 
    %Carry out the Case 1 tests.
    osc = 0;
    k = 0;
    if eqlnf < -2
       if cq < 0
           k = 1;
       end
       if cq > 0
           k = -1;
       end
    end
    if eqlnf == -2
       if abs(c1+c2+0.25)<= 1e-13 
             msg='  WARNING: borderline nonoscillatory/oscillatory classification. ';
             k = -1;
       else
         if c1+c2 < -0.25-1e-13 
             k = 1;
         end
         if c1+c2>-0.25
             k = -1;
         end
       end
    end
    if eqlnf > -2 
        if abs(c2+0.25) <= 1e-13
            c2 = -0.25;
        end
        if c2>=-0.25
            k = -1;
        end
    end 

   if k == 1 
      oscillatory = true;
   else
     if k == -1
         oscillatory = false;
     else
         msg='  NO INFORMATION on osc/nonosc class.';
     end
   end

   k=0;
   if eqlnf<-2 
    if cq>0
        k = -1;
    end
    if cq < 0
        k = 1;
    end
   end
   if eqlnf==-2 || abs(eqlnf+2)<1e-14 
     if abs(c1+c2-3/4)<=1e-13
           k = -1;
           msg='  WARNING: borderline Lc/Lp classification. ';
           uncertain=true;
     elseif (c1+c2>=3/4) 
         k = -1;
     elseif abs(c1+c2)<3/4-1e-13
         k = 1;
     elseif c1+c2 <-1e-13 
         k = 1;
     end 
   end
  if eqlnf > -2
     if abs(c2-3/4)<=1e-13
          k = -1;
           msg='  WARNING: borderline Lc/Lp classification. ';
           uncertain=true;
     end
     if c2>=3/4
         k = -1;
     end
     if abs(c2)<3/4-1e-13
         k = 1;
     end
     if c2<-1e-13
         k = 1;
     end
  end
  if (k == 1) 
            limit_circle = true;
  elseif k == -1
             limit_circle = false;
  else
             limit_circle = -1;
             msg='  NO INFORMATION on Lp/Lc class. ';
   end
end

%%%%%%%%%%%%%


     if sgn*eta(1)< 0 
        %Carry out the Case 2 tests
         osc = 0;
         if ((eqlnf > 0) && (cq < 0)) 
             osc = 1; %OSC
         elseif ((eqlnf > 0) && (cq > 0)) 
             osc = -1; %NONOSC
         end
         if (eqlnf == 0) %NONOSC/OSC with cutoff value 
            osc = -1;
            cev = cq; %then cutoff value= cq (Theorem 4)
            cont_spectrum = true;
            if (U*abs(cev) >= 1) 
                cont_spectrum = false;
            end
         end
         if (eqlnf < 0) %NONOSC/OSC with cutoff value
            osc = -1;
            cev = 0; %cutoff
            cont_spectrum = true;
         end
         if (osc == 1) 
            oscillatory = true;
         else
            if (osc == -1) 
               oscillatory = false;
            else
               msg='  NO INFORMATION on Osc/Nonosc class. ';
            end
         end
         lc = 0;
         if ((eqlnf > 2) && (cq > 0)) || eqlnf <=2  || cq ==0
             lc = -1; %limit point
         end
         if ((eqlnf > 2) && (cq < 0)) 
             lc = 1; %limit circle
         end
         if (lc == 1) 
            limit_circle = true;
         else
            if (lc == -1) 
               limit_circle = false;
            else
               limit_circle = -1; 
               msg=' NO INFORMATION on Lp/Lc class. ';
            end
         end
     end
   
    if eta(1)==0  %Carry out the Case 3 and 4 tests.
        if (sgn* (eq-ew)<0 && cq < 0)
            oscillatory = true;
            limit_circle = true;
        elseif (sgn* (eq-ew)<0 && cq>0)
            oscillatory = false;
            limit_circle = false;
        elseif eq==ew 
            oscillatory = false;
            limit_circle = false;
            cev = cq/cr + c3;
            cont_spectrum  = true;
            if U*abs(cev)>=1 
                cont_spectrum = false;
            end
        elseif (sgn* (eq-ew)>0 || cq==0)
            oscillatory = false;
            limit_circle = false;
            cev = c3;
            cont_spectrum = true;
             if U*abs(cev)>=1 
                cont_spectrum = false;
            end
        end
    end



      if (abs(cev) <= 1e-13) 
          cev = 0;
      end

%     Calculate the Friedrichs boundary condition (if appropriate).

      if cont_spectrum 
         BC(1) = 1;
         BC(2) = 0;
      end
      if (~cont_spectrum && ~oscillatory ) 
         if (sgn* (ew+ep)>0 && sgn* (eq+ep)>0)
            BC(1) = 0;
            BC(2) = 1;
         else
            if (sgn* (ew+ep)>0 && (eq+ep) == 0)
               BC(1) = sqrt(cp*abs(cq));
               BC(2) = 1;
               if (BC(1)>1)
                  BC(2) = 1/BC(1);
                  BC(1) = 1;
               end
            else
               if (sgn* (ew+ep)<0 || sgn*(eq+ep)< 0)
                  BC(1) = 1;
                  BC(2) = 0;
               end
            end
         end
      end

      if (~oscillatory) 
         if isinf(endpoint) && (eqlnf == -1) 
           %disp('infinite endpoint, eqlnf=-1')
           KCLASS=3;
         end
         i = ew - ep;
         if cq ~= 0
            k = eq - ep;
            if cq> 0
               if (k <i)
                   i = 0;
               end
            else
               i = min(i,k);
            end
         end

         if (~isinf(endpoint) && (i < 0)) 
%              Transform some nonoscillatory problems for which Tau
%              is unbounded near a finite singular endpoint.
               if irreg
                   KCLASS = 6;
               end
               if eta(1)>0 
                if irreg
               %    disp('"hard" irregular, eta(1)>0');
                   KCLASS = 10;
                end
%               Transform some nonoscillatory problems for which Tau
%               is unbounded near a finite singular endpoint.
                EMU = 0;
                D(1) = (eta(1)*sqrt(cp/cw))^(1/eta(1));
                D(2) = 1/sqrt(sqrt(cp*cw*D(1)^(ep+ew)));
                if (eqlnf==-2) || (c1==0)
                  if ~irreg
                      KCLASS = 9;
                  end
                  EMU = abs(0.25+c1+c2);
                  if EMU < 1e-13
                     EMU = 0.5;
                  else
                     EMU = 0.5 + sqrt(EMU);
                  end
                else
                  if ~irreg
                      KCLASS = 8;
                  end
                end
                eta(2) = EMU - 0.25* (ep+ew)/eta(1);
                if ((KCLASS==10) && (EMU == 0))
                  eta(2) = 0.5*(1-ep)/eta(1);
                  EMU = eta(2) + 0.25*(ep+ew)/eta(1);
                end
                D(3) = eta(2)* (eta(2)+ (ep-1)/eta(1));
                D(4) = D(3);
                if eqlnf == -2
                   D(4) = D(4) - c1;
                end
                if (abs(D(4))<=1e-13)
                   D(4) = 0;
                end
                D(4) = sqrt(abs(D(4)));
              end %eta(1)>0!!!!!
              if (KCLASS>=9)
                %warning('  This problem has unbounded [Ev*r(x)-q(x)]/p(x).')
                   if ((EMU>0) ||(ep+ew~=0)) 
                       %disp('A change of variables will be used near this
                       %endpoint.')
                   end
              end
            end

      if (~isinf(endpoint) && (~regular) && KCLASS == 0)
          KCLASS = 4;
      end
      if (qosc||posc||wosc) && isinf(endpoint) 
          ende = 99;
          if endpoint<0
              ende=-99;
          end
      end
      %disp(['  Classification type (KCLASS) is: ',num2str(KCLASS)]);
  end
   
     
%-------------------------------------------------------------------
function [EF,CF,OSC,exact,uncertain] = powerf(x,f,e,direction)
% Find the power function which "dominates" the tabled
% coefficient function.  The output is Cf and Ef such that
%           f(x)  is asymptotic to  Cf*x^Ef or Cf*(x-e)^Ef.
%     The vectors X(*) and F(*) hold the N input points:
%           F(I) = f(X(I)) I = 1,...,N.  
%     e=the value of the endpoint


uncertain=false;
OSC=false;

 %Estimate the exponent.
n = length(f);
if isinf(e)
  for k = 1:n-1
    if (f(k) ~= 0) && (f(k+1) ~= 0) 
         y(k) = log(abs(f(k+1)/f(k)))/log(abs(x(k+1)/x(k)));
    else
         y(k) = 0;
    end
  end   
 else
  for k = 1:n-1
    if (f(k) ~= 0) && (f(k+1) ~= 0)  
         y(k) = log(abs(f(k+1)/f(k)))/log(abs((x(k+1)-e)/(x(k)-e)));
    else
         y(k) = 0;
    end
  end
end
[EF,err] = aitken(n-1,y); % Use Aitken's algorithm to accelerate convergence of the sequence in y
k = EF+sign(EF)*0.5;
if EF == 0
    k = EF+0.5;
end
if abs(k-EF) <= 1e-6
    EF = k;
end
if abs(err) > 1e-13*max(1,abs(EF)) 
%       There is uncertainty in the exponent.
        uncertain = true;
        %disp('classification is not certain');
end
if (abs(EF) <= 10*eps) 
    EF = 0;
end

% Estimate the coefficient.
if isinf(e)
  for k = 1:n-1
   y(k) = f(k)/abs(x(k))^EF; 
  end  
else
  for k = 1:n-1
   y(k) = f(k)/abs(x(k)-e)^EF;
  end
end
[CF,err] = aitken(n-1,y);
if isinf(e)
    direction=direction*(-1);
    e=0;
end
% CF=CF
% EF=EF
% abs(f(n)-CF*(direction*(x(n)-e))^EF) 
% 20.0*1e-11*abs(f(n)) 
if (EF > 20) && (abs(CF) <= 1e-13) 
%   Coefficient probably has exponential behaviour.
       if CF >= 0
         CF = 1;
       else
         CF = -1;
       end
elseif ((abs(err) > 1e-13*max(1,abs(CF))) || ...
       ((abs(f(n)-CF*(direction*(x(n)-e))^EF) > 20.0*1e-11*abs(f(n))) && (EF ~= 0 ))) 
%      There is uncertainty in the coefficient; call such
%      cases oscillatory.
       uncertain=true;
       OSC = true;
end
if (abs(CF) > 1e7) 
    exact = false;
    return
end
k = CF+0.5;
if abs(k-CF) <= 1e-7  && k ~= 0 
    CF = k;
end
exact = true;
for k = 1:n
    if (abs(f(k)-CF*(direction*(x(k)-e))^EF) > 1e-11*abs(f(k))) 
        exact = false;
    end
end


%-------------------------------------------------------------
function [xlim,err] = aitken(n,x)
%     Use Aitken's algorithm to accelerate convergence of the sequence
%     in x(*).

if n <= 2 
   xlim = x(n);
   err = 0;
   return;
end
xold = 1e30;
for i = 1:n-2
   denom = x(i+2)-2*x(i+1)+x(i);
   if denom ~= 0 
        xlim = x(i)-(x(i+1)-x(i))^2/denom;
        err = xlim-xold;
   else
        err = x(i+2)-x(i+1);
        xlim = x(i+2);
   end
   if (abs(err) < max(1,abs(xlim))*1e-13) 
       return
   end
   xold = xlim;
end




function [theta0,theta1]=prufer(y0,y0p,y1,y1p,e,v,h,S)
if e > v %(i): e > Vi
    Si = sqrt(e-v); %local scale factor
    if y0p ~= 0
       theta0 = atan(Si*y0/y0p);
    elseif y0 > 0
         theta0 = pi/2;
    elseif y0 < 0
         theta0 = -pi/2;
    else
        theta0=0;
    end
    nphi = floor((theta0 + Si*h)/pi); 
    phibar = theta0 + Si*h - nphi*pi; %between -pi/2 and pi/2
    if y1p ~= 0
        phistar = atan(Si*y1/y1p);
    elseif y1 > 0
        phistar = pi/2;
    elseif y1 < 0
        phistar = -pi/2;
    else
        phistar = 0;
    end
    deltaphi = phistar - phibar;
    if deltaphi > pi/2
        deltaphi = deltaphi-pi; 
    elseif deltaphi < -pi/2
        deltaphi = deltaphi + pi;
    end
    theta1 = theta0 + Si*h + deltaphi;
    if Si == S
        return;
    end
%   else rescaling procedure
    sigma = S/Si;
    theta0 = theta0 + atan2((sigma-1)*sin(theta0)*cos(theta0),...
        1+(sigma-1)*sin(theta0)^2);
    theta1 = theta1 + atan2((sigma-1)*sin(theta1)*cos(theta1),...
        1+(sigma-1)*sin(theta1)^2);
    return;
%(ii): e <= Vi   
elseif y0p ~= 0
    theta0 = atan(S*y0/y0p);
elseif y0 > 0 
    theta0 = pi/2;
elseif y0 < 0
    theta0 = -pi/2;
else
    theta0=0;
end
if y1p ~= 0
    theta1 = atan(S*y1/y1p);   %ylp is soms Inf
elseif y1 > 0 
    theta1 = pi/2;
elseif y1 < 0
    theta1 = -pi/2;
else
    theta1=0;
end
if y0*y1 >= 0 
    if theta0 > 0 && theta1 < 0
        theta1 = theta1 + pi;
    elseif theta0 < 0 && theta1 > 0
        theta1 = theta1 - pi;
    end
    return
elseif theta0*theta1 > 0 
    theta1 = theta1 + pi;
    return
end


