function meshData = enlargeMesh(slp,mesh,tol,varargin)
%enlarges the mesh of an infinite problem

xmin = slp.xmin;
xmax = slp.xmax;
meshData.LNF=mesh.LNF;
meshData.halfRangeReduction=mesh.halfRangeReduction;
meshData.radial=mesh.radial;
if mesh.radial
    meshData.l=mesh.l;meshData.Vs=mesh.Vs; meshData.r0=mesh.r0; 
end
meshData.trunca=mesh.trunca;
meshData.truncb=mesh.truncb;
meshData.Infa=mesh.Infa;
meshData.Infb=mesh.Infb;
cutoff=slp.classification.cutoff;
contspectrum=slp.classification.cont_spectrum;
wkb=10;
if nargin>3
    e=varargin{1};
    if nargin>4
        wkb=max(varargin{2},10);
    end
else
    e=cutoff;
end

if meshData.LNF %problem in Liouville normal form
    pot = slp.q; %the potential
    %nodes/weights for Gauss-Legendre quadrature
    nodes = [0.095012509837637440185, 0.281603550779258913230,...
       0.458016777657227386342, 0.617876244402643748447,...
       0.755404408355003033895, 0.865631202387831743880,...
       0.944575023073232576078, 0.989400934991649932596];
    weights = [0.189450610455068496285, 0.182603415044923588867,...
       0.169156519395002538189, 0.149595988816576732081,...
       0.124628971255533872052, 0.095158511682492784810,...
       0.062253523938647892863, 0.027152459411754094852];
    P = 16;
    l1= legendre((1-nodes)/2,P);
    l2= legendre((1+nodes)/2,P);
    meshData.imatch=mesh.imatch;
    meshData.h=[];meshData.V0=[];
    cu=[]; cup=[];cv=[];cvp=[];
    rcu=[]; rcup=[];rcv=[];rcvp=[];
%     if slp.classification.cutoff<1e100
%        factor=2;
%     else
%        factor=1.5;
%     end
    if isinf(xmin)
        %add intervals to the front
%         if pot(mesh.Infa)>cutoff && contspectrum(1)
%            x=mesh.Infa*1.2; 
%         elseif e<cutoff && contspectrum(1)
%            x=mesh.Infa*1.1 ;
        if contspectrum(1)
            x=mesh.Infa*1.25;
        elseif pot(mesh.Infa)>cutoff && ~contspectrum(1) %interval probably already or almost large enough at the left side
           x=mesh.Infa*1.025; 
        else
           x = get_xmax(pot,mesh.Infa,cutoff,-1,wkb);%mesh.Infa*factor;
        end
        meshData.Infa=x;
        i=0;
        h=mesh.h(1);
        while x < mesh.Infa
            hnew = 2*h;
            i=i+1;
            it=0;
            while abs(hnew/h-1) > 0.1 && it < 20
                h=max(min([hnew,mesh.Infa-x,1e5]),eps*10);
                [V0,V] = calculateV(pot,x,h);          
                err = estimateError(V);
                hnew = 5*h;
                if err > eps
                 hnew = h * (tol / err)^(1/15);
                end
                if abs(h -(mesh.Infa-x))<eps && hnew>h, break; end
                it=it+1;
            end
            meshData.h(i)=h; meshData.V0(i)=V0;
            [cu(:,i),cup(:,i),cv(:,i),cvp(:,i)] = computeCoeffsBP(V); 
            [rcu(:,i),rcup(:,i),rcv(:,i),rcvp(:,i)] =computeCoeffsRP(V); 
            x = x + h;
            meshData.imatch=meshData.imatch+1;
        end
    end


    %add existing mesh
    meshData.h=[meshData.h mesh.h];
    meshData.V0=[meshData.V0 mesh.V0];
    meshData.cu=[cu'; mesh.cu];
    meshData.cv=[cv'; mesh.cv];
    meshData.cup=[cup'; mesh.cup];
    meshData.cvp=[cvp'; mesh.cvp];
    meshData.rcu=[rcu'; mesh.rcu];
    meshData.rcv=[rcv'; mesh.rcv];
    meshData.rcup=[rcup'; mesh.rcup];
    meshData.rcvp=[rcvp'; mesh.rcvp];

    if isinf(xmax)
         %add intervals to the end
        i=length(meshData.h);
        x = mesh.Infb;
%         if pot(mesh.Infb)>cutoff && contspectrum(2)
%            meshData.Infb=mesh.Infb*1.2;
%         elseif e<cutoff && contspectrum(2)
%            meshData.Infb=mesh.Infb*1.1;
        if contspectrum(2)
            meshData.Infb=mesh.Infb*1.25; 
        elseif pot(mesh.Infb)>cutoff && ~contspectrum(2) %interval probably already large enough at the right side
           meshData.Infb=mesh.Infb*1.025; 
        else
           meshData.Infb= get_xmax(pot,mesh.Infb,cutoff,1,wkb); %mesh.Infb*factor;
        end
        h=mesh.h(end);
        while x < meshData.Infb
            hnew = 2*h;
            i=i+1;
            it=0;
            while abs(hnew/h-1) > 0.1 && it < 20
                h=max(min([hnew,meshData.Infb-x,5e5]),eps*10);
                %[V0,V] = calculateV(pot,x,h);
                a=x; b= x+h;
                m = (b+a)/2;d = (b-a)/2;
                ff=arrayfun(pot,[m+d*nodes m-d*nodes]);
                fpm=ff(1:8);fmm=ff(9:end);
                qpot = d*((weights.*fmm)*l1+(weights.*fpm)*l2);
                V0 = qpot(1)/h;
                V = h*qpot(2:P+1).*(2*(1:P)+1);            
                err = estimateError(V);
                hnew = 5*h;
                if err > eps
                 hnew = h * (tol / err)^(1/15);
                end
                if abs(h -min(meshData.Infb-x,5e5))<eps && hnew>h, break; end
                it=it+1;
            end
            meshData.h(i)=h; meshData.V0(i)=V0; 
            if isinf(V0)|| isnan(V0)
                error(['Sorry, the potential can not be correctly evaluated between ' num2str(x) ' and ' num2str(x+h) ': an infinite or NaN value was obtained.'])
            end
            [cu,cup,cv,cvp] = computeCoeffsBP(V); 
            [rcu,rcup,rcv,rcvp] =computeCoeffsRP(V);
            meshData.cu=[meshData.cu; cu];
            meshData.cv=[meshData.cv; cv];
            meshData.cup=[meshData.cup; cup];
            meshData.cvp=[meshData.cvp; cvp];
            meshData.rcu=[meshData.rcu; rcu];
            meshData.rcv=[meshData.rcv; rcv];
            meshData.rcup=[meshData.rcup; rcup];
            meshData.rcvp=[meshData.rcvp; rcvp];
            x = x + h;
        end
    end
else %~LNF
    meshData.h=[];
    meshData.Q0=[];    meshData.Q1=[];    meshData.Q2=[];
    meshData.W0=[];    meshData.W1=[];    meshData.W2=[];
    meshData.P0=[];    meshData.P1=[];    meshData.P2=[];
    q=slp.q; w=slp.w;
    factor=1.5;
    if isinf(xmin)
        %add intervals to the front
        if q(mesh.Infa)/w(mesh.Infa)>slp.classification.cutoff
           x=mesh.Infa*1.2; 
        else
           x = mesh.Infa*factor;
        end
        meshData.Infa=x;
        i=0;
        h=mesh.h(1);
        while x < mesh.Infa
            hnew = 2*h;
            i=i+1;
            it=0;
            while abs(hnew/h-1) > 0.1 && it < 20
                h=max(min([hnew,mesh.Infa-x,1e5]),eps*10);
                [err,qq,ww,pp] = getErrorSLP(slp.q,slp.p,slp.w,h,x);
                err=h^2*err;        
                hnew = 5*h;
                if(err > eps) 
                    hnew = h * (tol / err)^(1/6);
                end
                if abs(h -(mesh.Infa-x))<eps && hnew>h, break; end
                it=it+1;
            end
            meshData.h(i)=h;
            meshData.Q0(i)=qq(1);meshData.Q1(i)=qq(2);meshData.Q2(i)=qq(3); %store data related to this mesh for later use
            meshData.W0(i)=ww(1);meshData.W1(i)=ww(2);meshData.W2(i)=ww(3);
            meshData.P0(i)=pp(1);meshData.P1(i)=pp(2);meshData.P2(i)=pp(3);
            x = x + h;
        end
    end


    %add existing mesh
    meshData.h=[meshData.h mesh.h];
    meshData.Q0=[meshData.Q0 mesh.Q0];
    meshData.Q1=[meshData.Q1 mesh.Q1];
    meshData.Q2=[meshData.Q2 mesh.Q2];
    meshData.P0=[meshData.P0 mesh.P0];
    meshData.P1=[meshData.P1 mesh.P1];
    meshData.P2=[meshData.P2 mesh.P2];
    meshData.W0=[meshData.W0 mesh.W0];
    meshData.W1=[meshData.W1 mesh.W1];
    meshData.W2=[meshData.W2 mesh.W2];

    if isinf(xmax)
         %add intervals to the end
        i=length(meshData.h);
        x = mesh.Infb;
        if q(mesh.Infb)/w(mesh.Infb)>slp.classification.cutoff
           meshData.Infb=mesh.Infb*1.2;
        else
           meshData.Infb=mesh.Infb*factor;
        end
        h=mesh.h(end);
        while x < meshData.Infb
            hnew = 2*h;
            i=i+1;
            it=0;
             while abs(hnew/h-1) > 0.1 && it < 20
                h=max(min([hnew,meshData.Infb-x,5e5]),eps*10);
                [err,qq,ww,pp] = getErrorSLP(slp.q,slp.p,slp.w,h,x);
                err=h^2*err;        
                hnew = 5*h;
                if(err > eps) 
                    hnew = h * (tol / err)^(1/6);
                end
                if abs(h -min(meshData.Infb-x,5e5))<eps && hnew>h, break; end
                it=it+1;
            end
            meshData.h(i)=h; 
            meshData.Q0(i)=qq(1);meshData.Q1(i)=qq(2);meshData.Q2(i)=qq(3); %store data related to this mesh for later use
            meshData.W0(i)=ww(1);meshData.W1(i)=ww(2);meshData.W2(i)=ww(3);
            meshData.P0(i)=pp(1);meshData.P1(i)=pp(2);meshData.P2(i)=pp(3);
            x = x + h;
        end
    end
    
end

function x=get_xmax(pot,x,cutoff,sign,wkb)
sums = 0;
vmin = pot(x);
he=0.01;
vx_prev=vmin;
x0=x;
while true
       x=x+sign*he;
       vx=pot(x);
       if isnan(vx) %te ver
           x=x-he;
           break;
       end
       if vx > vmin && (vx>vx_prev || (vx-vx_prev)/he >-1e-3)  
          sums=sums+sqrt(vx-vmin)*he;
       else
         vmin=min(vx,cutoff);
         sums=0;
       end
       if sums > wkb || abs(x)>abs(x0)*2 || (abs(vx-cutoff)<1e-15 && abs(vx-vx_prev)/he <1e-10 && x>x0+2*he) || x>x0+max(5e3,abs(x0)*1e2)
         break;
       end
       vx_prev=vx;
       he=min(he*1.2,max(abs(x+1)/10,he));
end
if abs(vx-pot(x+he))<1e-2
    x=x+sign*abs(x)/10;
end



function [V0,V]=calculateV(pot,xmin,h)
% [V0,V]=calculateV(xmin,h)
% this returns the values on [xmin,xmin+h] of V0 and Vbar_i, i=1..p
x = [0.095012509837637440185, 0.281603550779258913230,...
   0.458016777657227386342, 0.617876244402643748447,...
   0.755404408355003033895, 0.865631202387831743880,...
   0.944575023073232576078, 0.989400934991649932596];
w = [0.189450610455068496285, 0.182603415044923588867,...
   0.169156519395002538189, 0.149595988816576732081,...
   0.124628971255533872052, 0.095158511682492784810,...
   0.062253523938647892863, 0.027152459411754094852];
a=xmin; b= xmin+h;
m = (b+a)/2;
d = (b-a)/2;
ff=arrayfun(pot,[m+d*x m-d*x]);
fpm=ff(1:8);
fmm=ff(9:end);
P = 16;
l1= legendre((1-x)/2,P);
l2= legendre((1+x)/2,P);
qpot = d*((w.*fmm)*l1+(w.*fpm)*l2);
V0 = qpot(1)/h;
V = h*qpot(2:P+1).*(2*(1:P)+1);









