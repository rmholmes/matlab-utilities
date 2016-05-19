function [x,y,yp,sq] = computeEigenfunction(slpObject,eigenvalue,mesh,varargin)
% COMPUTEEIGENFUNCTION 
% [x,y,yp,sq] = computeEigenfunction(slpObject,eigenvalue,mesh[,evalPoints])
% returns eigenfunction corresponding to the eigenvalue e
%   Input arguments:
%                   eigenvalue: scalar value equal to the eigenvalue or 
%                   a vector with eigenvalue(1) equal to the index of the eigenvalue
%                   and eigenvalue(2) the actual value of the eigenvalue
%                   slpObject: the slp-problem for which an eigenfunction is
%                               computed (an instance of the class slp)
%                   mesh: a structure containing information on the mesh
%                         which was created and returned by
%                         computeEigenvalues
%                   evalPoints: optional argument, array containing extra
%                       x-values where the user wants the eigenfunction to be
%                       evaluated
%   Output arguments:
%                   x,y,yp: arrays of length equal to the number of
%                           meshpoints if evalPoints is empty or not provided
%                   x,y,yp: arrays of length equal to the number of
%                           meshpoints + the length of evalPoints
%                   x: x-values where the eigenfunction was evaluated
%                   y: evaluations of the eigenfunction in the x-values
%                      returned in the x-array
%                   yp: corresponding first order derivatives py'
%                   sq: (optional) estimation of the norm

e=eigenvalue(end); %actual value of the eigenvalue
if slpObject.liouvilleTransformed
    %problem has undergone a Liouville transformation
     slpObject.liouvilleTransformed=false;
     if nargin>3
        %convert "SLP"-evaluation points to the corresponding ones for the problem
        % in Liouville normal form
        tmp=sort(unique(varargin{1}));
        if tmp(1)<slpObject.x2r(slpObject.xmin)-1e-10 || tmp(end)>slpObject.x2r(slpObject.xmax)+1e-10
            error('computeEigenfunction: user-specified evaluation points need to be situated in the integration interval')
        end
        tmp(1)=max(tmp(1),slpObject.x2r(slpObject.xmin));
        tmp(end)=min(tmp(end),slpObject.x2r(slpObject.xmax));
        nodes=slpObject.r2x(tmp); 
        [x,y,yp,sq] = slpObject.computeEigenfunction(e,mesh,nodes');   
     else
        [x,y,yp,sq] = slpObject.computeEigenfunction(e,mesh);   
     end
     % compute the original Sturm-Liouville-eigenfunction
     % out of the Schrodinger-eigenfunction:
     [x,y,yp] = transformEigenfunction(slpObject,x,y,yp);
     tmp=find(~isnan(y));
     x=x(tmp); y=y(tmp); yp=yp(tmp);
     return;
end


if mesh.halfRangeReduction &&  slpObject.xmin==-slpObject.xmax
    %the problem is symmetric and half-range reduction has been applied
    % solving on the half range may make a problem more tractable
    
    
    %the half range problems:
    slpObject2=slp(slpObject.p,slpObject.q,slpObject.w,0,slpObject.xmax,0,1,slpObject.a1,slpObject.b1,slpObject.jumps(slpObject.jumps>0));
    slpObject3=slp(slpObject.p,slpObject.q,slpObject.w,0,slpObject.xmax,1,0,slpObject.a1,slpObject.b1,slpObject.jumps(slpObject.jumps>0));
    %is the eigenvalue even or odd?
    even=false;
    if length(eigenvalue)>1
        indx=eigenvalue(1);
        if rem(indx,2)==0
            even=true;
        end
    else
        pruf2=locateEigenvalue(slpObject2,e,mesh,true);
        pruf2= abs(round(pruf2)-pruf2);
        pruf3=locateEigenvalue(slpObject3,e,mesh,true);
        pruf3= abs(round(pruf3)-pruf3);
        if pruf2 < pruf3
            even=true;
        end
    end
    nodes=[];
    if nargin>3 && ~isempty(varargin{1})
        nodes=sort(unique(varargin{1}));
        nodes1=nodes(nodes>=0);
        nodes2=abs(nodes(nodes<=0));
        nodes=sort(unique([nodes1 nodes2]));
    end
    if even %even eigenvalue
        [x,y,yp,sq] = slpObject2.computeEigenfunction(e,mesh,nodes);
        % compose eigenfunction information over the full interval using
        % the computed eigenfunction of the half-range problem:
        if x(1)==0
            x=[-x(length(x):-1:1) x(2:end)];
            y=[y(length(y):-1:1) y(2:end)]/sqrt(2);
            yp=[-yp(length(yp):-1:1) yp(2:end)]/sqrt(2);
        else
            x=[-x(length(x):-1:1) x];
            y=[y(length(y):-1:1) y]/sqrt(2);
            yp=[-yp(length(yp):-1:1) yp]/sqrt(2);
        end
    else %odd eigenvalue
        [x,y,yp,sq] = slpObject3.computeEigenfunction(e,mesh,nodes);
        % compose eigenfunction information over the full interval using
        % the computed eigenfunction of the half-range problem:
        if x(1)==0
            x=[-x(length(x):-1:1) x(2:end)];
            y=[-y(length(y):-1:1) y(2:end)]/sqrt(2);
            yp=[yp(length(yp):-1:1) yp(2:end)]/sqrt(2);
        else
            x=[-x(length(x):-1:1) x];
            y=[-y(length(y):-1:1) y]/sqrt(2);
            yp=[yp(length(yp):-1:1) yp]/sqrt(2);
        end
    end
    for i=1:length(y) 
        if y(i) > 0.1
            break;
        end
        if y(i) < -0.1
            y = (-1).*y;
            yp = (-1).*yp;
            break;
        end
    end
    y=real(y);
    yp=real(yp);
    return;
end

if nargin>3 && ~isempty(varargin{1})
     tmp=sort(unique(varargin{1}));
     if tmp(1)<slpObject.xmin-1e-13 || tmp(end)>slpObject.xmax+1e-13
            error('computeEigenfunction: user-specified evaluation points need to be situated in the integration interval')
     end
     tmp(1)=max(tmp(1),slpObject.xmin);
     tmp(end)=min(tmp(end),slpObject.xmax);
     % meshpts = mesh points + user-specified evaluation points
    %we also use the original mesh points in order to make sure that the
    %step lengths in the new mesh are not too large
    meshpts=sort(unique([mesh.Infa mesh.Infa+cumsum(mesh.h) varargin{1}]));
    mesh2.radial=mesh.radial;
    mesh2.LNF=mesh.LNF;
    mesh2.halfRangeReduction=false;
    if mesh.LNF % problem in Liouville normal form
        pot = slpObject.q;
        %nodes and weigths for Gauss-Legendre quadrature:
        nodes = [0.095012509837637440185, 0.281603550779258913230,...
           0.458016777657227386342, 0.617876244402643748447,...
           0.755404408355003033895, 0.865631202387831743880,...
           0.944575023073232576078, 0.989400934991649932596];
        weigths = [0.189450610455068496285, 0.182603415044923588867,...
           0.169156519395002538189, 0.149595988816576732081,...
           0.124628971255533872052, 0.095158511682492784810,...
           0.062253523938647892863, 0.027152459411754094852];
        P = 16;
        % auxiliary data for Gauss-Legendre quadrature:
        l1= legendre((1-nodes)/2,P);
        l2= legendre((1+nodes)/2,P);
        % construct mesh2, i.e. a mesh with meshpoints formed by the
        % original meshpoints supplemented with the evaluation points
        % specified by the user (as last input argument):
        if mesh.radial
            mesh2.l=mesh.l;
            mesh2.r0=mesh.r0;
            mesh2.Vs=mesh.Vs;
            meshpts0=meshpts(meshpts<mesh2.r0);
            meshpts=meshpts(meshpts>=mesh2.r0);
            mesh2.Infa=meshpts(1);
        end
        mesh2.Infa=meshpts(1);

        for i=1:length(meshpts)-1
            x=meshpts(i);
            h=meshpts(i+1)-x;
            m = (x+x+h)/2; d = h/2;
            ff=arrayfun(pot,[m+d*nodes m-d*nodes]);
            fpm=ff(1:8); fmm=ff(9:end);
            qpot = d*((weigths.*fmm)*l1+(weigths.*fpm)*l2);
            mesh2.V0(i) = qpot(1)/h;
            mesh2.h(i)=h;
            V = h*qpot(2:P+1).*(2*(1:P)+1);
            [mesh2.cu(i,:),mesh2.cup(i,:),mesh2.cv(i,:),mesh2.cvp(i,:)] = computeCoeffsBP(V); 
        end
        [~,mesh2.imatch] = min(mesh2.V0);
    else % not in Liouville normal form
        for i=1:length(meshpts)-1
            x=meshpts(i);
            h=meshpts(i+1)-x;
            [~,qq,ww,pp] = getErrorSLP(slpObject.q,slpObject.p,slpObject.w,h,x);
            mesh2.h(i)=h;
            mesh2.Q0(i)=qq(1); mesh2.Q1(i)=qq(2); mesh2.Q2(i)=qq(3); %store data related to this mesh for later use
            mesh2.W0(i)=ww(1); mesh2.W1(i)=ww(2); mesh2.W2(i)=ww(3);
            mesh2.P0(i)=pp(1); mesh2.P1(i)=pp(2); mesh2.P2(i)=pp(3);
            mesh2.Infa=meshpts(1);
        end    
    end
    % compute the eigenfunction over the new mesh:
    [x,y,yp,sq] = slpObject.computeEigenfunction(e,mesh2);
    if length(meshpts)==length(x)
      x=meshpts; %to avoid that values in x are no longer exactly equal to user-defined meshpoints due to rounding errors;
    else
      x=[x(1) meshpts];  %radial problem
    end
    if mesh.radial
        %meshpoints may need to be added before mesh.r0
        g=@(x)arrayfun(slpObject.q,x).*x.^2;
        xt=[(1/3)*sqrt(5-2*sqrt(10/7)) (1/3)*sqrt(5+2*sqrt(10/7))]; %5 Legendre points in [0,r0]
        warning('off','all');
        yn=[];ynp=[];
        for mp=meshpts0
             m = mp/2;
             t=[m+m*xt m m-m*xt];
             p=polyfit(t,g(t),4);
             part.l=-1/2+sqrt(1+4*p(5))/2;
             part.Vs=fliplr(p(1:4));
             part.r0=mp;
             [y0,y0e]=initstep(slpObject,e,part);
             yn=[yn y0(1)]; ynp=[ynp y0(2)];
        end
        if x(1)<mesh.r0
            x=x(2:end); y=y(2:end); yp=yp(2:end);
        end
        %find scaling
        mp=x(1);
        m = mp/2;
        t=[m+m*xt m m-m*xt];
        p=polyfit(t,g(t),4);
        part.l=-1/2+sqrt(1+4*p(5))/2;
        part.Vs=fliplr(p(1:4));
        part.r0=mp;
        [y0,y0e]=initstep(slpObject,e,part);
        if abs(y0(1))>1e-10
            yn=yn*y(1)/y0(1); ynp=ynp*y(1)/y0(1);
        else
            yn=yn*yp(1)/y0(2); ynp=ynp*yp(1)/y0(2);
        end
        x=[meshpts0 x];
        y=[yn y];
        yp=[ynp yp];
        warning('on','all');
    end
    %remove mesh points from the original mesh, in order to retain only the
    %user-specified ones
    indices=[];
    for p=tmp
      indices=[indices find(abs(x-p)<eps)];
    end
    x=x(indices);
    y=real(y(indices));
    yp=real(yp(indices));
    return
end


%no extra evaluation points were passed as input argument

%initial values:
a0 = slpObject.a0;
b0 = slpObject.b0;
a1 = slpObject.a1;
b1 = slpObject.b1;
sq1 = sqrt(a0^2+b0^2);
sq2 = sqrt(a1^2+b1^2);
%initial value for left hand solution
if mesh.radial
    [yi,yei]=initstep(slpObject,e,mesh);
else
    yi = [-b0 / sq1;a0 / sq1];
    yei=[0; 0];
end
%initial value for right hand solution
yf = [-b1 / sq2;a1 / sq2];
%initial values for derivatives wrt E (eigenvalue):
yef=[0; 0];

%propagate left hand solution to matching point
[yl,yel] = propagateSolution(e,yi,yei,mesh,true,false);
%idem right hand solution
[yr,yer] = propagateSolution(e,yf,yef,mesh,false,false);
%computation of the norm
qr = [yer(2) -yer(1)]*yr;  %De formule blijft dezelfde voor SLP!!!
ql = [-yel(2) yel(1)]*yl;
  
%make sure left and right hand solution match:
if abs(yr(1)) < 1e-10 || (abs(yr(2))+abs(yl(2)) - (abs(yr(1))+abs(yl(1))) > (abs(yr(2))+abs(yl(2))))  
     gap = yl(2)/yr(2);
else
    gap = yl(1)/yr(1);
end
fnorm = ql+qr*gap*gap;
sqfnorm = 1;
if fnorm  <= 0 
	warning('warning:eigenfunction',' calculate_eigenfunction: There are problems with the normalization.')	
else	
    sqfnorm = sqrt(fnorm);%*sigmaL;  
end

yi = yi./sqfnorm;
yf = yf*gap/sqfnorm;

[ynL,qL] = propagateSolution2(e,yi,yei,mesh,true);
[ynR,qR] = propagateSolution2(e,yf,[0;0],mesh,false);
if ynL(2,end)*ynR(2,1)<0 && abs(ynL(2,end))>1e-1
    ynR=-ynR;
end
if ynL(1,end)*ynR(1,1)<0 && abs(ynL(1,end))>1e-1
    ynR=-ynR;
end
yn=[ynL ynR(:,2:end)];
if mesh.radial
  q=1;  %Besselnormalform0
else
  q=qL-qR;
end
if q <= 0
   sq = 1;
   warning('warning:eigenfunction1',' calculate_eigenfunction: There are problems with the normalization.' )
else	
   sq=1/sqrt(q);
end

yn = yn.*sq;
if mesh.radial
    x=[mesh.r0 mesh.r0+cumsum(mesh.h)];
    if mesh.r0>1e-5
    %add extra meshpoint near origin
        g=@(x)arrayfun(slpObject.q,x).*x.^2;
        r0=1e-5;
        m = r0/2;
        xt=[(1/3)*sqrt(5-2*sqrt(10/7)) (1/3)*sqrt(5+2*sqrt(10/7))]; %5 Legendre points in [0,r0]
        t=[m+m*xt m m-m*xt];
        warning('off','all');
        p=polyfit(t,g(t),4);
        warning('on','all');
        part.l=-1/2+sqrt(1+4*p(5))/2;
        part.Vs=fliplr(p(1:4));
        part.r0=r0;
        [y0,y0e]=initstep(slpObject,e,part);
        x=[r0 x];
        yn=[y0*sq./sqfnorm yn];
    end
    
else
    x=[mesh.Infa mesh.Infa+cumsum(mesh.h)];
end

for i=1:size(yn,2) 
    if yn(1,i) > 0.1
        break;
    end
    if yn(1,i) < -0.1
        yn = (-1).*yn;
        break;
    end
end

y = real(yn(1,:));
yp = real(yn(2,:));
end




function [yn,q] = propagateSolution2(e,y,ye,part,forward)
if part.LNF
    imatch=part.imatch;
else
    [~,imatch]=max((e*part.W0-part.Q0).*part.P0);
    if imatch==1 
        [~,imatch]=max((e*part.W0(2:end)-part.Q0(2:end)).*part.P0(2:end));
        if imatch>=length(part.h)-1
           imatch=ceil(length(part.h)/2);
        end
    end
    if imatch==length(part.h)
        [~,imatch]=max((e*part.W0(1:imatch-1)-part.Q0(1:imatch-1)).*part.P0(1:imatch-1));
        if imatch==1
           imatch=ceil(length(part.h)/2);
        end
    end
end    
if forward
 t = 1:imatch;
else
 t = length(part.h):-1:imatch+1;    
end

yn=zeros(2,length(t));
yn(:,1)=y;
q=[y(2) -y(1)]*ye;


hm=part.h;
%wght1(q+1) = Coefficient of Z^q in series expansion for eta_m.
%wght2(q+1) = Coefficient of Z^q in series expansion for eta_{m-1}.
wght1 =[1.527349308567059e-009 0.36365459727787e-10 0.00395276736172e-10 0.00002635178241e-10 0.00000012199899e-10 0.00000000042069e-10...
            0.00000000000113e-10 0];
wght2 = [2.901963686277412e-008  0.76367465428353e-9 0.00909136493195e-9 0.00006587945603e-9 0.00000032939728e-9 0.00000000121999e-9...
            0.00000000000351e-9 0.00000000000001e-9];
mmaxm=7;
if part.LNF
    cum=part.cu;
    cupm=part.cup;
    cvm=part.cv;
    cvpm=part.cvp;
    vbarm=part.V0;
    M = size(cum,2);
    eta=zeros(1,M+1);
    for i = 1:length(t)
        tt=t(i);
        cu=cum(tt,:);
        cup=cupm(tt,:);
        cv=cvm(tt,:);
        cvp=cvpm(tt,:);
        vbar=vbarm(tt);
        h=hm(tt);
        Z=(vbar-e)*h^2;
        %calculation of eta functions:
        if abs(Z) < 0.5 
           tmp=(Z.^(0:mmaxm))';
           eta(10) = wght1*tmp;
           eta(9) = wght2*tmp;
           for k=(8):-1:1
                 eta(k) = Z * eta(k+2) + (2*k+1) * eta(k+1);
           end
           csi = Z * eta(2) + eta(1);
        else
            if Z>1e3 
               %to avoid infinite values
               sq = sqrt(Z);
               csi = 1/2+exp(-2*sq)/2; %=cosh(sq)/exp(sq)    
               eta(1) = (1/2-exp(-2*sq)/2)/sq; %eta1/exp(sq) 
           elseif Z >= 0.5
            sq = sqrt(Z);
            csi = cosh(sq);    
            eta(1) = sinh(sq)/sq; 
           else
            sq = sqrt(abs(Z));
            csi = cos(sq);    
            eta(1) = sin(sq)/sq;
           end
           eta(2)=(csi-eta(1))/Z;
           for k=3:M+1
              eta(k)=(eta(k-2)-(2*k-3)*eta(k-1))/Z;
           end
        end    
       %end calculation eta functions.    
        tmp = eta(1:M)';
        eta1=eta(1);
        u  = csi + cu*tmp;
        up = (Z * eta1 + cup*tmp)/h;
        v = (eta1 + cv*tmp)*h;
        vp = csi + cvp*tmp;
        %using eq.(A.8):  
        tmp=eta(2:M+1)';
        ue = eta1 + cu*tmp;  
        uep = csi + eta1 + cup*tmp;
        ve = eta(2) + cv*tmp;
        vep = eta1 + cvp*tmp;
        h2=h*h;
        ue = -ue*h2/2; 
        uep = -uep*h/2; 
        ve = -ve*(h2*h)/2; 
        vep = -vep*h2/2; 
        if forward
            T=[u, v ; up, vp];
            Te=[ue, ve ; uep, vep];
        else
            T=[vp, -v ; -up, u];
            Te=[vep, -ve ; -uep, ue];
        end
        if e<vbar-eps && Z<1e3
            sigma = exp(sqrt(vbar-e)*h); %cfr.SLEDGE
        else
           sigma = 1;
        end      
        ye = Te*y/sigma; %ye=[0;0];
        y = T*y/sigma;
        yn=yn/sigma;
        yn(:,i+1)=y; 
        q=q/sigma^2;
        q = q + [yn(2,i+1) -yn(1,i+1)]*ye;
    end
else
    M=5;
    for i = 1:length(t)
            tt=t(i);
            h=hm(tt);
            P0=part.P0(tt); P1=part.P1(tt); P2=part.P2(tt);
            QW1=part.Q1(tt)-e*part.W1(tt); QW0=(part.Q0(tt)-e*part.W0(tt));
            Z=QW0*P0*h^2;
            if abs(Z) < 0.5          
               tmp=(Z.^(0:mmaxm))';
               eta(10) = wght1*tmp;
               eta(9) = wght2*tmp;
               for k=(8):-1:1
                 eta(k) = Z * eta(k+2) + (2*k+1) * eta(k+1);
               end
               csi = Z * eta(2) + eta(1);
            else
               if Z>1e3 
                   %to avoid infinite values
                   sq = sqrt(Z);
                   csi = 1/2+exp(-2*sq)/2; %=cosh(sq)/exp(sq)    
                   eta(1) = (1/2-exp(-2*sq)/2)/sq; %eta1/exp(sq) 
               elseif Z >= 0.5
                sq = sqrt(Z);
                csi = cosh(sq);    
                eta(1) = sinh(sq)/sq; 
               else
                sq = sqrt(abs(Z));
                csi = cos(sq);    
                eta(1) = sin(sq)/sq;
               end
               eta(2)=(csi-eta(1))/Z;
               for k=3:M+1
                  eta(k)=(eta(k-2)-(2*k-3)*eta(k-1))/Z;
               end
            end
            eta0=eta(1);  eta1=eta(2); eta2=eta(3); eta3=eta(4); eta4=eta(5); eta5=eta(6);

            R1=QW1; R0=QW0;
            U1=R1*P0-P1*R0;
            V1=R1*P0+P1*R0;

            wP1=P1/(h^2*P0);
          

            R2=part.Q2(tt)-e*part.W2(tt);
            U2=R2*P0-P2*R0;
            V2=R2*P0+P2*R0;
            wP2=P2/(h^2*P0);

            T(1,1)=csi+((U1*wP1)/12+(U2*wP2)/20)*eta0+(-(U1)/2-(wP1*(V2+U1)+wP2*V1)/4)*eta1+...
                ((3*wP1*V2+wP2*(7*V1-3*U2))/4-(U2*V2)/40-(U1*V1)/24)*eta2+(5*V1*(U2+2*V2)+V2*(5*U1-3*U2))/40*eta3;

            T(1,2)=h*P0*((1-(wP1^2)/2+wP2)*eta0+((3*(wP1^2+wP2^2))/2-3*wP2+(wP1*(U1+V1))/12+(wP2*(U2+V2))/20)*eta1-...
                ((V2)/2+(15*wP2^2)/2+(wP1*(U1+V1))/6+(wP2*(2*U2+7*V2))/10)*eta2+...
                (-(V2^2)/40-(V1^2)/24+(wP2*(9*U2+24*V2))/10)*eta3+(9*V2^2)/40*eta4);

            T(2,1)=(((wP1*U1)/12+(wP2*U2)/20)*csi+(Z+(U2)/2+(wP2*(V2-U2+U1))/4+(wP1*(V1-U2))/4)*eta0+...
                (-(3*U2)/2-(U2*V2)/40-(U1*V1)/24-(wP2*(9*V2-3*U2+3*U1))/4-(wP1*(3*V1+U1-3*U2))/4)*eta1+...
                (-(U1*V1)/24+(U2*V2)/10-(U1*(V1+V2)+U2*(V2-V1)+V1^2+V2^2)/8-(3*wP2*(U2-10*V2))/4)*eta2+((V2*(30*V2+27*U2))/40)*eta3)/h/P0;

            T(2,2)=(1+(wP1*(wP1+2*wP2-2))/2)*csi+((1-4*wP2)*wP1+((V1+U1)*wP1)/12+((V2+U2)*wP2)/20-(wP1^2-3*wP2^2)/2)*eta0+...
               (((1-wP1)*V1)/2-(9*wP2*(wP2-2*wP1))/2+(wP1*V2+wP2*U1)/2)*eta1-((3*wP2*(V2+2*V1))/2+(wP1*(V2+3*U2))/2+...
               (V2^2)/40+(V1^2)/24)*eta2-(V2*(3*V2+20*V1))/40*eta3;

            Te(1,1)=(eta0+((U1*wP1)/12+(U2*wP2)/20)*eta1+(-(U1)/2-(wP1*(V2+U1)+wP2*V1)/4)*eta2+...
                ((3*wP1*V2+wP2*(7*V1-3*U2))/4-(U2*V2)/40-(U1*V1)/24)*eta3+(5*V1*(U2+2*V2)+V2*(5*U1-3*U2))/40*eta4);

            Te(1,2)=h*P0*((1-(wP1^2)/2+wP2)*eta1+((3*(wP1^2+wP2^2))/2-3*wP2+(wP1*(U1+V1))/12+(wP2*(U2+V2))/20)*eta2-...
                ((V2)/2+(15*wP2^2)/2+(wP1*(U1+V1))/6+(wP2*(2*U2+7*V2))/10)*eta3+...
                (-(V2^2)/40-(V1^2)/24+(wP2*(9*U2+24*V2))/10)*eta4+(9*V2^2)/40*eta5);

            Te(2,1)=(((wP1*U1)/12+(wP2*U2)/20)*eta0+(Z+(U2)/2+(wP2*(V2-U2+U1))/4+(wP1*(V1-U2))/4)*eta1+2*eta0+...
                (-(3*U2)/2-(U2*V2)/40-(U1*V1)/24-(wP2*(9*V2-3*U2+3*U1))/4-(wP1*(3*V1+U1-3*U2))/4)*eta2+...
                (-(U1*V1)/24+(U2*V2)/10-(U1*(V1+V2)+U2*(V2-V1)+V1^2+V2^2)/8-(3*wP2*(U2-10*V2))/4)*eta3+((V2*(30*V2+27*U2))/40)*eta4)/h/P0;

            Te(2,2)=(1+(wP1*(wP1+2*wP2-2))/2)*eta0+((1-4*wP2)*wP1+((V1+U1)*wP1)/12+((V2+U2)*wP2)/20-(wP1^2-3*wP2^2)/2)*eta1+...
               (((1-wP1)*V1)/2-(9*wP2*(wP2-2*wP1))/2+(wP1*V2+wP2*U1)/2)*eta2-((3*wP2*(V2+2*V1))/2+(wP1*(V2+3*U2))/2+...
               (V2^2)/40+(V1^2)/24)*eta3-(V2*(3*V2+20*V1))/40*eta4;
         
            Te=-Te*part.W0(tt)*P0*h^2/2;
                       
            if QW0*P0 > eps
                   sigma = exp(sqrt(QW0*P0)*h); %cfr.SLEDGE
            else
                   sigma = 1;
            end

            if forward
             ye = Te*y/sigma; %ye=[0;0];
             y = T*y/sigma;
            else
             iT(1,1)=T(2,2); iT(1,2)=-T(1,2); iT(2,1)=-T(2,1); iT(2,2)=T(1,1);
             iTe(1,1)=Te(2,2); iTe(1,2)=-Te(1,2); iTe(2,1)=-Te(2,1); iTe(2,2)=Te(1,1); 
             ye = iTe*y/sigma;
             y=iT*y/sigma;            
            end
            yn=yn/sigma;
            yn(:,i+1)=y; 
            q=q/sigma^2;
            q = q + [yn(2,i+1) -yn(1,i+1)]*ye;     
    end
end
if ~forward
    yn=fliplr(yn);
end
end


function [r,yr,ypr] = transformEigenfunction(slpObject,x,y,yp)
% calculating the original Sturm-Liouville-eigenfunction 'z'
% out of the Schrodinger-eigenfunction 'y'

% x = vector with the values of the meshpoints
% r = idem but for original variable

r = slpObject.x2r(x);
p = slpObject.SLPp;
w = slpObject.SLPw;
m=zeros(1,length(x));

function x=mp(r)
    % using automatic differentiation
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
         mr=(p(r).*w(r))^(-1/4);
         x=-0.25*mr.^5*(dp.*tmpw(1)+tmpp(1).*dw); 
end
    

md=zeros(1,length(m));
mt=zeros(1,length(m));
for i = 1:length(m)
          ri=r(i);  
          pr=p(ri);
          m(i) = (pr*(w(ri)))^(-1/4);
          md(i) = mp(ri);
          mt(i) = m(i)*pr;
end
yr = y .* m;
ypr = arrayfun(p,r).*(md .* y + yp./(mt)); % py' en niet y' !!!
end


