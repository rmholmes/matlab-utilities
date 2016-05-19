function meshData = refineSingularMesh(slp,mesh,tol)
%is used to refine a mesh corresponding to a singular problem, this means
%that the truncated endpoints will move closer to the singular endpoint

meshData.LNF=mesh.LNF;
meshData.halfRangeReduction=mesh.halfRangeReduction;
meshData.radial=mesh.radial;

if meshData.LNF

    pot = slp.q; %the potential
    xmin = slp.xmin;
    xmax = slp.xmax;

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

    %refine first interval of mesh :
    x = xmin;
    meshData.h=[];
    meshData.V0=[];
    meshData.cu=[];
    meshData.cv=[];
    meshData.cup=[];
    meshData.cvp=[];
    meshData.rcu=[];
    meshData.rcv=[];
    meshData.rcup=[];
    meshData.rcvp=[];
    meshData.trunca=mesh.trunca;
    if mesh.trunca>0
        %split interval into two intervals
        i=1;
        h=mesh.h(1)/2;
        m = (x+x+h)/2; d = h/2;
        ff=arrayfun(pot,[m+d*nodes m-d*nodes]); %arrayfun: no problems due to forgotten vectorization in function definitions
        fpm=ff(1:8); fmm=ff(9:end);
        qpot = d*((weights.*fmm)*l1+(weights.*fpm)*l2);
        V0 = qpot(1)/h;
        V = h*qpot(2:P+1).*(2*(1:P)+1);
        meshData.h(i)=h; meshData.V0(i)=V0;
        [cu(:,i),cup(:,i),cv(:,i),cvp(:,i)] = computeCoeffsBP(V); 
        [rcu(:,i),rcup(:,i),rcv(:,i),rcvp(:,i)] =computeCoeffsRP(V);
        x=x+h;
        xend=x+h;
        while x < xend
            hnew = 2*h;
            i=i+1;
            it=0;
            while abs(hnew/h-1) > 0.1 && it < 20
                h=max(min(hnew,xend-x),eps*10);
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
                if abs(h -min(xend-x,5e5))<eps && hnew>h, break; end
                it=it+1;
            end
            meshData.h(i)=h; meshData.V0(i)=V0; 
            [cu(:,i),cup(:,i),cv(:,i),cvp(:,i)] = computeCoeffsBP(V); 
            [rcu(:,i),rcup(:,i),rcv(:,i),rcvp(:,i)]=computeCoeffsRP(V);
            x = x + h;
        end
        meshData.cu=cu';  %omdat kolommen toevoegen sneller gaat dan rijen toevoegen
        meshData.cup=cup'; meshData.cv=cv'; meshData.cvp=cvp';
        meshData.rcu=rcu';  meshData.rcup=rcup'; meshData.rcv=rcv'; meshData.rcvp=rcvp';
    end
    meshData.imatch=mesh.imatch+max(0,length(meshData.h)-1);
    %add middle part of mesh
    meshData.h=[meshData.h mesh.h(mesh.trunca+1:length(mesh.h)-mesh.truncb)];
    meshData.V0=[meshData.V0 mesh.V0(mesh.trunca+1:length(mesh.h)-mesh.truncb)];
    meshData.cu=[meshData.cu; mesh.cu(mesh.trunca+1:length(mesh.h)-mesh.truncb,:)];
    meshData.cv=[meshData.cv; mesh.cv(mesh.trunca+1:length(mesh.h)-mesh.truncb,:)];
    meshData.cup=[meshData.cup; mesh.cup(mesh.trunca+1:length(mesh.h)-mesh.truncb,:)];
    meshData.cvp=[meshData.cvp; mesh.cvp(mesh.trunca+1:length(mesh.h)-mesh.truncb,:)];
    meshData.rcu=[meshData.rcu; mesh.rcu(mesh.trunca+1:length(mesh.h)-mesh.truncb,:)];
    meshData.rcv=[meshData.rcv; mesh.rcv(mesh.trunca+1:length(mesh.h)-mesh.truncb,:)];
    meshData.rcup=[meshData.rcup; mesh.rcup(mesh.trunca+1:length(mesh.h)-mesh.truncb,:)];
    meshData.rcvp=[meshData.rcvp; mesh.rcvp(mesh.trunca+1:length(mesh.h)-mesh.truncb,:)];

    %truncb-part needs to refined and added
    x = xmin+sum(meshData.h);
    i=0;
    meshData.truncb=mesh.truncb;
    hs=[]; V0s=[]; cu=[]; cup=[]; cv=[]; cvp=[]; rcu=[]; rcup=[]; rcv=[]; rcvp=[];
    for g=length(mesh.h)-mesh.truncb+1:length(mesh.h)
        %split interval g into two intervals
         h=mesh.h(g)/2;
         xend=x+h;
         while x < xend
            hnew = 2*h;
            i=i+1;
            it=0;
            while abs(hnew/h-1) > 0.1 && it < 20
                h=max(min(hnew,xend-x),eps*10);
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
                if abs(h -min(xend-x,5e5))<eps && hnew>h, break; end
                it=it+1;
            end
            hs(i)=h; V0s(i)=V0; 
            [cu(:,i),cup(:,i),cv(:,i),cvp(:,i)] = computeCoeffsBP(V); 
            [rcu(:,i),rcup(:,i),rcv(:,i),rcvp(:,i)]=computeCoeffsRP(V);
            x = x + h;
        end
        i=i+1; 
        h=mesh.h(g)/2; 
        m = (x+x+h)/2; d = h/2;
        ff=arrayfun(pot,[m+d*nodes m-d*nodes]); %arrayfun: no problems due to forgotten vectorization in function definitions
        fpm=ff(1:8); fmm=ff(9:end);
        qpot = d*((weights.*fmm)*l1+(weights.*fpm)*l2);
        V0 = qpot(1)/h;
        V = h*qpot(2:P+1).*(2*(1:P)+1);
        hs(i)=h; V0s(i)=V0;
        [cu(:,i),cup(:,i),cv(:,i),cvp(:,i)] = computeCoeffsBP(V); 
        [rcu(:,i),rcup(:,i),rcv(:,i),rcvp(:,i)] =computeCoeffsRP(V);
    end
    meshData.h=[meshData.h hs];
    meshData.V0=[meshData.V0 V0s];
    meshData.cu=[meshData.cu; cu'];
    meshData.cv=[meshData.cv; cv'];
    meshData.cup=[meshData.cup; cup'];
    meshData.cvp=[meshData.cvp; cvp'];
    meshData.rcu=[meshData.rcu; rcu'];
    meshData.rcv=[meshData.rcv; rcv'];
    meshData.rcup=[meshData.rcup; rcup'];
    meshData.rcvp=[meshData.rcvp; rcvp'];
else
    
    xmin = slp.xmin;
   
    %refine first interval of mesh :
    x = xmin;
    meshData.h=[];
    meshData.Q0=[];    meshData.Q1=[];    meshData.Q2=[];
    meshData.W0=[];    meshData.W1=[];    meshData.W2=[];
    meshData.P0=[];    meshData.P1=[];    meshData.P2=[];
    meshData.trunca=mesh.trunca;
    if mesh.trunca>0
        %split interval into two intervals
        i=1;
        h=mesh.h(1)/2;
        [~,qq,ww,pp] = getErrorSLP(slp.q,slp.p,slp.w,h,x);
        meshData.h(i)=h; 
        meshData.Q0(i)=qq(1);meshData.Q1(i)=qq(2);meshData.Q2(i)=qq(3); %store data related to this mesh for later use
        meshData.W0(i)=ww(1);meshData.W1(i)=ww(2);meshData.W2(i)=ww(3);
        meshData.P0(i)=pp(1);meshData.P1(i)=pp(2);meshData.P2(i)=pp(3);        
        x=x+h;
        xend=x+h;
        while x < xend
            hnew = 2*h;
            i=i+1;
            it=0;
            while abs(hnew/h-1) > 0.1 && it < 20
                h=max(min(hnew,xend-x),eps*10);
                [err,qq,ww,pp] = getErrorSLP(slp.q,slp.p,slp.w,h,x);
                err=h^2*err;        
                hnew = 5*h;
                if(err > eps) 
                    hnew = h * (tol / err)^(1/6);
                end
                if abs(h -min(xend-x,5e5))<eps && hnew>h, break; end
                it=it+1;
            end
            meshData.h(i)=h;
            meshData.Q0(i)=qq(1);meshData.Q1(i)=qq(2);meshData.Q2(i)=qq(3); %store data related to this mesh for later use
            meshData.W0(i)=ww(1);meshData.W1(i)=ww(2);meshData.W2(i)=ww(3);
            meshData.P0(i)=pp(1);meshData.P1(i)=pp(2);meshData.P2(i)=pp(3);
            x = x + h;
        end
    end
    %add middle part of mesh
    meshData.h=[meshData.h mesh.h(mesh.trunca+1:length(mesh.h)-mesh.truncb)];
    meshData.Q0=[meshData.Q0 mesh.Q0(mesh.trunca+1:length(mesh.h)-mesh.truncb)];
    meshData.Q1=[meshData.Q1 mesh.Q1(mesh.trunca+1:length(mesh.h)-mesh.truncb)];
    meshData.Q2=[meshData.Q2 mesh.Q2(mesh.trunca+1:length(mesh.h)-mesh.truncb)];
    meshData.W0=[meshData.W0 mesh.W0(mesh.trunca+1:length(mesh.h)-mesh.truncb)];
    meshData.W1=[meshData.W1 mesh.W1(mesh.trunca+1:length(mesh.h)-mesh.truncb)];
    meshData.W2=[meshData.W2 mesh.W2(mesh.trunca+1:length(mesh.h)-mesh.truncb)];
    meshData.P0=[meshData.P0 mesh.P0(mesh.trunca+1:length(mesh.h)-mesh.truncb)];
    meshData.P1=[meshData.P1 mesh.P1(mesh.trunca+1:length(mesh.h)-mesh.truncb)];
    meshData.P2=[meshData.P2 mesh.P2(mesh.trunca+1:length(mesh.h)-mesh.truncb)];

    %truncb-part needs to refined and added
    x = xmin+sum(meshData.h);
    i=length(meshData.h);
    meshData.truncb=mesh.truncb;

    for g=length(mesh.h)-mesh.truncb+1:length(mesh.h)
        %split interval g into two intervals
         h=mesh.h(g)/2;
         xend=x+h;
         while x < xend
            hnew = 2*h;
            i=i+1;
            it=0;
            while abs(hnew/h-1) > 0.1 && it < 20
                h=max(min(hnew,xend-x),eps*10);
                [err,qq,ww,pp] = getErrorSLP(slp.q,slp.p,slp.w,h,x);
                err=h^2*err;        
                hnew = 5*h;
                if(err > eps) 
                    hnew = h * (tol / err)^(1/6);
                end
                if abs(h -min(xend-x,5e5))<eps && hnew>h, break; end
                it=it+1;
            end
            meshData.h(i)=h;
            meshData.Q0(i)=qq(1);meshData.Q1(i)=qq(2);meshData.Q2(i)=qq(3); %store data related to this mesh for later use
            meshData.W0(i)=ww(1);meshData.W1(i)=ww(2);meshData.W2(i)=ww(3);
            meshData.P0(i)=pp(1);meshData.P1(i)=pp(2);meshData.P2(i)=pp(3);
            x = x + h;
        end
        i=i+1; 
        h=mesh.h(g)/2; 
        [err,qq,ww,pp] = getErrorSLP(slp.q,slp.p,slp.w,h,x);
        meshData.h(i)=h;
        meshData.Q0(i)=qq(1);meshData.Q1(i)=qq(2);meshData.Q2(i)=qq(3); %store data related to this mesh for later use
        meshData.W0(i)=ww(1);meshData.W1(i)=ww(2);meshData.W2(i)=ww(3);
        meshData.P0(i)=pp(1);meshData.P1(i)=pp(2);meshData.P2(i)=pp(3);
    end
end

meshData.Infa=mesh.Infa;
meshData.Infb=mesh.Infb;
    



% 
% function [e_loc,qq,ww,pp] = getErrorSLP(qfun,pfun,wfun,h,x)
% [qq0,qq1,qq2]=getLegendrePolyn(qfun,x,h);
% [ww0,ww1,ww2]=getLegendrePolyn(wfun,x,h);
% [pp0,pp1,pp2]=getLegendrePolyn(@(x)1/pfun(x),x,h);
% qf=qfun(x+h); wf=wfun(x+h); pf=pfun(x+h);
% h2=h^2;
% qe=abs(qf-(qq0+qq1/h2+qq2/h2));%/max(1,abs(qf)); 
% we=abs(wf-(ww0+ww1/h2+ww2/h2));%/max(1,abs(wf));
% pe=abs(1/pf-(pp0+pp1/h2+pp2/h2));%/max(1,abs(1/pf));
% e_loc=max([qe we pe]);

% 
% function [V0,V1,V2]=getLegendrePolyn(fun,x,h)
% %Gauss Legendre (3 Legendre points)
% %realizes approximation by a series over shifted Legendre polynomials
% % the F_s coefficients are returned
% cc=sqrt(15)/5;
% c=h/2*[-cc+1 1 cc+1];
% f=arrayfun(fun,x+c);
% V0=1/18*(8*f(2)+5*f(1)+5*f(3));
% V1=(3*(-1/18*15^(1/2)*(f(1)-f(3))))*h^2;
% V2=5/9*(-2*f(2)+f(1)+f(3))*h^2;
