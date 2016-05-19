function meshData = refineMesh(slp,mesh)
%diffucties were detected during the eigenvalue computation, the mesh is
%refined by doubling the number of intervals;

%warning('refineMesh')

meshData.LNF=mesh.LNF;
meshData.halfRangeReduction=mesh.halfRangeReduction;
meshData.radial=mesh.radial;
meshData.Infa=mesh.Infa;
meshData.Infb=mesh.Infb;
if mesh.radial
    meshData.l=mesh.l;meshData.Vs=mesh.Vs; meshData.r0=mesh.r0; 
end
if meshData.LNF

    pot = slp.q; %the potential
    xmin = mesh.Infa;

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
    i=1;
    if mesh.radial
         meshData.r0=mesh.r0/2;
         meshData.Infa=meshData.r0;
         x=meshData.r0;
         h=mesh.r0/2;
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
         i=i+1;
    end
    for k=1:length(mesh.h);
        h=mesh.h(k)/2;
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
        i=i+1;
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
        i=i+1;
    end
    meshData.cu=cu';  %omdat kolommen toevoegen sneller gaat dan rijen toevoegen
    meshData.cup=cup'; meshData.cv=cv'; meshData.cvp=cvp';
    meshData.rcu=rcu';  meshData.rcup=rcup'; meshData.rcv=rcv'; meshData.rcvp=rcvp';    
    meshData.imatch=mesh.imatch*2;
    meshData.truncb=mesh.truncb;
   
else
    
    xmin = mesh.Infa;
   
    %refine each interval of mesh :
    x = xmin;
    meshData.h=[];
    meshData.Q0=[];    meshData.Q1=[];    meshData.Q2=[];
    meshData.W0=[];    meshData.W1=[];    meshData.W2=[];
    meshData.P0=[];    meshData.P1=[];    meshData.P2=[];
    meshData.trunca=mesh.trunca;
    i=1;
    if mesh.trunca
        h=mesh.h(1);
        [~,qq,ww,pp] = getErrorSLP(slp.q,slp.p,slp.w,h,x);
        meshData.h(i)=h; 
        meshData.Q0(i)=qq(1);meshData.Q1(i)=qq(2);meshData.Q2(i)=qq(3); %store data related to this mesh for later use
        meshData.W0(i)=ww(1);meshData.W1(i)=ww(2);meshData.W2(i)=ww(3);
        meshData.P0(i)=pp(1);meshData.P1(i)=pp(2);meshData.P2(i)=pp(3);
        x=x+h; i=i+1; 
    end
    for k=mesh.trunca+1:length(mesh.h)
        h=mesh.h(k)/2;
        [~,qq,ww,pp] = getErrorSLP(slp.q,slp.p,slp.w,h,x);
        meshData.h(i)=h; 
        meshData.Q0(i)=qq(1);meshData.Q1(i)=qq(2);meshData.Q2(i)=qq(3); %store data related to this mesh for later use
        meshData.W0(i)=ww(1);meshData.W1(i)=ww(2);meshData.W2(i)=ww(3);
        meshData.P0(i)=pp(1);meshData.P1(i)=pp(2);meshData.P2(i)=pp(3);
        x=x+h; i=i+1;
        [~,qq,ww,pp] = getErrorSLP(slp.q,slp.p,slp.w,h,x);
        meshData.h(i)=h; 
        meshData.Q0(i)=qq(1);meshData.Q1(i)=qq(2);meshData.Q2(i)=qq(3); %store data related to this mesh for later use
        meshData.W0(i)=ww(1);meshData.W1(i)=ww(2);meshData.W2(i)=ww(3);
        meshData.P0(i)=pp(1);meshData.P1(i)=pp(2);meshData.P2(i)=pp(3);
        i=i+1;
        x=x+h;
    end
    meshData.truncb=mesh.truncb;
end


    



