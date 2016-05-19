
function pruf=locateEigenvalue(slp,eold,part,reference)
% LOCATE_EIGENVALUE
% pruf=locate_eigenvalue(slp,eold,part,reference)
%   
% For given e (eold) the Schrodinger equation is integrated from the
% both ends up to the matching point to produce the Prufer's 
% theta(eold)/pi (pruf), which gives us information about the index of eold.
% reference is a boolean value which indicates if the basic method
% (CPM{16,14}) or the higher order reference method (CPM{18,16}) needs to
% be used


if part.LNF
    %matching point is fixed
    imatch=part.imatch;
    Vmatch=part.V0(imatch);
    part.S = 1;
    if((eold-Vmatch)> 1) 
          part.S = sqrt(eold-Vmatch);
    end
    ev0=eold-part.V0(1);
    ev1=eold-part.V0(end);
else
    %matching point depends on the value of eold
    tmp=part.trunca+1:length(part.h)-part.truncb;
    [mm,imatch]=max((eold*part.W0(tmp)-part.Q0(tmp)).*part.P0(tmp));
    imatch=imatch+part.trunca;
    part.imatch=imatch;
    part.S=sqrt(max(mm,1));
    ev0=(eold*part.W0(1)-part.Q0(1))*part.P0(1);
    ev1=(eold*part.W0(end)-part.Q0(end))*part.P0(end);

end
a0 = slp.a0;
b0 = slp.b0;
a1 = slp.a1;
b1 = slp.b1;
sq1 = sqrt(a0^2+b0^2);
sq2 = sqrt(a1^2+b1^2);
if part.radial
    [yi,yei,thforw]=initstep(slp,eold,part);
else
    yi = [-b0 / sq1;a0 / sq1]; %initial values
    [thforw]= prufer(yi(1),yi(2),yi(1),yi(2),ev0,0,part.S);
end
yf = [-b1 / sq2;a1 / sq2];
if thforw < 0
     thforw = pi+thforw; % thforw = theta_min 
end

theta = propagateSolution1(eold,yi,part,true,reference,thforw);

thforw = thforw + theta; % thforw = thforw + delta_i;

%backward
[thbackw] = prufer(yf(1),yf(2),yf(1),yf(2),ev1,0,part.S);
if thbackw > 0
    thbackw = thbackw-pi;
end

if length(part.h) > imatch
      theta = propagateSolution1(eold,yf,part,false,reference,-thbackw);
      thbackw = thbackw - theta; % thbackw = thbackw - delta_i; 
end
pruf = (thforw-thbackw)/pi;




function theta = propagateSolution1(e,y,part,forward,reference,thinit)
imatch=part.imatch;
S=part.S;
theta=0;

if forward
 t = 1:imatch;
else
 t = length(part.h):-1:imatch+1;    
end

hm=part.h;


%wght1(q+1) = Coefficient of Z^q in series expansion for eta_m. (voor M=9)
%wght2(q+1) = Coefficient of Z^q in series expansion for eta_{m-1}.
wght1 =[1.527349308567059e-009; 0.36365459727787e-10; 0.00395276736172e-10; 0.00002635178241e-10; 0.00000012199899e-10; 0.00000000042069e-10;...
            0.00000000000113e-10; 0];
wght2 = [2.901963686277412e-008;  0.76367465428353e-9; 0.00909136493195e-9; 0.00006587945603e-9; 0.00000032939728e-9; 0.00000000121999e-9;...
            0.00000000000351e-9; 0.00000000000001e-9];
mmaxm=7;

if part.LNF

    if reference
        cum=part.rcu;
        cupm=part.rcup;
        cvm=part.rcv;
        cvpm=part.rcvp;
    else
        cum=part.cu;
        cupm=part.cup;
        cvm=part.cv;
        cvpm=part.cvp;
    end
    vbarm=part.V0;
    M = size(cum,2);
    eta=zeros(1,M+1);

    for tt = t
        y1=y;
        cu=cum(tt,:);
        cup=cupm(tt,:);
        cv=cvm(tt,:);
        cvp=cvpm(tt,:);
        vbar=vbarm(tt);
        h=hm(tt);
        Z=(vbar-e)*h^2;
        %calculation of eta functions:
        if abs(Z) < 0.5          
           tmp=(Z.^(0:mmaxm));
           eta(10) = tmp*wght1;
           eta(9) = tmp*wght2;
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
        tmp = eta(1:M); tmp=tmp(:); %eta(1:M)';
        eta1=eta(1);
        u  = csi + cu*tmp;
        up = (Z * eta1 + cup*tmp)/h;
        v = (eta1 + cv*tmp)*h;
        vp = csi + cvp*tmp;
        if forward
            T=[u, v ; up, vp];
        else
            T=[vp, -v ; -up, u];
        end
        if e<vbar-eps && Z<1e3
           sigma = exp(sqrt(vbar-e)*h); %cfr.SLEDGE
           y = T*y/sigma;
        else
           y = T*y;
        end
        if forward
              [theta0,theta1] = prufer(y1(1),y1(2),y(1),y(2),e-vbar,h,S);
        else
              [theta0,theta1] = prufer(y(1),y(2),y1(1),y1(2),e-vbar,h,S);
        end
        thetaprev=theta;
        theta=theta+theta1-theta0;
        if isnan(theta) || isinf(theta)
            error('Difficulties detected during propagation due to overflow')
        end
        if ((tt==part.trunca) || (tt==length(part.h) && part.truncb)) && theta-thetaprev>pi
            %singular problem (not radial)
            %initial step may be still too big, be carefull with
            %theta-results on this interval
            theta=thetaprev;
        end
        if floor((thinit+thetaprev)/pi)> floor((thinit+theta)/pi) && abs(thetaprev-theta)>1e-3 && thinit+thetaprev>1e-3% can never be decreasing through multiple of pi
                theta=theta+pi;
        end
        if y1(1)*y(1)<0 && theta<thetaprev% wavefunction goes through zero, then theta should increase through multiple of pi
                theta=theta+pi;
        end
    end
else %not in LNF
    M=4;
    for tt = t
            y1=y;
            h=hm(tt);
            P0=part.P0(tt); P1=part.P1(tt); P2=part.P2(tt);
            QW1=part.Q1(tt)-e*part.W1(tt); QW0=(part.Q0(tt)-e*part.W0(tt));
            Z=QW0*P0*h^2;
            R1=QW1; R0=QW0;
            U1=R1*P0-P1*R0;
            V1=R1*P0+P1*R0;
            wP1=P1/(h^2*P0);
            R2=part.Q2(tt)-e*part.W2(tt);
            if reference
                P2=0;
                R2=0;
            end
            U2=R2*P0-P2*R0;
            V2=R2*P0+P2*R0;
            wP2=P2/(h^2*P0);
            
            if abs(Z) < 0.5          
               tmp=(Z.^(0:mmaxm));
               eta(10) = tmp*wght1;
               eta(9) = tmp*wght2;
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
            eta0=eta(1);  eta1=eta(2); eta2=eta(3); eta3=eta(4); eta4=eta(5);



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

           
            if QW0*P0 > eps && Z <1e3
                   sigma = exp(sqrt(QW0*P0)*h); %cfr.SLEDGE
             else
                   sigma = 1;
            end

           if forward
             y=T*y/sigma; 
           else
             iT(1,1)=T(2,2); iT(1,2)=-T(1,2); iT(2,1)=-T(2,1); iT(2,2)=T(1,1);
             y=iT*y/sigma;  
           end        
            if forward
                [theta0,theta1] = prufer(y1(1),y1(2),y(1),y(2),-QW0*P0,h,S);
                 if tt==part.trunca                      
                     while theta1>theta0+pi/2
                         theta1=theta1-pi;
                     end
                 end
            else
                [theta0,theta1] = prufer(y(1),y(2),y1(1),y1(2),-QW0*P0,h,S);
                 if tt>length(part.h)-part.truncb
                    while theta1>theta0+pi/2
                         theta1=theta1-pi;
                    end
                end
            end
            thetaprev=theta;
            theta=theta+theta1-theta0;
            if (tt~=part.trunca) && y1(1)*y(1)<0 && (theta<thetaprev || (thetaprev==0 && abs(theta-thetaprev)<1e-10))  % wavefunction goes through zero, then theta should increase through multiple of pi
                %so theta cannot be decreasing
                 theta=theta+pi;
            end
            if floor((thinit+thetaprev)/pi)> floor((thinit+theta)/pi) && abs(thetaprev-theta)>1e-3 && thinit+thetaprev>1e-3 % can never be decreasing through multiple of pi
                theta=theta+pi;
            end
            if thetaprev-theta>2.5 %thetaprev is larger than theta and difference is almost pi
                %such a sudden drop of almost pi in theta is probably not correct, and explained
                %by a close call between pi/2 and -pi/2 somewhere in the
                %Prufer procedure
                %pryce31
                theta=theta+pi;
            end
    end
end



function [theta0,theta1]=prufer(y0,y0p,y1,y1p,ev,h,S)
if ev > 0 %(i): e > Vi
    Si = sqrt(ev); %local scale factor
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
    %end
    if theta1<theta0 && y0*y1<0
        theta1=theta1+pi; %het kan niet dat in bepaald interval er een zero in y zit maar dat theta toch daalt,
        %dan zou het door veelvoud van pi moeten gaan en het stijgt altijd
        %door veelvoud van pi
    end
    if Si == S
        return;
    end
%   else rescaling procedure;
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



