
function [y,ye] = propagateSolution(e,y,ye,part,forward,reference)
%propagation of the solution from left to right if forward = true, and from
%right to left if forward =false
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


hm=part.h;
%wght1(q+1) = Coefficient of Z^q in series expansion for eta_m.
%wght2(q+1) = Coefficient of Z^q in series expansion for eta_{m-1}.
wght1 =[1.527349308567059e-009 0.36365459727787e-10 0.00395276736172e-10 0.00002635178241e-10 0.00000012199899e-10 0.00000000042069e-10...
            0.00000000000113e-10 0];
wght2 = [2.901963686277412e-008  0.76367465428353e-9 0.00909136493195e-9 0.00006587945603e-9 0.00000032939728e-9 0.00000000121999e-9...
            0.00000000000351e-9 0.00000000000001e-9];
mmaxm=7;
if part.LNF %Liouville normal form -> CPM{16,14} algorithm is used
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
           dsigma = 1/2*(h/(sigma*sqrt(vbar-e)));
        else
           sigma = 1;
           dsigma=0;
        end
        ye = Te*y/sigma+T*(ye/sigma+y*dsigma);
        y = T*y/sigma;
    end
else %problem in full SLP-form -> sixth order CP method for SLPs is used
    M=5;
    for i = 1:length(t)
            tt=t(i);
            h=hm(tt);
            P0=part.P0(tt); P1=part.P1(tt); P2=part.P2(tt);
            W0=part.W0(tt);
            QW1=part.Q1(tt)-e*part.W1(tt); QW0=(part.Q0(tt)-e*part.W0(tt));
            R2=part.Q2(tt)-e*part.W2(tt);
            R1=QW1; R0=QW0;
            dU2=-part.W2(tt)*P0+P2*W0;
            dV2=-part.W2(tt)*P0-P2*W0;
            if reference
                %lower order method (approximation of coeff.functions by first degree
                %polynomial) is used for reference
                P2=0; 
                R2=0;
                dU2=0;
                dV2=0;
            end
            U1=R1*P0-P1*R0;
            V1=R1*P0+P1*R0;
            wP1=P1/(h^2*P0);            
            U2=R2*P0-P2*R0;
            V2=R2*P0+P2*R0;
            wP2=P2/(h^2*P0);
            dU1= -part.W1(tt)*P0+P1*part.W0(tt); %derivative of U1 wrt E
            dV1= -part.W1(tt)*P0-P1*W0;
                        
            
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
            while csi > 1e100
                csi = csi/1e140;
                eta = eta./1e140;
            end
            eta0=eta(1);  eta1=eta(2); eta2=eta(3); eta3=eta(4); eta4=eta(5); eta5=eta(6);

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

            Te=-Te*W0*P0*h^2/2;


            Te(1,1)=Te(1,1)+((dU1*wP1)/12+(dU2*wP2)/20)*eta0+(-(dU1)/2-(wP1*(dV2+dU1)+wP2*dV1)/4)*eta1+...
                ((3*wP1*dV2+wP2*(7*dV1-3*dU2))/4-(dU2*V2)/40-(U2*dV2)/40-(dU1*V1)/24-(U1*dV1)/24)*eta2+...
                (5*dV1*(U2+2*V2)+5*V1*(dU2+2*dV2)+dV2*(5*U1-3*U2)+V2*(5*dU1-3*dU2))/40*eta3;

            Te(1,2)=Te(1,2)+h*P0*(+((wP1*(dU1+dV1))/12+(wP2*(dU2+dV2))/20)*eta1-...
                ((dV2)/2+(wP1*(dU1+dV1))/6+(wP2*(2*dU2+7*dV2))/10)*eta2+...
                (-(2*V2*dV2)/40-(2*V1*dV1)/24+(wP2*(9*dU2+24*dV2))/10)*eta3+(9*2*V2*dV2)/40*eta4);

            Te(2,1)=Te(2,1)+(((wP1*dU1)/12+(wP2*dU2)/20)*csi+((dU2)/2+(wP2*(dV2-dU2+dU1))/4+(wP1*(dV1-dU2))/4)*eta0+...
                (-(3*dU2)/2-(dU2*V2)/40-(U2*dV2)/40-(dU1*V1)/24-(U1*dV1)/24-(wP2*(9*dV2-3*dU2+3*dU1))/4-(wP1*(3*dV1+dU1-3*dU2))/4)*eta1+...
                (-(dU1*V1)/24-(U1*dV1)/24+(dU2*V2)/10+(U2*dV2)/10-(dU1*(V1+V2)+U1*(dV1+dV2)+dU2*(V2-V1)+U2*(dV2-dV1)+...
                2*V1*dV1+2*V2*dV2)/8-(3*wP2*(dU2-10*dV2))/4)*eta2+((dV2*(30*V2+27*U2)+V2*(30*dV2+27*dU2))/40)*eta3)/h/P0;

            Te(2,2)=Te(2,2)+(((dV1+dU1)*wP1)/12+((dV2+dU2)*wP2)/20)*eta0+...
               (((1-wP1)*dV1)/2+(wP1*dV2+wP2*dU1)/2)*eta1-((3*wP2*(dV2+2*dV1))/2+(wP1*(dV2+3*dU2))/2+...
               (2*V2*dV2)/40+(2*V1*dV1)/24)*eta2-(dV2*(3*V2+20*V1)+V2*(3*dV2+20*dV1))/40*eta3;
                    
                       
            if QW0*P0 > eps
                   sigma = exp(sqrt(QW0*P0)*h); %cfr.SLEDGE, to avoid overlow
                   dsigma = 1/2*(h*W0*P0/(sigma*sqrt(QW0*P0))); %derivative of 1/sigma
                   %dsigma moet in rekening gebracht worden, zie pryce31
                   %over [0.0001,100]
             else
                   sigma = 1;
                   dsigma = 0;
            end
            
            
            if forward
             ye = Te*y/sigma+T*(ye/sigma+y*dsigma); %ook afgeleide van sigma gebruiken?
             y = T*y/sigma;
           else
             iT(1,1)=T(2,2); iT(1,2)=-T(1,2); iT(2,1)=-T(2,1); iT(2,2)=T(1,1);
             iTe(1,1)=Te(2,2); iTe(1,2)=-Te(1,2); iTe(2,1)=-Te(2,1); iTe(2,2)=Te(1,1);
             ye = iTe*y/sigma+iT*(ye/sigma+y*dsigma);
             y=iT*y/sigma;
            end
    end
end