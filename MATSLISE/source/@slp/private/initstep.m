function [varargout]=initstep(slpObject,e,part)
%specially tuned algorithm consistent with the singular nature around the
%origin
max_number_corrections=15;
Vs=part.Vs;
Vm1=Vs(1); V0=Vs(2); V1=Vs(3); V2=Vs(4);
l=part.l;
r=part.r0;
y0=r^(l+1);
y0p=(l+1)*r^l;
Vbar=V0-e;
Vs=[Vm1 Vbar V1 V2]';
Vsd=[0 -1 0 0]'; %derivative w.r.t. E
y=y0;
yp=y0p;
ye=0;
yep=0;
B=zeros(60,15);
Be=zeros(60,15);
for p=1:4
       B(p,1)=Vs(p)/(p*(2*l+p+1));
       Be(p,1)=Vsd(p)/(p*(2*l+p+1));
       y=y+B(p,1)*r^(l+p+1);
       ye=ye+Be(p,1)*r^(l+p+1);
       yp=yp+B(p,1)*(l+p+1)*r^(l+p);
       yep=yep+Be(p,1)*(l+p+1)*r^(l+p);
end

A=zeros(60,15);
Ae=zeros(60,15);
for q=2:max_number_corrections
    yprev=y;
    for p=q:4*q
        if p<4*q-2
            A(p,q)=Vm1*B(p-1,q-1);
            Ae(p,q)=Vm1*Be(p-1,q-1);
        end
        if p>q && p<4*q-1
            A(p,q)=A(p,q)+Vbar*B(p-2,q-1);
            Ae(p,q)=Ae(p,q)-B(p-2,q-1)+Vbar*Be(p-2,q-1);
        end
        if p>q+1 && p < 4*q
            A(p,q)=A(p,q)+V1*B(p-3,q-1);
            Ae(p,q)=Ae(p,q)+V1*Be(p-3,q-1);
        end
        if p>q+2
            A(p,q)=A(p,q)+V2*B(p-4,q-1);
            Ae(p,q)=Ae(p,q)+V2*Be(p-4,q-1);
        end
        B(p,q)=A(p,q)/(p*(1+p+2*l));
        Be(p,q)=Ae(p,q)/(p*(1+p+2*l));
        y=y+B(p,q)*r^(l+p+1);
        yp=yp+B(p,q)*(l+p+1)*r^(l+p);
        ye=ye+Be(p,q)*r^(l+p+1);
        yep=yep+Be(p,q)*(l+p+1)*r^(l+p);
    end
    if abs(y-yprev)/abs(y)<1e-15  %tkan zijn dat y zelf zeer klein is dus treshold moet misschien van y afhangen?
              break;
    end
end

varargout{1}=[y;yp];
varargout{2}=[ye;yep];


if nargout>2 %prufer phase also requested
    sn=part.S;
    step=r/10;
    r=step;
    
    parttemp.l=l;
    parttemp.Vs=part.Vs;
    parttemp.r0=r;
    yt=initstep(slpObject,e,parttemp); 
    y1=yt(1);
    
    nzeros=0;
    for i=2:10
        r=i*step;
        y_prev=y1;
        %y1=r^(l+1);
        parttemp.r0=r;
        y1=initstepFast(e,parttemp);
        if y1*y_prev<=0
            nzeros=nzeros+1;
        end
    end
    % Generation of the Prufer phase
    tmp = atan(sn*y/yp);
    if (-1)^nzeros*yp < 0
      tmp = pi + tmp;
    end
    varargout{3} = nzeros*pi + tmp;  
end
end


function y=initstepFast(e,part)
max_number_corrections=10;
Vs=part.Vs;
Vm1=Vs(1); V0=Vs(2); V1=Vs(3); V2=Vs(4);
l=part.l;
r=part.r0;
y0=r^(l+1);
Vbar=V0-e;
Vs=[Vm1 Vbar V1 V2]';
y=y0;
B=zeros(40,10);
for p=1:4
       B(p,1)=Vs(p)/(p*(2*l+p+1));
       y=y+B(p,1)*r^(l+p+1);
end

A=zeros(40,10);
for q=2:max_number_corrections
    yprev=y;
    for p=q:4*q
        if p<4*q-2
            A(p,q)=Vm1*B(p-1,q-1);
        end
        if p>q && p<4*q-1
            A(p,q)=A(p,q)+Vbar*B(p-2,q-1);
        end
        if p>q+1 && p < 4*q
            A(p,q)=A(p,q)+V1*B(p-3,q-1);
        end
        if p>q+2
            A(p,q)=A(p,q)+V2*B(p-4,q-1);
        end
        B(p,q)=A(p,q)/(p*(1+p+2*l));
        y=y+B(p,q)*r^(l+p+1);
    end
    if abs(y-yprev)/abs(y)<1e-13  %tkan zijn dat y zelf zeer klein is dus treshold moet misschien van y afhangen?
              break;
    end
end

end
% 
% function [varargout]=initstep(slpObject,e,part)
% %specially tuned algorithm consistent with the singular nature around the
% %origin
% max_number_corrections=15;
% Vs=part.Vs;
% Vm1=Vs(1); V0=Vs(2); V1=Vs(3); V2=Vs(4);
% l=part.l;
% r=part.r0;
% y0=r^(l+1);
% y0p=(l+1)*r^l;
% Vbar=V0-e;
% i=1;
% Vs=[Vm1 Vbar V1 V2]';
% Vsd=[0 -1 0 0]'; %derivative w.r.t. E
% y=y0;
% ye=0;
% maxn=max_number_corrections*4;
% tmpB=zeros(1,maxn);
% tmpBe=zeros(1,maxn);
% for p=i:i*4
%        tmp=zeros(1,4);
%        tmp(p)=1;
%        tmpB(p)=(tmp*Vs)/((l+p)*(l+p+1)-l*(l+1));
%        tmpBe(p)=(tmp*Vsd)/((l+p)*(l+p+1)-l*(l+1));
%        y=y+tmpB(p)*r^(l+p+1);
%        ye=ye+tmpBe(p)*r^(l+p+1);
%        %yep=yep+(l+p+1)*Be(i,p)*r^(l+p);
% end
% y1p=Vm1*r^(l+2)*(l+2)/r/(2*l+2)+Vbar*r^(l+3)*(l+3)/r/(6+4*l)+V1*r^(l+4)*(l+4)/r/(6*l+12)+V2*r^(l+5)*(l+5)/r/(8*l+20);
% yp=y0p+y1p;
% yep=-r^(l+3)*(l+3)/r/(4*l+6);
% 
% Vmatrix=diag(Vm1*ones(4,1));
% Vmatrix(1,2:4)=[Vbar V1 V2];
% Vmatrix(2,3:4)=[Vbar V1];
% Vmatrix(3,4)=Vbar;
% Vmatrixe=[0 -1 0 0; 0 0 -1 0; 0 0 0 -1; 0 0 0 0];
% ll=l*(l+1);
% tmpA=(1+l+(1:maxn))';
% rrs=(r.^tmpA);
% rrsd=(tmpA.*rrs)/r;
% 
% for i=2:max_number_corrections
%           Bs=[tmpB(i-1:(i-1)*4) 0 0 0 0];
%           Bse=[tmpBe(i-1:(i-1)*4) 0 0 0 0];
%           tmpB=zeros(1,maxn);
%           tmpBe=zeros(1,maxn);
%           %voor p=i
%           tmpA = Vm1*Bs(1);
%           tmpAe = Vm1*Bse(1);
%           lp=l+i;
%           tn=(lp*(lp+1)-ll);
%           tmpB(i)=tmpA/tn;
%           tmpBe(i)=tmpAe/tn;
%           for p=i+1:i*4
%             lp=lp+1;  
%             tmpi=max(1,p-i-2):p-i+1;  
%             tmp = [Bs(tmpi) 0 0 0 0];
%             tmpV = Vmatrix(:,min(p-i+1,4));
%             tmpVe = Vmatrixe(:,min(p-i+1,4));
%             tmpA = tmp(1:4)*tmpV;
%             tmp2= [Bse(tmpi) 0 0 0 0];
%             %tmpAe = tmp2(1:4)*tmpV-tmp(2);
%             tmpAe = tmp2(1:4)*tmpV+tmp(1:4)*tmpVe;
%             tn=(lp)*(lp+1)-ll;
%             tmpB(p) = tmpA/tn;
%             tmpBe(p) = tmpAe/tn;
%           end
%           ytmp=tmpB*rrs;
%           y=y+ytmp;
%           yp=yp+tmpB*rrsd;
%           ye=ye+tmpBe*rrs;
%           yep=yep+tmpBe*rrsd;
%           if abs(ytmp)/abs(y)<1e-15  %tkan zijn dat y zelf zeer klein is dus treshold moet misschien van y afhangen?
%               break;
%           end
% end
% varargout{1}=[y;yp];
% varargout{2}=[ye;yep];
% 
% if nargout>2
%     imatch=part.imatch;
%     Vmatch=part.V0(imatch);
%     if e-Vmatch > 1
%       sn = sqrt(e-Vmatch);
%     else
%       sn = 1;
%     end
%     step=r/10;
%     r=step;
%     y1=r^(l+1);
%     nzeros=0;
%     for i=2:10
%         r=i*step;
%         y_prev=y1;
%         y1=r^(l+1);
%         if y1*y_prev<=0
%             nzeros=nzeros+1;
%         end
%     end
%     % Generation of the Prufer phase
%     tmp = atan(sn*y/yp);
%     if (-1)^nzeros*yp < 0
%       tmp = pi + tmp;
%     end
%     varargout{3} = nzeros*pi + tmp;  
% end