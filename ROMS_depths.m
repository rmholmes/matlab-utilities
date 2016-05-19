function [z_rho,z_w] = ROMS_depths(hname,zeta,h)
%This function is a quick ROMS depths calculating function where zeta
%and h are input, and the history filename is used to find the
%stretching functions/parameters
%
% Note: The last dimension of the outputs z_rho and z_w are the
% depth dimensions.

%Get parameters:
Vtransform = ncread(hname,'Vtransform');
Vstretching = ncread(hname,'Vstretching');
hc = ncread(hname,'hc');
Cs_r = ncread(hname,'Cs_r');
Cs_w = ncread(hname,'Cs_w');
s_rho = ncread(hname,'s_rho');
s_w   = ncread(hname,'s_w');

N = length(s_rho);
Nw = N+1;

z_rho = zeros([size(zeta) N]);
z_w   = zeros([size(zeta) Nw]);

if (length(size(z_rho)) == 2)
    str = '(:,k)';
elseif (length(size(z_rho)) == 3)
    str = '(:,:,k)';
elseif (length(size(z_rho)) == 4)
    str = '(:,:,:,k)';
end

if (Vtransform == 1)
    for k=1:N
        z0=(s_rho(k)-Cs_r(k))*hc + Cs_r(k).*h;
        eval(['z_rho' str ' = z0 + zeta.*(1.0 + z0./h);']);
    end    
        eval(['z_w' strrep(str,'k','1') '=-h;']);
    for k = 2:Nw
        z0=(s_w(k)-Cs_w(k))*hc + Cs_w(k).*h;
        eval(['z_w' str ' = z0 + zeta.*(1.0 + z0./h);']);
    end
else
    for k=1:N
        z0=(hc.*s_rho(k)+Cs_r(k).*h)./(hc+h);
        eval(['z_rho' str ' =zeta+(zeta+h).*z0;']);
    end
    for k=1:Nw
        z0=(hc.*s_w(k)+Cs_w(k).*h)./(hc+h);
        eval(['z_w' str '=zeta+(zeta+h).*z0;']);
    end
end

end


