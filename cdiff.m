function dF = cdiff(L,F,TYPE);
%This function calculates the central difference derivative of field F
%with respect to variable L in direction TYPE with locations L. L
%should be 2D (L=x or L=y) if x and y or L should be 1D (L=z) matching
%the respective dimension size of F.
%type = 'x'; take x derivative F 2D or 3D.
%type = 'y'; take y derivative F 2D or 3D.
%Type = 'z'  ; take z derivative with F being 3D.

switch TYPE
  case 'x'
    NFD = length(size(F));
    NLD = length(size(L));
    if (NFD == 2)
        [xL,yL] = size(F);
        cdm = cdiffm(xL);
        if (NLD == 2)
            dF = (cdm*F)./(cdm*L);
        else
            dF = (cdm*F)./(cdm*L(:,:,1));
        end
    else
        [xL,yL,zL] = size(F);
        cdm = cdiffm(xL);
        dF = zeros(size(F));
        if (NLD == 3)
            for z = 1:zL
                dF(:,:,z) = (cdm*F(:,:,z))./(cdm*L(:,:,z));
            end 
        else
            for z = 1:zL
                dF(:,:,z) = (cdm*F(:,:,z))./(cdm*L);
            end 
        end
    end
    
  case 'y'
    NFD = length(size(F));
    NLD = length(size(L));
    if (NFD == 2)
        [xL,yL] = size(F);
        cdm = cdiffm(yL);
        if (NLD == 2)
            dF = ((cdm*F')./(cdm*L'))';
        else
            dF = ((cdm*F')./(cdm*L(:,:,1)'))';
        end
    else
        [xL,yL,zL] = size(F);
        cdm = cdiffm(yL);
        dF = zeros(size(F));
        if (NLD == 3)
            for z = 1:zL
                dF(:,:,z) = ((cdm*F(:,:,z)')./(cdm*L(:,:,z)'))';
            end 
        else
            for z = 1:zL
                dF(:,:,z) = ((cdm*F(:,:,z)')./(cdm*L'))';
            end 
        end
    end
    
  case 'z'
    if (length(size(F))==3)
    [xL,yL,zL] = size(F);
    [B,C,Z] = meshgrid(1:yL,1:xL,L);
    dF = zeros(size(F));
    cdm = cdiffm(zL);
    for j=1:yL
        dF(:,j,:) = permute((cdm*permute(F(:,j,:),[3 1 2]))./ ...
                            (cdm*permute(Z(:,j,:),[3 1 2])),[2 3 1]);
    end
    else
        [xL,yL,zL,tL] = size(F);
        [B,C,Z] = meshgrid(1:yL,1:xL,L);
        dF = zeros(size(F));
        cdm = cdiffm(zL);
        for t=1:tL
            for j=1:yL
            dF(:,j,:,t) = permute((cdm*permute(F(:,j,:,t),[3 1 2]))./ ...
                                (cdm*permute(Z(:,j,:),[3 1 2])),[2 3 1]);
            end
        end 
    end
end
end

function cdm = cdiffm(n)
%This function creates the central difference matrix for
%calculating first-order derivatives on an array:
D = ones(n,1);
cdm = spdiags([-0.5.*D 0.5.*D], [-1 1],n,n);
cdm(1,1)=-1;
cdm(1,2)=1;
cdm(n,n)=1;
cdm(n,n-1)=-1;
end
