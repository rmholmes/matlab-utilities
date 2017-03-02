function [anom,clima] = climatology(dnum,field,dim)
% This function subtracts monthly climatological averages from a
% given field and returns the anomalies and climatology.
    
    sz = size(field);
    szCLI = sz;
    szCLI(dim) = 12;
    clima = zeros(szCLI);
    anom = zeros(sz);
    
    dvec = datevec(dnum);
    
    if (dim == 1)
        for i=1:12
            inds = find(dvec(:,2) == i);
            clima(i,:,:,:) = mean(field(inds,:,:,:),1);
        end
        for i = 1:length(dnum)
            month = dvec(i,2);
            anom(i,:,:,:) = field(i,:,:,:) - clima(month,:,:,:);
        end
    
    elseif (dim == 2)
        for i=1:12
            inds = find(dvec(:,2) == i);
            clima(:,i,:,:) = mean(field(:,inds,:,:),2);
        end
        for i = 1:length(dnum)
            month = dvec(i,2);
            anom(:,i,:,:) = field(:,i,:,:) - clima(:,month,:,:);
        end

    elseif (dim == 3)
        for i=1:12
            inds = find(dvec(:,2) == i);
            clima(:,:,i,:) = mean(field(:,:,inds,:),3);
        end
        for i = 1:length(dnum)
            month = dvec(i,2);
            anom(:,:,i,:) = field(:,:,i,:) - clima(:,:,month,:);
        end

    elseif (dim == 4)
        for i=1:12
            inds = find(dvec(:,2) == i);
            clima(:,:,:,i) = mean(field(:,:,:,inds),4);
        end
        for i = 1:length(dnum)
            month = dvec(i,2);
            anom(:,:,:,i) = field(:,:,:,i) - clima(:,:,:,month);
        end
    end
end

