
%This function takes a central derivative of a field F defined at
%the locations x,y in the direction specified by type.
function dF = PerfDerv(x,y,F,type)

dF = NaN*zeros(size(F));

if (strcmp(type,'x'))
    dF(2:(end-1),:) = (F(3:end,:)-F(1:(end-1),:))./(x(3:end,:)-x(1:(end- ...
                                                      1),:));
elseif (strcmp(type,'y'))
    dF(:,2:(end-1)) = (F(:,3:end)-F(:,1:(end-1)))./(y(:,3:end)-y(:,1:(end- ...
                                                      1)));
end
end
