
%This function takes a central derivative of a field F defined on a
%regular grid.
function dF = Dsimp(dx,F,type)

    dF = NaN*zeros(size(F));
    if (strcmp(type,'x'))
        dF(2:(end-1),2:(end-1)) = (F(3:end,2:(end-1))-F(1:(end-2),2:(end-1)))/2/ ...
            dx;
    else
        dF(2:(end-1),2:(end-1)) = (F(2:(end-1),3:end)-F(2:(end-1),1:(end-2)))/2/ ...
            dx;
    end
end
