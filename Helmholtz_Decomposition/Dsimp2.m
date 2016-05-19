
%This function takes a central derivative on a staggered grid of a
%field F defined on a regular grid.
function dF = Dsimp2(dx,F,type)

    if (strcmp(type,'x'))
        dF = (F(2:end,:)-F(1:(end-1),:))/dx;
    else
        dF = (F(:,2:end)-F(:,1:(end-1)))/dx;
    end
end
