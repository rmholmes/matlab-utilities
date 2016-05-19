function x= log(x)
% LOG(X) for  AD object

    x.tc = [log(x.tc(1)),(1/(x.tc(1)))*x.tc(2),(-1/(x.tc(1))^2)*(x.tc(2))^2+(1/(x.tc(1)))*x.tc(3)];

end


