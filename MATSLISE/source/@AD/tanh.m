function x = tanh(x)
% Tanh for AD objects.
x1=tanh(x.tc(1));
 x.tc = [x1,(1-x1^2)*x.tc(2),-2*x1*(1-x1^2)*x.tc(3)+(1-x1^2)*x.tc(3)];
