function x = exp(x)
e=exp(x.tc(1));
x.tc=[e,e*x.tc(2),x.tc(3)*e+x.tc(2)^2*e];
