function res = legendre(gammam,p)
% LEGENDRE  res = legendre(gamma,p)
% For input gamma (between 0 and 1) this returns the value of the shifted
% Legendre polynomials L_i(gamma), i=0..p
for j=1:length(gammam)
    gamma=gammam(j);
    res(j,1) = 1;
    res(j,2) = -1 + 2*gamma;
    x = 2* gamma -1;
    for i = 2:p
       res(j,i+1) = ((2*(i-1)+1)*x*res(j,i)-(i-1)*res(j,i-1))/(i);
    end
end