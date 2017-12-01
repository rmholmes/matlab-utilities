function [cs,Ss,Cs] = Bmodes(z,N2,m)
%This script solves the Sturm-Louiville equation for the baroclinic
%modes given an input N2 on z grid points. the largest m wave speeds
%are returned in cs, with accompanying b/w normalized
%eigenfunctions Ss and u/v/p normalized eigenfunctions Cs.
%
% Ryan Holmes

g = 9.81;
N = length(z);

%Matrices:
D = diff(z);
Q = spdiags(2./(N2(2:(end-1)).*D(1:(end-1)).*D(2:end).*...
                (D(1:(end-1))+D(2:end))),0,N-2,N-2);

A = spdiags([-(D(1:(end-1))+D(2:end)) [0; D(1:(end-2))] [D(3:end); 0]],[0 1 ...
                    -1],N-2,N-2);

[EV,ED] = eigs(Q*A,m,'sm');
cs = sqrt(-1./diag(ED));

%%%Calculate vertical eigenfunctions:
Ss = [zeros(1,m); EV; zeros(1,m);];
%Normalize:
A = sqrt(g./sum(avg(Ss,1).*avg(Ss,1).*repmat(avg(N2).*D,[1 m]),1));
Ss = Ss.*repmat(A,[N 1]);

%conjugate:
Cs = repmat(cs',[N-1 1]).^2.*diff(Ss,[],1)./repmat(diff(z,[],1),[1 m])/g;
%Normalize:
A = sqrt(cs'.^2./sum(Cs.*Cs.*repmat(D,[1 m]),1)/g);
Cs = Cs.*repmat(A,[N-1 1]);

%Plotting:
if (0)
    h=subplot(1,3,1);
    plot(N2,z);
    ylabel('Depth (m)');
    xlabel('$N^2$');
    h1=subplot(1,3,2);
    plot(Cs,avg(z));
    xlabel('$C$');
    h2=subplot(1,3,3);
% $$$     plot(repmat(N2,[1 m]).*Ss,z);
    plot(Ss,z);
    xlabel('$S$');
    linkaxes([h h1 h2],'y');
end

end
