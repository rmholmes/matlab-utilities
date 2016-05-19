function Klotter
% test problem 3 in SLTSTPAK
o=slp(@p,@q,@w,8/7,8,1,0,1,0);
%display information returned by the automatic classification algorithm:
disp(classificationInfo(o))

%compute first eigenvalues:
[E,meshData] = o.computeEigenvalues(0,9,1e-12,true);
disp('The first 11 eigenvalues:')
line='-----------------------------------------';
disp(line)
disp('  k     E_k           estimated error')
disp(line)
for i=1:length(E.eigenvalues)
    fprintf('%4.0f  %16.12f  %+5.2e\n', E.indices(i),E.eigenvalues(i), E.errors(i))
end
disp(line)

%eigenfunction computations:
e=E.eigenvalues(5);
figure
%only in meshpoints:
[x,y,yp] = o.computeEigenfunction(e,meshData);
hold on
plot(x,y,'b*')
%user-defined mesh is passed as last input-argument:
[x,y,yp] = o.computeEigenfunction(e,meshData,linspace(8/7,8,200));
plot(x,y,'r')
hold off
xlabel('x')
title('Eigenfunction')

end

function r=p(x)
r=1;
end

function r=q(x)
r=3/(4*x^2);
end

function r=w(x)
r=64*pi^2/(9*x^6);
end