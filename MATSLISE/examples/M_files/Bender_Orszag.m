function Bender_Orszag
%a problem with a continuous spectrum and a limited number of discrete
%eigenvalues

%define slp-object:
o=slp(@p,@q,@w,-inf,inf,1,0,1,0);
%display information returned by the automatic classification algorithm:
disp(classificationInfo(o))

%compute first 11 eigenvalues:
[E,meshData] = computeEigenvalues(o,0,10,1e-12,true);
%display results:
disp('The first 11 eigenvalues:')
line='-----------------------------------------';
disp(line)
disp('  k     E_k           estimated error')
disp(line)
for i=1:length(E.eigenvalues)
    fprintf('%4.0f  %16.12f  %+5.2e\n', E.indices(i),E.eigenvalues(i), E.errors(i))
end
disp(line)
end

function r=p(x)
r=1;
end

function r=q(x)
m=10; 
r=-m*(m+1)/cosh(x)^2;
end

function r=w(x)
r=1;
end