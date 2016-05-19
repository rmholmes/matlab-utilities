function Airy
%Airy equation (problem 27 in test set in book Pryce)

%construct slp-object:
o=slp(@p,@q,@w,0,inf,1,0,1,0);
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
r=x;
end

function r=w(x)
r=1;
end