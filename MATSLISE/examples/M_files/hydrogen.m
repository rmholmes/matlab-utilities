function hydrogen
%hydrogen equation
o=slp(@p,@q,@w,0,inf,1,0,1,0);
%display information returned by the automatic classification algorithm:
disp(classificationInfo(o))

%compute first 11 eigenvalues:
[E,meshData] = computeEigenvalues(o,0,10,1e-10,true);
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

%computing a high eigenvalue:
E = computeEigenvalues(o,1000,1000,1e-10,true)
end

function r=p(x)
r=1;
end

function r=q(x) 
r=-1/x+2/x^2;
end

function r=w(x)
r=1;
end