function Legendre
%Legendre equation

%o=slp(@p,@q,@w,-1,1,1,0,1,0);

%display information returned by the automatic classification algorithm:
%disp(classificationInfo(o))

%call constructor of the class slp:
slpObject=slp(@p,@q,@w,-1,1,1,0,1,0);

%compute first eigenvalues with indices from 0 to 10
tic
E = computeEigenvalues(slpObject,0,10,1e-10,true);
toc

%display results:
disp('The first 11 eigenvalues:')
line='-----------------------------------------';
disp(line)
disp('  k     E_k             error')
disp(line)
for i=1:length(E.eigenvalues)
    exactError= E.eigenvalues(i)-(i*(i-1)+1/4);
    fprintf('%4.0f  %16.12f  %+5.2e\n', E.indices(i),E.eigenvalues(i),exactError)
end
disp(line)

%call constructor of the class slp:
slpObject=slp(@p,@q,@w,-1,1,1,0,1,0);

%compute some high eigenvalues
E = computeEigenvalues(slpObject,100,110,1e-14,true);

%display results:
disp('Some high eigenvalues:')
line='-----------------------------------------';
disp(line)
disp('  k     E_k             error')
disp(line)
for i=1:length(E.eigenvalues)
    n=E.indices(i);
    exactError= E.eigenvalues(i)-(n*(n+1)+1/4);
    fprintf('%4.0f  %16.10f  %+5.2e\n', E.indices(i),E.eigenvalues(i),exactError)
end
disp(line)
end

function r=p(x)
r=1-x^2;
end

function r=q(x) 
r=1/4;
end

function r=w(x)
r=1;
end


