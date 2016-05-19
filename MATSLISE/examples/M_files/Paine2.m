function Paine2
%illustration of the use of the MATSLISE-functions
%Paine problem 2

%creation of the slp-object:
o=slp(@p,@q,@w,0,pi,1,0,1,0);

%display classification information:
disp(classificationInfo(o))

%compute eigenvalues with indices between 0 and 10:
[E,meshData] = o.computeEigenvalues(0,10,1e-12,true);

%display eigenvalue results:
disp('The first 11 eigenvalues:')
line='-----------------------------------------';
disp(line)
disp('  k     E_k           estimated error')
disp(line)
for i=1:length(E.eigenvalues)
    fprintf('%4.0f  %16.12f  %+5.2e\n', E.indices(i),E.eigenvalues(i), E.errors(i))
end
disp(line)

%compute eigenvalues between 100 and 200:
[E2,meshData] = o.computeEigenvalues(100,200,1e-12,false);


%display eigenvalue results:
disp('Eigenvalues between 100 and 200:')
disp(line)
disp('  k     E_k           estimated error')
disp(line)
for i=1:length(E2.eigenvalues)
    fprintf('%4.0f  %16.12f  %+5.2e\n', E2.indices(i),E2.eigenvalues(i), E2.errors(i))
end
disp(line)

%eigenfunction computations:
e=E.eigenvalues(1);
figure
%only in meshpoints:
[x,y,yp] = o.computeEigenfunction(e,meshData);
hold on
plot(x,y,'b*')
%user-defined mesh is passed as last input-argument:
[x,y,yp] = o.computeEigenfunction(e,meshData,linspace(0,pi,50));
plot(x,y,'r')
hold off

%makes a plot of the mesh used to obtain the results:
plotMesh(o,meshData);


end

function r=p(x)
r=1;
end

function r=q(x)
r=1/(x+0.1)^2;
end

function r=w(x)
r=1;
end