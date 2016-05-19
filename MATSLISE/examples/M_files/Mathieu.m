function Mathieu
%illustration of the use of the MATSLISE-functions
%Mathieu equation

o=slp(@p,@q,@w,0,pi,1,0,1,0);

%display classification information:
disp(classificationInfo(o))

%compute eigenvalues with indices between 0 and 10:
[E,meshData] = o.computeEigenvalues(0,30,1e-12,true);

%first eigenfunction evaluated in the mesh points:
[x,y,yp] = o.computeEigenfunction(E.eigenvalues(1),meshData);
figure
plot(x,y,x,yp)
xlabel('x')
legend('y','yp')

%10th eigenfunction evaluated in a range of user-specified points:
[x,y,yp] = o.computeEigenfunction(E.eigenvalues(21),meshData,linspace(0,pi,500));
figure
%plot(x,y,x,yp)
plot(x,y)
legend('y_{20}')
xlabel('x')
axis tight
%legend('y','yp')


plotMesh(o,meshData);

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
r=2*cos(2*x);
end

function r=w(x)
r=1;
end