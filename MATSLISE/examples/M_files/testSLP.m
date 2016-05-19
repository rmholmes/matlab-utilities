function testSLP
%test problem for the SLCPM12 package (Paine problem)
o=slp(@p,@q,@w,0,-sqrt(0.2)+sqrt(0.2+2*pi),1,0,1,0);

%display information returned by the automatic classification algorithm:
disp(classificationInfo(o))

%compute first eigenvalues:
[E,meshData] = o.computeEigenvalues(0,10,1e-13,true);
disp('Number of mesh intervals:')
numberOfMeshIntervals=length(meshData.h) %number of mesh intervals

%compute eigenfunction:
e=E.eigenvalues(1);
%in the meshpoints:
[x1,y,yp1,sq] = o.computeEigenfunction(e,meshData);
figure
plot(x1,y,'b*')
hold on
%in a more dense set of x-values:
[x,y,yp2,sq] = o.computeEigenfunction(e,meshData,0:0.01:-sqrt(0.2)+sqrt(0.2+2*pi));
plot(x,y,'g')
hold off
title('eigenfunction y_0')
figure
plot(x1,yp1,'r*:') %derivative of eigenfunction: py'
title('eigenfunction derivative py''')

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
r=(sqrt(0.2)+x)^3;
end

function r=q(x) 
r=4*(sqrt(0.2)+x);
end

function r=w(x)
r=(sqrt(0.2)+x)^5;
end