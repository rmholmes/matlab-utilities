function CoffeyEvans
%Coffey-Evans problem

%construct slp-object
o=slp(@p,@q,@w,-pi/2,pi/2,1,0,1,0);

%display information returned by the automatic classification algorithm:
disp(classificationInfo(o))

%compute first 11 eigenvalues:
[E,meshData] = computeEigenvalues(o,0,10,1e-12,true);

disp('The first 11 eigenvalues:')
line='-----------------------------------------';
disp(line)
disp('  k     E_k           estimated error')
disp(line)
for i=1:length(E.eigenvalues)
    fprintf('%4.0f  %16.12f  %+5.2e\n', E.indices(i),E.eigenvalues(i), E.errors(i))
end
disp(line)


%compute third eigenfunction
index=2;
e=E.eigenvalues(index+1); %value of the eigenvalue with index 2(as returned by computeEigenvalues)
figure
%evaluate eigenfunction in meshpoints:
[x1,y1,yp] = o.computeEigenfunction([index e],meshData);
%evaluate eigenfunction in a more dense set of x-values:
[x2,y2,yp] = o.computeEigenfunction([index e],meshData,linspace(-pi/2,pi/2,2^9));
plot(x1,y1,'k.',x2,y2,'k')
xlabel('x')
title('eigenfunction y_2')
% %compute eigenfunction with index 3:
index=3;
e=E.eigenvalues(index+1);
figure
[x1,y1,yp] = o.computeEigenfunction([index e],meshData);
[x2,y2,yp] = o.computeEigenfunction([index e],meshData,linspace(-pi/2,pi/2,2^9));
plot(x1,y1,'k.',x2,y2,'k')
xlabel('x')
title('eigenfunction y_{3}')
end

function r=p(x)
r=1;
end

function r=q(x)
b=30;
r=-2*b*cos(2*x)+b^2*sin(2*x)^2;
end

function r=w(x)
r=1;
end