function harmonic_oscillator
%harmonic oscillator
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


%compute first eigenfunction
e=E.eigenvalues(1); %value of the first eigenvalue (as returned by computeEigenvalues)
figure
%evaluate eigenfunction in meshpoints:
[x1,y1,yp] = o.computeEigenfunction(e,meshData);
%evaluate eigenfunction in a more dense set of x-values:
[x2,y2,yp] = o.computeEigenfunction(e,meshData,linspace(x1(1),x1(end),2^9));
plot(x1,y1,'k.',x2,y2,'k')
xlabel('x')
title('eigenfunction y_0')

%compute eigenfunction with index 10
e=E.eigenvalues(11); %value of the first eigenvalue (as returned by computeEigenvalues)
figure
%evaluate eigenfunction in meshpoints:
[x1,y1,yp] = o.computeEigenfunction(e,meshData);
%evaluate eigenfunction in a more dense set of x-values:
[x2,y2,yp] = o.computeEigenfunction(e,meshData,linspace(x1(1),x1(end),2^9));
plot(x1,y1,'k.',x2,y2,'k')
xlabel('x')
title('eigenfunction y_{10}')

end

function r=p(x)
r=1;
end

function r=q(x) 
r=x^2;
end

function r=w(x)
r=1;
end