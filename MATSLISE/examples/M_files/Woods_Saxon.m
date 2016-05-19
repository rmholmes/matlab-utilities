function Woods_Saxon
%Woods-Saxon potential

%construct slp-object:
o=slp(@p,@q,@w,0,inf,1,0,1,0);

%display classification information:
disp(classificationInfo(o))

%compute first eigenvalues:
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


%first eigenfunction evaluated in the mesh points:
[x,y,yp] = o.computeEigenfunction(E.eigenvalues(1),meshData);
figure
plot(x,y,x,yp)
xlabel('x')
legend('y','yp')

%10th eigenfunction evaluated in a range of user-specified points:
[x,y,yp] = o.computeEigenfunction(E.eigenvalues(10),meshData,linspace(x(1),x(end),500));
figure
plot(x,y,x,yp)
xlabel('x')
legend('y','yp')
end

function r=p(x)
r=1;
end

function r=q(x) 
r=-50*(1-5*(exp((x-7)/0.6))/(3*(1+(exp((x-7)/0.6)))))/(1+(exp((x-7)/0.6)));
end

function r=w(x)
r=1;
end