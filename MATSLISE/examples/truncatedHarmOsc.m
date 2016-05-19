function example1
%illustration of the use of the MATSLISE-functions
% Harmonic oscillator problem over the interval [-10,10]

%creation of the slp-object:
o=slp(@p,@q,@w,-10,10,1,0,1,0);

%display classification information:
disp(classificationInfo(o))

%compute eigenvalues with indices between 0 and 10:
[E,meshData] = computeEigenvalues(o,0,10,1e-14,true);

%display eigenvalue results:
disp('The first 11 eigenvalues of the (truncated) harmonic oscillator problem:')

line='-----------------------------------------';
disp(line)
disp('  k     E_k           estimated error')
disp(line)
for i=1:length(E.eigenvalues)
    fprintf('%4.0f  %16.12f  %+5.2e\n', E.indices(i),E.eigenvalues(i), E.errors(i))
end
disp(line)


%makes a plot of the mesh used to obtain the results:
plotMesh(o,meshData);
%note: on this symmetric problem, half range reduction was applied
%automatically

%eigenfunction computations:
figure
%only in meshpoints:
[x,y,yp] = computeEigenfunction(o,E.eigenvalues(9),meshData);
plot(x,y,'b')
hold on
%user-defined mesh is passed as last input-argument:
[x,y,yp] = computeEigenfunction(o,E.eigenvalues(9),meshData,-10:0.1:10);
plot(x,y,'r')
hold off
end




%coefficient functions
function r=p(~)
r=1;
end

function r=q(x)
r=x.^2;
end

function r=w(~)
r=1;
end