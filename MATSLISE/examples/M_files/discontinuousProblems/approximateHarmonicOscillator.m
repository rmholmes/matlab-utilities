function approximateHarmonicOscillator
%SLTSTPAK #53
%%q(x) = piecewise linear function joining the points
%(x, x^2) for x an integer

jumps=-10:10;  %vector with location of the jumps in the coefficient function q
o=slp(@p,@q,@w,-10,10,1,0,1,0,jumps); 
disp(classificationInfo(o))

%compute the first eigenvalues:
[E,meshData] = o.computeEigenvalues(0,10,1e-12,true);
%compute the first eigenfunction:
e=E.eigenvalues(1);
figure
%eigenfunction in the meshpoints:
[x,y,yp] = o.computeEigenfunction(e,meshData);
hold on
plot(x,y,'b.')
%eigenfunction on a more dense set of x-values, to produce smoother plot:
[x,y,yp] = o.computeEigenfunction(e,meshData,linspace(-10,10,512));
plot(x,y,'r')
xlabel('x')
hold off

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

function y=q(x)
%q(x) = piecewise linear function joining the points
%(x, x^2) for x an integer
x1= floor(x);%lower integer
x2= x1+1;%upper integer
y1=x1^2;
y2=x2^2;
y=(y2-y1)*(x-x1)+y1;
end

function r=w(x)
r=1;
end