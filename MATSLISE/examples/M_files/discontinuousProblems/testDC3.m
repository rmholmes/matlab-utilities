function testDC3
%problem with 2 jumps
jumps=[-1 1];
o=slp(@p,@q,@w,-30,30,1,0,1,0,jumps);
disp(classificationInfo(o))

[E,meshData] = o.computeEigenvalues(0,10,1e-12,true);
e=E.eigenvalues(1);
figure
[x,y,yp] = o.computeEigenfunction(e,meshData);
hold on
plot(x,y,'b*')
[x,y,yp] = o.computeEigenfunction(e,meshData,-30:0.1:30);
plot(x,y,'r')
hold off

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
%a problem with two jumps (one in -1 and one in 1)
% definition of the three potential parts:
if abs(x)>1
    r=0;
else
    r=-1;
end
end

function r=w(x)
r=1;
end