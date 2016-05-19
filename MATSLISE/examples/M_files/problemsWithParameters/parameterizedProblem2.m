function parameterizedProblem2
%Illustration on how one can compare several problems: 
%e.g. study how small changes in the boundary conditions or in the potential function, influence the problem. 

%The problem studied here is the Coffey-Evans equation with a parameter in the potential function
%as in the test problem number 7 in the Pryce test set.

disp('Computing eigenvalues and eigenfunction....')
Eigenvalues=[];Eigenfunction=[];
for b=10:2.5:30
    %construction of slp-object for each parameter value b:
   o=slp(@p,@(x) q(b,x),@w,-pi/2,pi/2,1,0,1,0);
   %compute first 21 eigenvalues for each b-value:
   [E,meshData] = o.computeEigenvalues(0,20,1e-12,true);
   Eigenvalues=[Eigenvalues E.eigenvalues'];
   %compute eigenvalue with index 2 for each b-value:
   [x,y,yp] = o.computeEigenfunction([2 E.eigenvalues(3)],meshData,linspace(-pi/2,pi/2,100));
   Eigenfunction=[Eigenfunction y'];
end

figure
set(gcf,'DefaultAxesColorOrder',hsv(10))
plot(0:20,Eigenvalues,'LineStyle',':','Marker','.')
xlabel('k')
ylabel('E_k')
legend('b=10','b=12.5','b=15','b=17.5','b=20','b=22.5','b=25','b=27.5','b=30');
title('Eigenvalues')

figure
set(gcf,'DefaultAxesColorOrder',hsv(10))
plot(x,Eigenfunction)
xlabel('x')
ylabel('y_2')
legend('b=10','b=12.5','b=15','b=17.5','b=20','b=22.5','b=25','b=27.5','b=30');
axis tight
title('Eigenfunction')

end

function r=p(x)
r=1;
end

function r=q(b,x)
r=-2*b*cos(2*x)+b^2*sin(2*x)^2;
end

function r=w(x)
r=1;
end