function parameterizedProblem
%Illustration on how one can compare several problems: 
%e.g. study how small changes in the boundary conditions or in the potential function, influence the problem. 

%The problem studied here is the Mathieu equation with a parameter c in the potential function
%as in the test problem number 5 in the Pryce test set.

%we look how the first eigenvalues change with c
%and then how the first eigenfunction varies with c
disp('Computing eigenvalues and eigenfunctions....')
Eigenvalues=[];Eigenfunction=[];
for c=1:5  %the parameter c takes the values 0,1,2,3,4,5
    %construct slp-object for each c-value:
   o=slp(@p,@(x) q(c,x),@w,0,40,1,0,1,0);
   %compute first 21 eigenvalues for each c-value:
   [E,meshData] = o.computeEigenvalues(0,20,1e-12,true);  
   Eigenvalues=[Eigenvalues E.eigenvalues'];
   %compute first eigenfunction for each c-value
   [x,y,yp] = o.computeEigenfunction(E.eigenvalues(1),meshData,0:0.1:40);
   Eigenfunction=[Eigenfunction y'];
end

%make plots of the results:
figure
plot(0:20,Eigenvalues,'LineStyle',':','Marker','.')
xlabel('k')
ylabel('E_k')
legend('c=1','c=2','c=3','c=3','c=5');
title('Eigenvalues')

figure
plot(x,Eigenfunction)
xlabel('x')
ylabel('y_0')
legend('c=1','c=2','c=3','c=3','c=5');
title('First eigenfunction')

end

function r=p(x)
r=1;
end

function r=q(c,x)
r=c*cos(x);
end

function r=w(x)
r=1;
end