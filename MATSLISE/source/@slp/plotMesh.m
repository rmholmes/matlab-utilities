function plotMesh(slp,mesh)
%SLP: plotmesh(slpObject,meshData)
% generates a visualization of the mesh which was constructed by Matslise
% for the slp-object passed as first input argument. The input argument
% meshData is a structure returned by the method computeEigenvalues.
a=mesh.Infa;b=mesh.Infb;
if liouvilleNormalForm(slp)
    %Schrodinger problem
    step = min(abs(b-a)/5000,0.05);
    while (b-a)/step >10000
        step=step*2;
    end
    w = a:step:b;
    if mesh.halfRangeReduction
      x = [0 cumsum(mesh.h)];   
    else
      x = [a a + cumsum(mesh.h)];
    end
    w = sort([w x]);
    Y=arrayfun(slp.q,w);  
    legendStr{1} = 'meshpoints';
    figure;
    my=min(Y);
    phandles= plot(x,my*ones(1,length(x)),'bo','MarkerSize',5,'MarkerFaceColor','b');
    hold on;
    legendStr{2} = 'x_{match}';
    phandles = [phandles plot(x(mesh.imatch+1),my,'rs','MarkerSize',5,'MarkerFaceColor','r')];
    phandles = [phandles plot(w,Y,'k:')];
    legendStr{3} = 'Q(x)';
    legendStr{4} = 'constant approx.';
    cvbar=[];
    cx=[];
    for i=1:length(mesh.V0)
        cx=[cx x(i) x(i+1)];
        cvbar=[cvbar mesh.V0(i) mesh.V0(i)];
    end
    phandles= [phandles plot(cx,cvbar,'b')];
    hold off;
    legend(phandles(1:end),legendStr,'Location','Best');
    xlabel('x');
    title('The reference potential function + the meshpoints');
else
    disp('This functionality has not been implemented yet for a problem not in Liouville normal form')
end
axis 'tight'; 