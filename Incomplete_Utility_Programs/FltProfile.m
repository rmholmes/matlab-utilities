function [T,D,VarP] = FltProfile(flt_pos,hname,dname,VarOp,domain,h_filter,v_filter,tposvec)
%This function plots an along front-depth slice of a field
%following the flt position vector pos, with the floats on top.
%
%T = along track distance (m)
%
%D = depth
%
%VarP = Variable
%
%hname = filename of history file (variables without underscores)
%
%dname = filename of diagnostics file (variables with underscores)
%
%flt_pos = time series of float positions
%((flt)x(time)x(time,Xgrid,Ygrid,depth,lon,lat)) (positions must be at history field times)
%
%VarOp = pcolPlotNC field operation cell.
%
%domain = [xmin ymin zmin tmin; xmax ymax zmax tmax]; All total
%history file dimension min and max to extract. 
%
%h_filter = radius of circular horizontal box filter window. (2
%elements for two spatial filters)
%
%v_filter = # pts of vertical box filter window.
%
%tposvec = times (indices of flt_pos time dimension) to plot float positions.
%
%filenames:
ncid = netcdf.open(hname,'NC_NOWRITE');
ncidD = netcdf.open(dname,'NC_NOWRITE');
fID = netcdf.inqVarID(ncid,'f');
lonID = netcdf.inqVarID(ncid,'lon');
latID = netcdf.inqVarID(ncid,'lat');

svec = domain(2,:)-domain(1,:)+1;
xL = svec(1);yL=svec(2);zL=svec(3);tL=svec(4);
fL = length(flt_pos(:,1,1));

domain4D = ['[' num2str(domain(1,:)-1) '],[' num2str([xL yL zL tL]) ']'];
domain2D = ['[' num2str(domain(1,1:2)-1) '],[' num2str([xL yL]) ']'];
zvec = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'z'),[domain(1,3)-1],[zL],'double');

%%%%%Construct field string:
if (length(strfind(VarOp{end}, '(1)')) == 0)
    numvars = length(VarOp);
    VarOp{numvars+1} = '(1)';
    for i=2:numvars
        VarOp{numvars+1} = [VarOp{numvars+1} ' + (' num2str(i) ...
                            ')'];
    end
end
GetVarStr = VarOp{end};
for i=1:(length(VarOp)-1)
    %Get Variable i:
    
    %%%Un-built variables;%%%%
    if (length(strfind(VarOp{i}, '_')) == 0) %History variable
        GetVar = ['(netcdf.getVar(ncid,netcdf.inqVarID(ncid,''' VarOp{i} '''),' domain4D ...
                  ',''double''))'];
    else %Diagnostic variable
        GetVar = ['(netcdf.getVar(ncidD,netcdf.inqVarID(ncidD,''' ...
                  VarOp{i} '''),' domain4D ',''double''))'];
    end             
    %Replace String:
    GetVarStr = strrep(GetVarStr,['(' num2str(i) ')'],GetVar);
end
    
GetVarStr = strrep(GetVarStr,'x_rho',['(pi*6371000/180' ...
                    '*(netcdf.getVar(ncid,lonID,' domain2D ',''double'')+132).*cos(pi/180' ...
                    '*netcdf.getVar(ncid,latID,' domain2D ',''double'')))']);
GetVarStr = strrep(GetVarStr,'y_rho',['(pi*6371000/180*' ...
                    'netcdf.getVar(ncid,latID,' domain2D ',''double''))']);
    
%%If want f in calculations:
GetVarStr = strrep(GetVarStr,'f_cor',['repmat((netcdf.getVar(ncid,fID,' ...
                    domain2D ',''double'')),[1 1 ' num2str([zL tL]) ...
                    '])']);

%Get variable:
status = 'extracting variable...'
eval(['Var = squeeze(' GetVarStr ');']);
%Get rho:
eval(['rho = (netcdf.getVar(ncid,netcdf.inqVarID(ncid,''rho''),' domain4D ...
                  ',''double''));']);
%Filter in vertical and horizontal:
status = 'filtering variable...'
for t = 1:tL
        Var(:,:,:,t) = filter_field(filter_field(filter_field(Var(:,:,:,t),h_filter(1),['-' ...
                            's']),h_filter(2),'-s'),v_filter,'-v');
        rho(:,:,:,t) = filter_field(filter_field(rho(:,:,:,t),h_filter(1),['-' ...
                            's']),h_filter(2),'-s');
end
COM = [round(mean(flt_pos(:,:,2)-domain(1,1)+1,1))' ...
       round(mean(flt_pos(:,:,3)-domain(1,2)+1,1))'];
[D,T] = meshgrid(zvec,flt_pos(1,:,1)');
VarP = zeros(size(D));
rhoP = zeros(size(D));
for t = 1:tL
    VarP(t,:) = permute(Var(COM(t,1),COM(t,2),:,t),[1 3 2 4]);
    rhoP(t,:) = permute(rho(COM(t,1),COM(t,2),:,t),[1 3 2 4]);
end
pcolPlot(interp2(T,2),interp2(D,2),interp2(VarP,2));
hold on;
contour(interp2(T,2),interp2(D,2),interp2(rhoP,2),[20:0.2:30],['-' ...
                    'k'],'LineWidth',1);
hold on;
plot(flt_pos(1,tposvec,1),flt_pos(:,tposvec,4),'.m','MarkerSize', ...
     15);
hold on;
plot(flt_pos(1,tposvec,1),mean(flt_pos(:,tposvec,4),1),'.k','MarkerSize',25);
xlabel('Day');
ylabel('Depth');

end

