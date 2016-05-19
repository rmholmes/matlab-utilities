function [Var,C,Z,Hsbl,Corners] = PlotAFront(fname_VARS,fname_DIAS,VarOp,Front,W,Res,flip,plotting)
%This function plots an along-front average in C,Z (across
%front-depth) space of the variable defined by VarOp (see
%pcolPlotNC.nc) for the front defined by Front = [lon_w lat_w lon_e
%lat_e depth time] (western and eastern ends of front + time of front)
%(the depth is not used). W is the half-width of the along-front
%direction and Res = [N_g N_c N_z] is the resolution in the three
%directions.
%
%flip = flip direction of front before calculating.
%plotting = 1 -> plot, otherwise don't. 
    
%VarOp form:
%VarOp = {'Var1' 'Var2' '{(1)::C::(2)}'} - plot cross-front
%component of rotated vectors with u-comp = (1) and v-comp =
%(2). '{(1)::A::(2)}' = along front component.
%VarOp = {'Var1'} - no rotation.
%VarOp = {'Var1' 'Var2'} - defaults to along front
% VarOp = {'u' 'v' 'w'} - plot along front vel. in color and other's
% as vectors.
%
%Can use 'x_rho','y_rho','f_cor'
    
%%%Fliping of x-axis:
    if (flip == 1)
        tmp1 = Front(1);
        tmp2 = Front(2);
        Front(1) = Front(3);
        Front(2) = Front(4);
        Front(3) = tmp1;
        Front(4) = tmp2;
    end
    
%%%%IN-FILE options%%%%%
    zst = 10; %Starting z-level.
    dmn_bdy = 0.5; %buffer length (degrees) on sides of rotated
                   %domain. 
%Plotting across-front vectors;
    scale = 10000;
    ypts = 25;
    zpts = 35;
    
%%%%%%%%Load FILE and Variable IDs%%%%%%%%
    ncid = netcdf.open(fname_VARS,'NC_NOWRITE');
    if (fname_DIAS ~= 0)
    ncidD = netcdf.open(fname_DIAS,'NC_NOWRITE');
    end
    if ((length(fname_VARS) >= 12) && strcmp(fname_VARS((end-11):end), ...
        'ocean_his.nc')) 
                                                    %IF USING
                                                    %ORIGINAL ROMS
                                                    %OUTPUT (ONLY
                                                    %TRY WITH rho
                                                    %VARS). 
    %Dimensions:
    xiID = netcdf.inqDimID(ncid,'xi_rho');
    etaID = netcdf.inqDimID(ncid,'eta_rho');
    sID = netcdf.inqDimID(ncid,'s_rho');
    timeID = netcdf.inqDimID(ncid,'ocean_time');
    %Lon,Lat,z:
    lonID = netcdf.inqVarID(ncid,'lon_rho');
    latID = netcdf.inqVarID(ncid,'lat_rho');
    zID = netcdf.inqVarID(ncid,'s_rho');
    rhoID = netcdf.inqVarID(ncid,'rho');
    fID = netcdf.inqVarID(ncid,'f');
    else                                           %IF USING
                                                   %COMPILED DATA
                                                   %FILES
                                                   %(make_vars.m) 
    %Dimensions:
    xiID = netcdf.inqDimID(ncid,'xiD');
    etaID = netcdf.inqDimID(ncid,'etaD');
    sID = netcdf.inqDimID(ncid,'sD');
    timeID = netcdf.inqDimID(ncid,'timeD');
    %Lon,Lat,z:
    lonID = netcdf.inqVarID(ncid,'lon');
    latID = netcdf.inqVarID(ncid,'lat');
    zID = netcdf.inqVarID(ncid,'z');
    rhoID = netcdf.inqVarID(ncid,'rho');
    %HsblID = netcdf.inqVarID(ncid,'Hsbl');
    fID = netcdf.inqVarID(ncid,'f');
    end
    
% $$$     %Get zeta for PG calculation:
% $$$     ncidzeta = netcdf.open(['/mnt/Data2/ryan/eqROMS/EQROMSout/eq20/' ...
% $$$                         'run2/ocean_his.nc'],'NC_NOWRITE');
% $$$     zetaID = netcdf.inqVarID(ncidzeta,'zeta');

    %Clean:
    [temp, xiL]=netcdf.inqDim(ncid,xiID);
    [temp, etaL] = netcdf.inqDim(ncid,etaID);
    [temp, zL] = netcdf.inqDim(ncid,sID);
    [temp, timeL] = netcdf.inqDim(ncid,timeID);
    clear temp xiID etaID sID timeID;


%%%%%%%%%%%Get operation string and variable strings%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
comp = '2D';
if (length(strfind(VarOp{end}, '(1)')) == 0)
        numvars = length(VarOp);
        if (strcmp(VarOp{1},'u') & strcmp(VarOp{2},'v') & ...
            strcmp(VarOp{3},'w') )
            comp = '3D';
            VarOp{numvars+1} = '(1):::(2):;:(3)';
        elseif (numvars == 2)
            VarOp{3} = '{(1)::A::(2)}';
        else
        VarOp{numvars+1} = '(1)';
        for i=2:numvars

            VarOp{numvars+1} = [VarOp{numvars+1} ' + (' num2str(i) ...
                                ')'];
        end
        end
    end

%%%%Get VarNAME for field variables:
    IDs = zeros(length(VarOp)-1,1);
    VarNAME = VarOp{end};
    for i=1:(length(VarOp)-1)
        VarNAME = strrep(VarNAME,['(' num2str(i) ')'],VarOp{i});
    end

%%%%%%%Make Rotated Grid%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CF = pi*6371000/180; %conversion factor from lat/lon to distance.
lon_w = Front(1);
lat_w = Front(2);
lon_e = Front(3);
lat_e = Front(4);%Frontal ends
x_w = CF*(lon_w+132)*cos(pi/180*lat_w);
y_w = CF*lat_w;
x_e = CF*(lon_e+132)*cos(pi/180*lat_e);
y_e = CF*lat_e;%Frontal ends cartesian
N_c = Res(2);
N_g = Res(1);%Resolutions
if (length(Res) == 3)
    N_z = Res(3);
else
    N_z = 0;
end
lc = 2*W/N_c;
lg = sqrt((x_e-x_w)^2+(y_e-y_w)^2)/N_g; %length of space-steps
theta = atan2(y_e-y_w,x_e-x_w); %angle of front to eastward
lon_rot = zeros(N_g+1,N_c+1); 
lat_rot = lon_rot; %initialize lon and lat of rotated grid.
for i=1:(N_g+1)
    for j=1:(N_c+1)
phi = atan2((j-(N_c/2+1))*lc,(i-1)*lg); %angle of point from
                                        %frontal line. 
leng = sqrt(((j-(N_c/2+1))*lc)^2+((i-1)*lg)^2); %distance from
                                                %point to east
                                                %end. 
lat_rot(i,j) = (y_w+sin(phi+theta)*leng)/CF; %lat of point
lon_rot(i,j) = (x_w+cos(phi+theta)*leng)/CF/ ...
    cos(pi/180*lat_rot(i,j))-132; %lon of point
    end
end
Corners = [lon_rot(1,1) lat_rot(1,1); lon_rot(1,end) lat_rot(1,end);...
           lon_rot(end,end) lat_rot(end,end); lon_rot(end,1) lat_rot(end,1)];


%%%%%%%%%%%Get Field to plot%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Get smaller domain:
lon = netcdf.getVar(ncid,lonID,'double');
lat = netcdf.getVar(ncid,latID,'double');
t = Front(end);
[tmp mnln] = min(abs(lon(:,1)-min(min(lon_rot))+dmn_bdy));
[tmp mxln] = min(abs(lon(:,1)-max(max(lon_rot))-dmn_bdy));%find minimum
                                                      %and maximum
                                                      %lon/lat
                                                      %indicies
                                                      %needed.
[tmp mnlt] = min(abs(lat(1,:)-min(min(lat_rot))+dmn_bdy));
[tmp mxlt] = min(abs(lat(1,:)-max(max(lat_rot))-dmn_bdy));
xiL = mxln-mnln+1;
etaL = mxlt-mnlt+1; %small domain lengths
domain = '[mnln-1 mnlt-1 zst-1 t-1],[xiL etaL zL-zst+1 1]';
domain2D = '[mnln-1 mnlt-1],[xiL etaL]';%domain vectors
domain3D = '[mnln-1 mnlt-1 t-1],[xiL etaL 1]';

lon = lon(mnln:mxln,mnlt:mxlt);
lat = lat(mnln:mxln,mnlt:mxlt);%shrink lon/lat.

%Get Variable:
    GetVarStr = VarOp{end};
    for i=1:(length(VarOp)-1)
        %Get Variable i:
        
        %%%Un-built variables;%%%%
        if (strcmp(VarOp{i},'Z'))
            GetVar = '(Z)';
        elseif (strcmp(VarOp{i},'u_g')) %Geostrophic velocity
            GetVar = ['(bsxfun(@times,netcdf.getVar(ncidD,' ...
            'netcdf.inqVarID(ncidD,''v_prsgrd''),' domain ...
            ',''double''),1./(netcdf.getVar(ncidD,fID,' domain2D ...
             ',''double''))))'];
        elseif (strcmp(VarOp{i},'v_g')) %Geostrophic velocity
            GetVar = ['(bsxfun(@times,netcdf.getVar(ncidD,' ...
            'netcdf.inqVarID(ncidD,''u_prsgrd''),' domain ...
            ',''double''),1./(-netcdf.getVar(ncidD,fID,' domain2D ...
             ',''double''))))'];
        elseif (strcmp(VarOp{i},'u_ag')) %Ageostrophic velocity
            GetVar = ['(bsxfun(@times,(netcdf.getVar(ncidD,' ...
            'netcdf.inqVarID(ncidD,''v_prsgrd''),' domain ...
            ',''double'')+netcdf.getVar(ncidD,' ...
            'netcdf.inqVarID(ncidD,''v_cor''),' domain ...
            ',''double'')),1./(-netcdf.getVar(ncidD,fID,' domain2D ...
             ',''double''))))'];
        elseif (strcmp(VarOp{i},'v_ag')) %Ageostrophic velocity
            GetVar = ['(bsxfun(@times,(netcdf.getVar(ncidD,' ...
            'netcdf.inqVarID(ncidD,''u_prsgrd''),' domain ...
            ',''double'')+netcdf.getVar(ncidD,' ...
            'netcdf.inqVarID(ncidD,''u_cor''),' domain ...
            ',''double'')),1./(netcdf.getVar(ncid,fID,' domain2D ...
             ',''double''))))'];
% $$$         elseif (strcmp(VarOp{i},'zeta'))
% $$$             GetVar = ['repmat(netcdf.getVar(ncidzeta,zetaID,' domain2D ...
% $$$                       ',''double''),[1 1 zL-zst+1])'];
            %%%Normal variables;%%%%
        else
            if (length(strfind(VarOp{i}, '_')) == 0) %History variable
        GetVar = ['(netcdf.getVar(ncid,netcdf.inqVarID(ncid,''' VarOp{i} '''),' domain ...
                  ',''double''))'];
            else %Diagnostic variable
                GetVar = ['(netcdf.getVar(ncidD,netcdf.inqVarID(ncidD,''' ...
                          VarOp{i} '''),' domain ',''double''))'];
            end             
        end
        %Replace String:
        GetVarStr = strrep(GetVarStr,['(' num2str(i) ')'],GetVar);
    end
    %%If want x and y and z in calculations:
    zvec = netcdf.getVar(ncid,zID,[zst-1],[zL-zst+1],'double');
    GetVarStr = strrep(GetVarStr,'x_rho',['repmat((pi*6371315/180' ...
             '*netcdf.getVar(ncid,lonID,' domain2D ',''double'').*cos(pi/180' ...
             '*netcdf.getVar(ncid,latID,' domain2D ',''double''))),[1 1 zL-zst+1])']);
    GetVarStr = strrep(GetVarStr,'y_rho',['repmat((pi*6371315/180*' ...
                        'netcdf.getVar(ncid,latID,' domain2D ',''double'')),[1 1 ' ...
                        'zL-zst+1])']);
    
    %%If want f in calculations:
    GetVarStr = strrep(GetVarStr,'f_cor',['repmat((netcdf.getVar(ncid,fID,' ...
                        domain2D ',''double'')),[1 1 zL-zst+1])']);

    
    %%%%%%Do 2D rotation replacements:
    Asplit = strfind(GetVarStr,'::A::'); %Find along front
                                         %projections.
    while (length(Asplit) ~= 0)
    
        Asplitst = strfind(GetVarStr(1:(Asplit(1))),'{');
        Asplitst = Asplitst(end);
        Asplitend = strfind(GetVarStr((Asplit(1)):end),'}');
        Asplitend = Asplitend(1) + Asplit(1)-1;
        uvar = GetVarStr((Asplitst+1):(Asplit(1)-1));
        vvar = GetVarStr((Asplit(1)+5):(Asplitend-1));
        GetVarStr = strrep(GetVarStr,['{' uvar '::A::' vvar '}'], ...
                           ['-(cos(atan2(' vvar ',' uvar ')-theta).*sqrt((' ...
                            uvar ').^2+(' vvar ').^2))']);
        Asplit = strfind(GetVarStr,'::A::');
    end    
    
    Csplit = strfind(GetVarStr,'::C::'); %Find across front
                                         %projections.
    while (length(Csplit ~= 0))

        Csplitst = strfind(GetVarStr(1:(Csplit(1))),'{');
        Csplitst = Csplitst(end);
        Csplitend = strfind(GetVarStr((Csplit(1)):end),'}');
        Csplitend = Csplitend(1) + Csplit(1)-1;
        uvar = GetVarStr((Csplitst+1):(Csplit(1)-1));
        vvar = GetVarStr((Csplit(1)+5):(Csplitend-1));
        GetVarStr = strrep(GetVarStr,['{' uvar '::C::' vvar '}'], ...
                        ['(sin(atan2(' vvar ',' uvar ')-theta).*sqrt((' ...
                         uvar ').^2+(' vvar ').^2))']);
        Csplit = strfind(GetVarStr,'::C::');
    
    end

%%%%%%Get Variable and axis on Rotated Grid%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~strcmp(comp,'3D'))
eval(['VarIN = squeeze(' GetVarStr ');']);
Var = zeros(N_c+1,zL-zst+1);
for z=1:(zL-zst+1)
    Var(:,z) = mean(interp2(lon',lat',VarIN(:,:,z)',lon_rot',lat_rot'),2);
end

%%%vectors for 2nd and third component:
else
    splitter1 = strfind(GetVarStr,':::');
    splitter2 = strfind(GetVarStr,':;:');
eval(['VarU = squeeze(' GetVarStr(1:(splitter1-1)) ');']);
eval(['VarV = squeeze(' GetVarStr((splitter1+3):(splitter2-1)) ');']);
eval(['VarWIN = squeeze(' GetVarStr((splitter2+3):end) ');']);
VarGIN = cos(atan2(VarV,VarU)-theta).*sqrt(VarU.^2+VarV.^2);
VarCIN = sin(atan2(VarV,VarU)-theta).*sqrt(VarU.^2+VarV.^2);
Var = zeros(N_c+1,zL-zst+1);
VarC = Var;
VarW = Var;
for z=1:(zL-zst+1)
    Var(:,z) = mean(interp2(lon',lat',VarGIN(:,:,z)',lon_rot', lat_rot'),2);
    VarC(:,z) = mean(interp2(lon',lat',VarCIN(:,:,z)',lon_rot', ...
                             lat_rot'),2);    
    VarW(:,z) = mean(interp2(lon',lat',VarWIN(:,:,z)',lon_rot', ...
                             lat_rot'),2);    
end
end

%if (plotting == 1)
%%%Get density for isopycnals and mixed-layer depth%%%%%%
eval(['rho = netcdf.getVar(ncid,rhoID,' domain ',''double'');']);
%eval(['Hsbl = netcdf.getVar(ncid,HsblID,' domain3D ',''double'');']);
%Hsbl = mean(interp2(lon',lat',Hsbl',lon_rot',lat_rot'),2);
rhoR = zeros(N_c+1,zL-zst+1);
for z=1:(zL-zst+1)
    rhoR(:,z) = mean(interp2(lon',lat',rho(:,:,z)',lon_rot',lat_rot'),2);
end
%end

%Condition on rho for FrontForces 9-4-13:
% $$$ Var = Var.*(rhoR<24.6);

%%%Create C,Z vectors and do z-interp if wanted%%%%%
[Z,C]=meshgrid(netcdf.getVar(ncid,zID,[zst-1],[zL-zst+1],'double'),(-W:lc:W));
if (N_z > 0) %doing Z-interp
    zstval = netcdf.getVar(ncid,zID,[zst-1],[1],'double');
    lz = -zstval/N_z;
    [Zn,Cn]=meshgrid((zstval+lz/2):lz:(-lz/2),-W:lc:W);
    rhoR = interp2(C',Z',rhoR',Cn',Zn')';
    Var = interp2(C',Z',Var',Cn',Zn')';
% $$$     if (~strcmp(comp,'N'))
% $$$         VarC = interp2(C',Z',VarC',Cn',Zn')';
% $$$     end
    if (strcmp(comp,'3D'))
        VarW = interp2(C',Z',VarW',Cn',Zn')';
        VarC = interp2(C',Z',VarC',Cn',Zn')';
    end
    Z = Zn;
    C = Cn;
end

if (plotting == 1)
%%%DO PLOTTING:%%%%%%
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
pcolor(C,Z,Var);
xlab = xlabel('Across Front distance (m)');
ylab = ylabel('Depth (m)');
tit = title(['Front Slice of ' VarNAME ' ' comp ' Front; Lon = ' ...
             num2str((lon_w+lon_e)/2) ', Lat = ' num2str((lat_w+lat_e)/2) ...
             ', t = ' num2str(t) ' )']);
colormap(jet);
colorbar;
shading flat;
coloraxis=caxis;
hold on;
contour(C,Z,rhoR,20:0.2:27.4,'-k','LineWidth',1);
contour(C,Z,rhoR,[23 24.6],'-k','LineWidth',2);
caxis(coloraxis);
hold on;
%plot(C(:,1),Hsbl,'Color',[1 1 1],'LineWidth',2);
if (strcmp(comp,'3D')) %If rotated; plot cross-front vector field:
    for yp=1:ypts
        for zp=1:zpts
            zpos = floor(zp*length(Z(1,:))/(zpts+1));
            ypos = floor(yp*length(C(:,1))/(ypts+1));
            vel = VarC(ypos,zpos);
            velW = VarW(ypos,zpos);
            plot([C(ypos,zpos); C(ypos,zpos)+vel*scale],[Z(ypos,zpos); ...
                                Z(ypos,zpos)+velW*scale],'-','Color',[0 ...
                                0 0],'LineWidth',1);
            plot(C(ypos,zpos),Z(ypos,zpos),'*k','Color',[0 0 0],'MarkerSize',4);
        end
    end
end
% $$$ set(gca,'FontSize',20);
% $$$ set(xlab,'FontSize',20);
% $$$ set(ylab,'FontSize',20);
% $$$ set(tit,'FontSize',20);
%    axis([-210 -100 -15 15]);
%daspect([15 15 1]);
end
netcdf.close(ncid);
if (fname_DIAS ~= 0)
netcdf.close(ncidD);
end
Hsbl = 0;

end
