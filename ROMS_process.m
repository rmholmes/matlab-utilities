function ROMS_process(hname)
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%
% This function processes a raw ROMS output _his file and creates an
% extras file that adds various useful/common variables, including
% z_rho, to add more... The filename and location is always hname
% _ext.nc and only contains rho variables and dimensions, although
% the dimensions etc. should only be read from hname.
%
% INPUTS:
%
% hname = _his.nc history file generated by ROMS.
%
%
%---------------------------------------------------------------------
%
% Dependencies; none. 
%
% Ryans ROMS Matlab and netcdf Utilities 24/7/13
%
%---------------------------------------------------------------------
%---------------------------------------------------------------------

status = ['This is ROMS_process running on the file ' hname]

%Processed file information -----------------------------------------
ncid = netcdf.open(hname,'NC_NOWRITE');
%Get dimensions:
xiID = netcdf.inqDimID(ncid, 'xi_rho');
etaID = netcdf.inqDimID(ncid, 'eta_rho');
sID = netcdf.inqDimID(ncid, 's_rho');
swID = netcdf.inqDimID(ncid,'s_w');
tID = netcdf.inqDimID(ncid, 'ocean_time');

%Lengths of file:
[tmp,xLhis] = netcdf.inqDim(ncid,xiID);
[tmp,yLhis] = netcdf.inqDim(ncid,etaID);
[tmp,zLhis] = netcdf.inqDim(ncid,sID);
[tmp,tLhis] = netcdf.inqDim(ncid,tID);
netcdf.close(ncid);

Outfilename = [hname(1:(end-3)) '_ext.nc'];
if (exist(Outfilename) == 0) 
%Create new extras file ----------------------------------------------
status = ['Creating new extras file...']
ncidOUT = netcdf.create(Outfilename,'NETCDF4');

xiID_OUT = netcdf.defDim(ncidOUT,'xi_rho',xLhis);
etaID_OUT = netcdf.defDim(ncidOUT,'eta_rho',yLhis);
sID_OUT = netcdf.defDim(ncidOUT,'s_rho',zLhis);
swID_OUT = netcdf.defDim(ncidOUT,'s_w',zLhis+1);
tID_OUT = netcdf.defDim(ncidOUT,'ocean_time',tLhis);
D4D_OUT = [xiID_OUT etaID_OUT sID_OUT tID_OUT]; %quick 4D dims.
D4Dw_OUT = [xiID_OUT etaID_OUT swID_OUT tID_OUT];

%define new variables -------------------------------------
%----------------------------------------------------------

%z_rho:
z_rhoID_OUT = netcdf.defVar(ncidOUT,'z_rho','double',D4D_OUT);

%z_w:
z_wID_OUT = netcdf.defVar(ncidOUT,'z_w','double',D4Dw_OUT);

% $$$ %Sf:
% $$$ SFID_OUT = netcdf.defVar(ncidOUT,'SF','double',[xiID_OUT etaID_OUT ...
% $$$                     tID_OUT]);
% $$$ 
% $$$ %TIVtop and TIVbot:
% $$$ TIVtopID_OUT = netcdf.defVar(ncidOUT,'TIVtop','double',[xiID_OUT etaID_OUT ...
% $$$                     tID_OUT]);
% $$$ TIVbotID_OUT = netcdf.defVar(ncidOUT,'TIVbot','double',[xiID_OUT etaID_OUT ...
% $$$                     tID_OUT]);

%----------------------------------------------------------
%----------------------------------------------------------

%Write some descriptions/global attributes:
netcdf.putAtt(ncidOUT,netcdf.getConstant('NC_GLOBAL'), ...
              'Description',['This is a netcdf extras file created ' ...
                    'by ROMS_process.m']);
netcdf.putAtt(ncidOUT,netcdf.getConstant('NC_GLOBAL'), ...
              'Base file',hname);
netcdf.putAtt(ncidOUT,netcdf.getConstant('NC_GLOBAL'), ...
              'File_creation_data_time',clock);

%End define mode adding a 20000 byte buffer on the end of the header:
netcdf.endDef(ncidOUT,20000,4,0,4);
status = ['Finished creating new file']
elseif (exist(Outfilename) == 2)
%Append to file new variables------------------------------
status = ['Defining new variables in extras file...']
ncidOUT = netcdf.open(Outfilename,'NC_WRITE');

xiID_OUT = netcdf.inqDimID(ncidOUT,'xi_rho');
etaID_OUT = netcdf.inqDimID(ncidOUT,'eta_rho');
sID_OUT = netcdf.inqDimID(ncidOUT,'s_rho');
swID_OUT = netcdf.inqDimID(ncidOUT,'s_w');
tID_OUT = netcdf.inqDimID(ncidOUT,'ocean_time');
D4D_OUT = [xiID_OUT etaID_OUT sID_OUT tID_OUT]; %quick 4D dims.
D4Dw_OUT = [xiID_OUT etaID_OUT swID_OUT tID_OUT];

%define new variables -------------------------------------
%----------------------------------------------------------
netcdf.reDef(ncidOUT);
z_wID_OUT = netcdf.defVar(ncidOUT,'z_w','double',D4Dw_OUT);
netcdf.endDef(ncidOUT,20000,4,0,4);
% $$$ SFID_OUT = netcdf.inqVarID(ncidOUT,'SF');
% $$$ TIVtopID_OUT = netcdf.inqVarID(ncidOUT,'TIVtop');
% $$$ TIVbotID_OUT = netcdf.inqVarID(ncidOUT,'TIVbot');

%----------------------------------------------------------
%----------------------------------------------------------
status = ['Finished appending to file']
end

%Put new variables-----------------------------------------
%----------------------------------------------------------

%%%%%Z_RHO:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
status = ['Getting and putting z_rho...']
%Get z_rho:
for t = 1:tLhis
    status = ['time indicy = ' num2str(t) ' of ' num2str(tLhis)]
    z_rho = get_depths(hname,1,0,t,0);
    netcdf.putVar(ncid,z_rhoID_OUT,[0 0 0 t-1],[xLhis yLhis zLhis 1], ...
                  z_rho);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Z_W:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
status = ['Getting and putting z_w...']
%Get z_w:
for t = 1:tLhis
    status = ['time indicy = ' num2str(t) ' of ' num2str(tLhis)]
    z_w = get_depths(hname,5,0,t,0);
    netcdf.putVar(ncid,z_wID_OUT,[0 0 0 t-1],[xLhis yLhis zLhis+1 1], ...
                  z_w);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%SF:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ status = ['Getting and putting streamfunction...']
% $$$ lon = ncread(hname,'lon_rho');
% $$$ lat = ncread(hname,'lat_rho');
% $$$ lonrange = 141:601;%-205 to -90
% $$$ latrange = 66:188;%-15 to 15
% $$$ lon = lon(lonrange,latrange);
% $$$ lat = lat(lonrange,latrange);
% $$$ SF = NaN*zeros(xLhis,yLhis);
% $$$ 
% $$$ for t = 1:tLhis
% $$$     status = ['time indicy = ' num2str(t) ' of ' num2str(tLhis)]
% $$$     u = GetVar(hname,0,{'u'},{[lonrange(1) lonrange(end)],[latrange(1) ...
% $$$                         latrange(end)],[20 20],[t t]});
% $$$     v = GetVar(hname,0,{'v'},{[lonrange(1) lonrange(end)],[latrange(1) ...
% $$$                         latrange(end)],[20 20],[t t]});
% $$$     [SF(lonrange,latrange),tmp] = HelmD(lon,lat,u,v);
% $$$     netcdf.putVar(ncid,SFID_OUT,[0 0 t-1],[xLhis yLhis 1],SF);
% $$$ end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%TIVtop and TIVbot:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ status = ['Getting and putting TIVtop/TIVbot...']
% $$$ rhotop = 23;
% $$$ rhobot = 24.6;
% $$$ for t = 1:tLhis
% $$$     status = ['time indicy = ' num2str(t) ' of ' num2str(tLhis)]
% $$$     TIVtop = GetVar(hname,0,{'z_rho'},{[1 xLhis],[1 yLhis],rhotop,[t t]});
% $$$     TIVbot = GetVar(hname,0,{'z_rho'},{[1 xLhis],[1 yLhis],rhobot,[t t]});
% $$$     netcdf.putVar(ncid,TIVtopID_OUT,[0 0 t-1],[xLhis yLhis 1],TIVtop);
% $$$     netcdf.putVar(ncid,TIVbotID_OUT,[0 0 t-1],[xLhis yLhis 1],TIVbot);
% $$$ end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------------------------------------------
%----------------------------------------------------------
netcdf.close(ncidOUT);
status = 'Done!'

end
