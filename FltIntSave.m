
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%
% This script interpolates a variable onto all float tracks in a
% float file and saves back into that float file.
%
%---------------------------------------------------------------------
%---------------------------------------------------------------------

fname = ['/mnt/Data4/ryan/PEQ_MRZV/MRZV_Ri_ROFF/' ...
         'MRZV_Ri_roff_floats.nc'];
hname = MRZV_his;
dname = 0;
h_filter = 0;
v_filter = 0;
t_filter = 0;
VarOp = {'v' 'Dz(1)'};
name = 'vSh';
neworold = 1; %1 = new variable, 0 = overwrite old variable.

%Get number of floats and setup out file.
ncid = netcdf.open(fname,'NC_WRITE');
drID = netcdf.inqDimID(ncid,'drifter');
%tID = netcdf.inqDimID(ncid,'ocean_time');
tID = netcdf.inqDimID(ncid,'ftime');
[tmp nflt] = netcdf.inqDim(ncid,drID);
[tmp tL] = netcdf.inqDim(ncid,tID);


'Getting flt_pos...'
flt_pos = FltPos(fname,1:nflt,0,0);

'Interpolating variable...'
VarOUT = FltInterp(flt_pos,hname,dname,VarOp,h_filter,v_filter, t_filter);

%Re-reverse time order:
VarOUT = fliplr(VarOUT);

%Make NaNs large:
VarOUT(isnan(VarOUT)) = 1e35;

'Putting variable...'
%Put variable:
if (neworold == 0)
    %Old variable:
    VarID = netcdf.inqVarID(ncid,name);
elseif (neworold == 1)
    %New variable:
    netcdf.reDef(ncid);
    VarID = netcdf.defVar(ncid,name,'double',[drID tID]);
    netcdf.putAtt(ncid,VarID,'Filtering',['h_filter = ' num2str(h_filter) ...
                        ' v_filter = ' num2str(v_filter)]);
    netcdf.endDef(ncid);
end

netcdf.putVar(ncid,VarID,VarOUT);
netcdf.close(ncid);
'Done!'