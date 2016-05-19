function A = FmovieNCFLT(fname_VARS,fname_FLTS,VarOp,slice,clim,Flts)
%This function makes and plays a movie of the Variable given by VAROp
%(see pcolPlotNC.m) taken from a netcdf file fname_VARS. The inputs
%are the same as pcolPlotNC.m (except the t slice entry is
%irrelevant/unnecessary). clim is a two element vector with the caxis
%limits. In addition, the function plots the positions of the
%floats FLTS contained in fname_FLTS.

%%%%%IN-FILE options%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IF WANT TO RUN BACKWARDS IN TIME
RBACKWARDS = 0;

%Initial time:
ti = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%Load files and Flts Vars%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncidFLT = netcdf.open(fname_FLTS,'NC_NOWRITE');
NFID = netcdf.inqDimID(ncidFLT,'drifter');
[temp, NF]=netcdf.inqDim(ncidFLT,NFID);
ftimeID = netcdf.inqVarID(ncidFLT,'ocean_time');
ftime = netcdf.getVar(ncidFLT,ftimeID,'double');
lonFID = netcdf.inqVarID(ncidFLT,'lon');
latFID = netcdf.inqVarID(ncidFLT,'lat');
depthFID = netcdf.inqVarID(ncidFLT,'depth');

%Make his_time_step:
ncid = netcdf.open(fname_VARS,'NC_NOWRITE');
timeID = netcdf.inqDimID(ncid,'timeD');
[temp, timeL] = netcdf.inqDim(ncid,timeID);
HtimeID = netcdf.inqVarID(ncid,'time');
Htime = netcdf.getVar(ncid,HtimeID,'double');
netcdf.close(ncid);

if (Flts ~= 0)
NF = length(Flts);
end

%Make lon/lat/dep arrays:
lonarr = zeros(timeL,NF);
latarr = lonarr;
deparr = lonarr;
for i = 1:timeL
    ['processing floats for time' num2str(i)]
        %Closest float time to his time:
    [temp tf] = min(abs(ftime - Htime(i)));
    if (Flts == 0)
        lonarr(i,:) = netcdf.getVar(ncidFLT,lonFID,[0 tf-1],[NF 1], ...
                                            'double');
        latarr(i,:) = netcdf.getVar(ncidFLT,latFID,[0 tf-1],[NF 1], ...
                                            'double');
        deparr(i,:) = netcdf.getVar(ncidFLT,depthFID,[0 tf-1],[NF ...
                            1], 'double');
    else
        for f = 1:NF
            fn = Flts(f);
            lonarr(i,f) = netcdf.getVar(ncidFLT,lonFID,[fn-1 tf-1],[1 ...
                                1], 'double');
            latarr(i,f) = netcdf.getVar(ncidFLT,latFID,[fn-1 tf-1],[1 ...
                                1], 'double');
            deparr(i,f) = netcdf.getVar(ncidFLT,depthFID,[fn-1 tf- ...
                                1],[1 1], 'double');
        end
    end
end
netcdf.close(ncidFLT);
%Do movie:
% $$$ fig1 = figure;%('visible','off');
% $$$ set(fig1, 'Position', get(0,'Screensize')); % Maximize figure. 
% $$$ winsize = get(fig1,'Position');

A = moviein(timeL);%,fig1,winsize);
for td=ti:timeL
    if (RBACKWARDS == 1)
        t = timeL-td+1;
    else
        t = td;
    end
    clf;
    pcolPlotNC(fname_VARS,VarOp,[slice t],1);
    caxis(clim);
    hold on;
    if (length(slice) == 1) %Isopycnal
        plot(lonarr(t,:),latarr(t,:),'*k','MarkerSize',20);
    elseif (slice(1) == 0 && slice(2) == 0)
        plot(lonarr(t,:),latarr(t,:),'*k','MarkerSize',20);
    elseif (slice(1) == 0 && slice(3) == 0)
        plot(lonarr(t,:),deparr(t,:),'*k','MarkerSize',20);
    elseif (slice(2) == 0 && slice(3) == 0)
        plot(latarr(t,:),deparr(t,:),'*k','MarkerSize',20);
    else
        status = 'BAD INPUT!'
    end
    %    axis([-220 -90 -20 20]);
    hold off;
    A(:,td-ti+1)=getframe(gcf);%(gcf);
end
end
