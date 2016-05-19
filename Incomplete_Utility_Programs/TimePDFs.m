
%This script compiles floats into sets corresponding to floats that
%enter a particular vortex at a particular time. It then can
%calculate propertie PDFs for these floats at different times
%relative to this time, and can calculate the track of the COM of
%the floats.
%
%Intended for use on basin and eq20 backwards offline floats
%calculations. 
%
%Ryan Holmes 17-11-12.

%%%%%%%%%%%%%%%%%%%%%%%%%%OPTIONS%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%BASIN RUN%%%%%%%%%%%
%Filenames:
fname = basinr1F;
fnameH = basinr1V;
%In-Vortex conditions:
NoV = 4;
thSF = [0.2 0.22 0.22 0.22]*1e6;
mxradius = 5.8e5;
lonHigh = 0;
%cPos:
load('basincPos');
cPos(cPos == 0) = NaN;
%Plot Settings:
PDFtimeSHF = [-100:1:0];
caxs = [-100 0];
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%EQ20 RUN%%%%%%%%%%%
%Filenames:
% $$$ fname = eq20r2F;
% $$$ fnameH = eq20r2V;
% $$$ %In-Vortex conditions:
% $$$ NoV = 3;
% $$$ thSF = [0.24 0.24 0.24]*1e6;
% $$$ mxradius = 5.8e5;
% $$$ lonHigh = -121.8;
% $$$ %cPos:
% $$$ load('eq20cPos');
% $$$ cPos(cPos == 0) = NaN;
% $$$ cPos(abs(cPos)>1e15) = NaN;
% $$$ %Plot Settings:
% $$$ PDFtimeSHF = [-40:1:0];
% $$$ caxs = [-40 0];
% $$$ Ri = fliplr(ncread(fname,'Ri')); %For Richardson number Evolution.
% $$$ FGF = fliplr(ncread(fname,'FGF'));
% $$$ FGF(FGF<0) = 0;
% $$$ FGF = log10(FGF);
%%%%%%%%%%%%%%%%%%%%%%%%

%Threshold number of floats to count as a PDF:
thFLTS = 150;
thPDFs = 10;

%PDF PROPERTIES:
Bins = 300;
Filter = 9;

%Colormap:
% $$$ cc = fliplr(jet(length(PDFtimeSHF)-3));
% $$$ cc(1,:) = 0.0;
% $$$ ccext = zeros(length(PDFtimeSHF),3);
% $$$ ccext(4:end,:) = cc;
% $$$ ccext(3,:) = 0.4;
% $$$ ccext(2,:) = 0.6;
% $$$ ccext(1,:) = 0.8;
% $$$ cc = ccext;
% $$$ cc = flipud(cc);

cc = fliplr(jet(length(PDFtimeSHF)));
cc(1,:) = 0.0;
% $$$ ccext = zeros(length(PDFtimeSHF),3);
% $$$ ccext(4:end,:) = cc;
% $$$ ccext(3,:) = 0.4;
% $$$ ccext(2,:) = 0.6;
% $$$ ccext(1,:) = 0.8;
% $$$ cc = ccext;
cc = flipud(cc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%WORKING CODE%%%%%%%%%%%%
%Load variables: (assumed all backwards)
'loading variables'

%Times:
timeF = flipud(ncread(fname,'ocean_time'));
tLF = length(timeF); %Total float times
timeH = ncread(fnameH,'time');
tLH = length(timeH); %Total History times

%Variables:
temp = -fliplr(ncread(fname,'temp'));
salt = -fliplr(ncread(fname,'salt'));
depth = fliplr(ncread(fname,'depth'));
lon = fliplr(ncread(fname,'lon'));
lat = fliplr(ncread(fname,'lat'));
PVv = fliplr(ncread(fname,'PVv'));
PVh = fliplr(ncread(fname,'PVh'));
rho = fliplr(ncread(fname,'rho'));
DBDX = fliplr(ncread(fname,'DBDX'));
DBDY = fliplr(ncread(fname,'DBDY'));
SF = fliplr(ncread(fname,'SF'));
Fnum = length(SF(:,1)); %Number of floats

% $$$ [lon(4,:)' lat(4,:)' depth(4,:)' temp(4,:)' salt(4,:)' PVv(4,:)' ...
% $$$  rho(4,:)' SF(4,:)' ocean_time/86400]

%Initialization positions:
lonINI = lon(:,end);
latINI = lat(:,end);
depthINI = depth(:,end);
%Invalidate bad values:
INVALID = abs(temp)>100;
temp(INVALID) = NaN;
salt(INVALID) = NaN;
depth(INVALID) = NaN;
lon(INVALID) = NaN;
lat(INVALID) = NaN;
PVv(INVALID) = NaN;
PVh(INVALID) = NaN;
rho(INVALID) = NaN;
DBDX(INVALID) = NaN;
DBDY(INVALID) = NaN;
SF(INVALID) = NaN;
PV = PVv+PVh;

%X and Y positions:
CF = pi*6371000/180;
x = CF*(lon+132).*cos(pi/180*lat);
y = CF*lat;
xINI = CF*(lonINI+132).*cos(pi/180*latINI);
yINI = CF*latINI;

%Get cPos and cPos positions: (vortex centers). 
xCpos = squeeze(CF*(cPos(1,:,:)+132).*cos(pi/180*cPos(2,:,:)));
yCpos = squeeze(CF*cPos(2,:,:));

%Identify vortex being entered:
'Calculating intialization time...'
INVALID = isnan(temp);
REDUCED = zeros(size(INVALID));
REDUCED(:,2:end) = INVALID(:,2:end) & (~INVALID(:,1:(end-1)));
TIME = repmat((1:tLF)',[1 Fnum]);
timeINIf = TIME(logical(REDUCED')); %Initialization time indicy
                                    %floats time.
if (timeF(end) == timeH(end) & timeF(1) == timeH(1))
timeINIh = timeINIf; %Initialization time indicy history time. 
same = 0;
elseif (timeF(end) == timeH(end) & length(timeF) == length(timeH)-1)
    timeINIh = timeINIf+1;
    same = 1;
else
    'DODDDGGGYYYY TIMEEEINNNGGG!'
end

radINI = zeros(Fnum,NoV); %Intialization radii from each vortex.
for V=1:NoV
radINI(:,V) = sqrt((xINI-xCpos(V,timeINIh)').^2+(yINI-yCpos(V,timeINIh)').^2);
end
[temporary Vortex] = min(radINI,[],2); %Vortex number of each
                                        %float.

%Make time series of radius, angle and SF diff from ultimately
%entered vortex for each float:
'Calculating radii, angle and SF diff...'
SFr = zeros(Fnum,tLF); 
SFa = zeros(Fnum,tLF);
SFd = zeros(Fnum,tLF);
for tf = 1:tLF
    th = tf+same;
    SFr(:,tf) = sqrt((x(:,tf)-xCpos(Vortex,th)).^2+(y(:,tf)- ...
                                                    yCpos(Vortex,th)).^2);
    SFd(:,tf) = SF(:,tf)-permute(cPos(3,Vortex,th),[2 1 3]);
    SFa(:,tf) = atan2(y(:,tf)-yCpos(Vortex,th),x(:,tf)-xCpos(Vortex,th));
end
%Radial buoyancy gradient:
DBDR = DBDX.*cos(SFa)+DBDY.*sin(SFa);
Gradhb = sqrt(DBDX.^2+DBDY.^2);
Vortex = repmat(Vortex,[1 tLF]);

'Calculating conditions...'
%%%%%%%%%%%%%%In or out of TIV:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cond = rho >= 23 & rho <= 24.6 & SFr < mxradius;
if (NoV == 3)
    cond = cond & ((Vortex == 1 & SFd < thSF(1)) | ... 
                   (Vortex == 2 & SFd < thSF(2)) | ...
                   (Vortex == 3 & SFd < thSF(3)));
elseif (NoV == 4)
    cond = cond & ((Vortex == 1 & SFd < thSF(1)) | ... 
                   (Vortex == 2 & SFd < thSF(2)) | ...
                   (Vortex == 3 & SFd < thSF(3)) | ...
                   (Vortex == 4 & SFd < thSF(4)));
end
%1 if in TIV

%%%%ENTRANCE TIMES:%%%%%%
condENT = zeros(size(cond));
condENT(:,2:end) = cond(:,2:end) & ~cond(:,1:(end-1)) & ...
    ~isnan(rho(:,1:(end-1))) & ~isnan(SFr(:,1:(end-1)));
%1 if just entered TIV correctly (TIV 'appearing' doesn't count).

condENTFIRST = condENT;
got1 = find(condENTFIRST(:,1));
gotENT = find(0);
for t = 1:tLF
    foundENT = ~(condENTFIRST(got1,t));
    ls = length(gotENT);
    lsf = length(foundENT);
        if (lsf ~= 0)
    gotENT((ls+1):(ls+lsf)) = got1(foundENT);
        end
    condENTFIRST(gotENT,t) = 0;
    got1 = find(condENTFIRST(:,t));
end
%1 if just entered TIV correctly for the FIRSTTTTT time.

[timeENTf,FloatENTf] = find(condENTFIRST');
VorENTf = Vortex(FloatENTf,1); %Vortices of these floats.

%timeENTf = time indicy just after a float entered a TIV correctly.
%valENTf = Corresponding float number.

%%%%COM calculation:

%Full COM:
COM = zeros(2*tLF,5);
COM(:,1) = (-tLF+1):1:tLF;
for t = 1:2*tLF
    time = COM(t,1);
    valtime = timeENTf+time;
    valfloat = FloatENTf(valtime>0 & valtime<=tLF);
    valtime = valtime(valtime>0 & valtime<=tLF);
    COMRad = SFr(sub2ind(size(SFr),valfloat,valtime));
    COMAng = SFa(sub2ind(size(SFa),valfloat,valtime));
    COM(t,4) = length(find(~isnan(COMRad)));
    COMAng = COMAng(~isnan(COMAng));
    COMRad = COMRad(~isnan(COMRad));
    COMx = COMRad.*cos(COMAng);
    COMy = COMRad.*sin(COMAng);
    COMx = mean(COMx);
    COMy = mean(COMy);
    COM(t,3) = atan2(COMy,COMx);
    COM(t,2) = sqrt(COMx^2+COMy^2);
    COMDep = depth(sub2ind(size(depth),valfloat,valtime));
    COM(t,5) = mean(COMDep(~isnan(COMDep)));
end

%COM EUC and NECC:
NECCFlts = find(0);
cnt = length(NECCFlts);
for t = -60:1:-15
    valtime = timeENTf+t;
    valfloat = FloatENTf(valtime>0 & valtime<=tLF);
    valtime = valtime(valtime>0 & valtime<=tLF);
    SplitLat = lat(sub2ind(size(lat),valfloat,valtime));
    num = length(find(SplitLat>5));
    if (num>0)
        NECCFlts((cnt+1):(cnt+num)) = valfloat(SplitLat > 5);
    end
    cnt = length(NECCFlts);
end

NECCFlts = sort(NECCFlts)';
identicals = [false;diff(NECCFlts)==0];
NECCFlts = NECCFlts(~identicals);
condNECC = zeros(size(condENTFIRST(:,1)));
condNECC(NECCFlts) = 1;
condNECC = repmat(condNECC,[1 length(condENTFIRST(1,:))]);
condEUC = ~condNECC;

[timeENTEUC,FloatENTEUC] = find(condENTFIRST' & condEUC');
VorENTEUC = Vortex(FloatENTEUC,1); %Vortices of these floats.
[timeENTNECC,FloatENTNECC] = find(condENTFIRST' & condNECC');
VorENTNECC = Vortex(FloatENTNECC,1); %Vortices of these floats.

%EUC COM:
COMEUC = zeros(2*tLF,5);
COMEUC(:,1) = (-tLF+1):1:tLF;
for t = 1:2*tLF
    time = COMEUC(t,1);
    valtime = timeENTEUC+time;
    valfloat = FloatENTEUC(valtime>0 & valtime<=tLF);
    valtime = valtime(valtime>0 & valtime<=tLF);
    COMEUCRad = SFr(sub2ind(size(SFr),valfloat,valtime));
    COMEUCAng = SFa(sub2ind(size(SFa),valfloat,valtime));
    COMEUC(t,4) = length(find(~isnan(COMEUCRad)));
    COMEUCAng = COMEUCAng(~isnan(COMEUCAng));
    COMEUCRad = COMEUCRad(~isnan(COMEUCRad));
    COMEUCx = COMEUCRad.*cos(COMEUCAng);
    COMEUCy = COMEUCRad.*sin(COMEUCAng);
    COMEUCx = mean(COMEUCx);
    COMEUCy = mean(COMEUCy);
    COMEUC(t,3) = atan2(COMEUCy,COMEUCx);
    COMEUC(t,2) = sqrt(COMEUCx^2+COMEUCy^2);
    COMEUCDep = depth(sub2ind(size(depth),valfloat,valtime));
    COMEUC(t,5) = mean(COMEUCDep(~isnan(COMEUCDep)));
end

%NECC COMNECC:
COMNECC = zeros(2*tLF,5);
COMNECC(:,1) = (-tLF+1):1:tLF;
for t = 1:2*tLF
    time = COMNECC(t,1);
    valtime = timeENTNECC+time;
    valfloat = FloatENTNECC(valtime>0 & valtime<=tLF);
    valtime = valtime(valtime>0 & valtime<=tLF);
    COMNECCRad = SFr(sub2ind(size(SFr),valfloat,valtime));
    COMNECCAng = SFa(sub2ind(size(SFa),valfloat,valtime));
    COMNECC(t,4) = length(find(~isnan(COMNECCRad)));
    COMNECCAng = COMNECCAng(~isnan(COMNECCAng));
    COMNECCRad = COMNECCRad(~isnan(COMNECCRad));
    COMNECCx = COMNECCRad.*cos(COMNECCAng);
    COMNECCy = COMNECCRad.*sin(COMNECCAng);
    COMNECCx = mean(COMNECCx);
    COMNECCy = mean(COMNECCy);
    COMNECC(t,3) = atan2(COMNECCy,COMNECCx);
    COMNECC(t,2) = sqrt(COMNECCx^2+COMNECCy^2);
    COMNECCDep = depth(sub2ind(size(depth),valfloat,valtime));
    COMNECC(t,5) = mean(COMNECCDep(~isnan(COMNECCDep)));
end

%save('basinCOM.mat','COM','COMEUC','COMNECC');

%%Time Shifted PDFs at different
%%times:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GROUPS = cell(NoV,tLF);
NumENTf = zeros(NoV,tLF);
for V = 1:NoV
    for t = 1:tLF
        GROUPS{V,t} = FloatENTf(VorENTf == V & timeENTf == t);
        NumENTf(V,t) = length(GROUPS{V,t});
    end
end
%GROUPS{V,t} = Float numbers that enter vortex V at time t
%correctly.
%NumENTf = number of floats in that group.

'Plotting...'
fig1 = figure;
set(fig1, 'Position', get(0,'Screensize')); % Maximize figure. 
set(fig1,'Color',[1 1 1]);

% $$$ %%%%%Temperature:
% $$$ subplot(3,2,1);
% $$$ % $$$ PDF = FloatPDF([1 2 3 4],[10 33],250,1);
% $$$ % $$$ PDFtempSavedx = PDF(:,1);
% $$$ % $$$ PDFtempSaved = zeros(length(PDF(:,1)),LTSVEC);
% $$$ PDF = FloatPDF([10:33],[10 33],Bins,Filter);
% $$$ PDFx = PDF(:,1);
% $$$ PDFy = 0*PDF(:,2);
% $$$ for TS = 1:length(PDFtimeSHF)
% $$$     PDFy = 0*PDFy;
% $$$     cntPDFs = 0;
% $$$     for V = 1:NoV
% $$$         for t = 1:tLF
% $$$             timeshft = t+PDFtimeSHF(TS);
% $$$             if (NumENTf(V,t)>thFLTS & timeshft>0 & timeshft < tLF)
% $$$             
% $$$             Var = temp(sub2ind(size(temp),GROUPS{V,t},timeshft*ones(length(GROUPS{V,t}),1)));
% $$$             if (length(find(~isnan(Var)))>thFLTS)
% $$$             PDF = FloatPDF(Var(~isnan(Var)),[10 33],Bins,Filter);
% $$$             PDFy = PDFy + PDF(:,2);
% $$$             cntPDFs = cntPDFs+1;
% $$$             end
% $$$             end
% $$$         end
% $$$     end
% $$$     if (cntPDFs > thPDFs)
% $$$         %Normalize:
% $$$         Total = sum((PDFy(1:(end-1))+PDFy(2:end))/2.*(PDFx(2:end)-PDFx(1: ...
% $$$                                                   (end-1))));
% $$$         PDFy = PDFy/Total;
% $$$     plot(PDFx,PDFy,'Color',cc(TS,:),'LineWidth',3);
% $$$     hold on;
% $$$     Mean = sum((PDFy(1:(end-1)).*PDFx(1:(end-1))+PDFy(2:end).*PDFx(2:end))/2.*(PDFx(2:end)-PDFx(1: ...
% $$$                                                   (end-1))));
% $$$     plot([Mean Mean],[0.45 0.5],'Color',cc(TS,:),'LineWidth',2);
% $$$     end
% $$$ end
% $$$ %cb = colorbar;
% $$$ %ylabel(cb,'Days to TIV entrance');
% $$$ colorbar('off');
% $$$ caxis(caxs);
% $$$ colormap(cc);
% $$$ %title('basin TIV entrance Float evolution');
% $$$ xlabel('Temperature (\circ C)');
% $$$ ylabel('Probability Density');
% $$$ xlim([10 30]);
% $$$ ylim([0 0.5]);
% $$$ %set(subplot(3,2,1),'Position',[0.05 0.74 0.4 0.2]);
% $$$ 
% $$$ %%%%%%Density:
% $$$ subplot(3,2,2);
% $$$ % $$$ PDF = FloatPDF([1 2 3 4],[10 33],250,1);
% $$$ % $$$ PDFtempSavedx = PDF(:,1);
% $$$ % $$$ PDFtempSaved = zeros(length(PDF(:,1)),LTSVEC);
% $$$ PDF = FloatPDF([20:28],[20 28],Bins,Filter);
% $$$ PDFx = PDF(:,1);
% $$$ PDFy = 0*PDF(:,2);
% $$$ for TS = 1:length(PDFtimeSHF)
% $$$     PDFy = 0*PDFy;
% $$$     cntPDFs = 0;
% $$$     for V = 1:NoV
% $$$         for t = 1:tLF
% $$$             timeshft = t+PDFtimeSHF(TS);
% $$$             if (NumENTf(V,t)>thFLTS & timeshft>0 & timeshft < tLF)
% $$$             
% $$$             Var = rho(sub2ind(size(rho),GROUPS{V,t},timeshft*ones(length(GROUPS{V,t}),1)));
% $$$             if (length(find(~isnan(Var)))>thFLTS)
% $$$             PDF = FloatPDF(Var(~isnan(Var)),[20 28],Bins,Filter);
% $$$             PDFy = PDFy + PDF(:,2);
% $$$             cntPDFs = cntPDFs+1;
% $$$             end
% $$$             end
% $$$         end
% $$$     end
% $$$     if (cntPDFs > thPDFs)
% $$$         %Normalize:
% $$$         Total = sum((PDFy(1:(end-1))+PDFy(2:end))/2.*(PDFx(2:end)-PDFx(1: ...
% $$$                                                   (end-1))));
% $$$         PDFy = PDFy/Total;
% $$$     plot(PDFx,PDFy,'Color',cc(TS,:),'LineWidth',3);
% $$$     hold on;
% $$$     Mean = sum((PDFy(1:(end-1)).*PDFx(1:(end-1))+PDFy(2:end).*PDFx(2:end))/2.*(PDFx(2:end)-PDFx(1: ...
% $$$                                                   (end-1))));
% $$$     plot([Mean Mean],[1.35 1.5],'Color',cc(TS,:),'LineWidth',2);
% $$$     end
% $$$ end
% $$$ cb = colorbar;
% $$$ ylabel(cb,'Days to TIV entrance');
% $$$ caxis(caxs);
% $$$ colormap(cc);
% $$$ %title('basin TIV entrance Float evolution');
% $$$ xlabel('Density (kg m^{-3})');
% $$$ %ylabel('Probability Density');
% $$$ xlim([20 28]);
% $$$ ylim([0 1.5]);
% $$$ %set(subplot(3,2,1),'Position',[0.05 0.74 0.4 0.2]);
% $$$ 
% $$$ %%%%%%Salinity:
% $$$ subplot(3,2,3);
% $$$ % $$$ PDF = FloatPDF([1 2 3 4],[10 33],250,1);
% $$$ % $$$ PDFtempSavedx = PDF(:,1);
% $$$ % $$$ PDFtempSaved = zeros(length(PDF(:,1)),LTSVEC);
% $$$ PDF = FloatPDF([33:36],[33 36],Bins,Filter);
% $$$ PDFx = PDF(:,1);
% $$$ PDFy = 0*PDF(:,2);
% $$$ for TS = 1:length(PDFtimeSHF)
% $$$     PDFy = 0*PDFy;
% $$$     cntPDFs = 0;
% $$$     for V = 1:NoV
% $$$         for t = 1:tLF
% $$$             timeshft = t+PDFtimeSHF(TS);
% $$$             if (NumENTf(V,t)>thFLTS & timeshft>0 & timeshft < tLF)
% $$$             
% $$$             Var = salt(sub2ind(size(salt),GROUPS{V,t},timeshft*ones(length(GROUPS{V,t}),1)));
% $$$             if (length(find(~isnan(Var)))>thFLTS)
% $$$             PDF = FloatPDF(Var(~isnan(Var)),[33 36],Bins,Filter);
% $$$             PDFy = PDFy + PDF(:,2);
% $$$             cntPDFs = cntPDFs+1;
% $$$             end
% $$$             end
% $$$         end
% $$$     end
% $$$     if (cntPDFs > thPDFs)
% $$$         %Normalize:
% $$$         Total = sum((PDFy(1:(end-1))+PDFy(2:end))/2.*(PDFx(2:end)-PDFx(1: ...
% $$$                                                   (end-1))));
% $$$         PDFy = PDFy/Total;
% $$$     plot(PDFx,PDFy,'Color',cc(TS,:),'LineWidth',3);
% $$$     hold on;
% $$$     Mean = sum((PDFy(1:(end-1)).*PDFx(1:(end-1))+PDFy(2:end).*PDFx(2:end))/2.*(PDFx(2:end)-PDFx(1: ...
% $$$                                                   (end-1))));
% $$$     plot([Mean Mean],[5.4 6],'Color',cc(TS,:),'LineWidth',2);
% $$$     end
% $$$ end
% $$$ colorbar('off');
% $$$ %cb = colorbar;
% $$$ %ylabel(cb,'Days to TIV entrance');
% $$$ caxis(caxs);
% $$$ colormap(cc);
% $$$ %title('basin TIV entrance Float evolution');
% $$$ xlabel('Salinity (psu)');
% $$$ ylabel('Probability Density');
% $$$ xlim([33.5 35.5]);
% $$$ ylim([0 6]);%
% $$$             %set(subplot(3,2,1),'Position',[0.05 0.74 0.4 0.2]);
% $$$ 
% $$$ %%%%%%%%%Latitude:
% $$$ subplot(3,2,4);
% $$$ % $$$ PDF = FloatPDF([1 2 3 4],[10 33],250,1);
% $$$ % $$$ PDFtempSavedx = PDF(:,1);
% $$$ % $$$ PDFtempSaved = zeros(length(PDF(:,1)),LTSVEC);
% $$$ PDF = FloatPDF([-5:15],[-5 15],Bins,Filter);
% $$$ PDFx = PDF(:,1);
% $$$ PDFy = 0*PDF(:,2);
% $$$ for TS = 1:length(PDFtimeSHF)
% $$$     PDFy = 0*PDFy;
% $$$     cntPDFs = 0;
% $$$     for V = 1:NoV
% $$$         for t = 1:tLF
% $$$             timeshft = t+PDFtimeSHF(TS);
% $$$             if (NumENTf(V,t)>thFLTS & timeshft>0 & timeshft < tLF)
% $$$             
% $$$             Var = lat(sub2ind(size(lat),GROUPS{V,t},timeshft*ones(length(GROUPS{V,t}),1)));
% $$$             if (length(find(~isnan(Var)))>thFLTS)
% $$$             PDF = FloatPDF(Var(~isnan(Var)),[-5 15],Bins,Filter);
% $$$             PDFy = PDFy + PDF(:,2);
% $$$             cntPDFs = cntPDFs+1;
% $$$             end
% $$$             end
% $$$         end
% $$$     end
% $$$     if (cntPDFs > thPDFs)
% $$$         %Normalize:
% $$$         Total = sum((PDFy(1:(end-1))+PDFy(2:end))/2.*(PDFx(2:end)-PDFx(1: ...
% $$$                                                   (end-1))));
% $$$         PDFy = PDFy/Total;
% $$$     plot(PDFx,PDFy,'Color',cc(TS,:),'LineWidth',3);
% $$$     hold on;
% $$$     Mean = sum((PDFy(1:(end-1)).*PDFx(1:(end-1))+PDFy(2:end).*PDFx(2:end))/2.*(PDFx(2:end)-PDFx(1: ...
% $$$                                                   (end-1))));
% $$$     plot([Mean Mean],[0.45 0.5],'Color',cc(TS,:),'LineWidth',2);
% $$$     end
% $$$ end
% $$$ cb = colorbar;
% $$$ ylabel(cb,'Days to TIV entrance');
% $$$ caxis(caxs);
% $$$ colormap(cc);
% $$$ %title('basin TIV entrance Float evolution');
% $$$ xlabel('Latitude (\circ N)');
% $$$ %ylabel('Probability Density');
% $$$ xlim([-5 15]);
% $$$ ylim([0 0.5]);
% $$$ %set(subplot(3,2,1),'Position',[0.05 0.74 0.4 0.2]);
% $$$ 
% $$$ %%%%%%%PV:
% $$$ subplot(3,2,5);
% $$$ % $$$ PDF = FloatPDF([1 2 3 4],[10 33],250,1);
% $$$ % $$$ PDFtempSavedx = PDF(:,1);
% $$$ % $$$ PDFtempSaved = zeros(length(PDF(:,1)),LTSVEC);
% $$$ PDF = FloatPDF([-0.5:0.1:2]*1e-8,[-0.5 2]*1e-8,Bins,Filter);
% $$$ PDFx = PDF(:,1);
% $$$ PDFy = 0*PDF(:,2);
% $$$ for TS = 1:length(PDFtimeSHF)
% $$$     PDFy = 0*PDFy;
% $$$     cntPDFs = 0;
% $$$     for V = 1:NoV
% $$$         for t = 1:tLF
% $$$             timeshft = t+PDFtimeSHF(TS);
% $$$             if (NumENTf(V,t)>thFLTS & timeshft>0 & timeshft < tLF)
% $$$             
% $$$             Var = PV(sub2ind(size(PV),GROUPS{V,t},timeshft*ones(length(GROUPS{V,t}),1)));
% $$$             if (length(find(~isnan(Var)))>thFLTS)
% $$$             PDF = FloatPDF(Var(~isnan(Var)),[-0.5 2]*1e-8,Bins,Filter);
% $$$             PDFy = PDFy + PDF(:,2);
% $$$             cntPDFs = cntPDFs+1;
% $$$             end
% $$$             end
% $$$         end
% $$$     end
% $$$     if (cntPDFs > thPDFs)
% $$$         %Normalize:
% $$$         Total = sum((PDFy(1:(end-1))+PDFy(2:end))/2.*(PDFx(2:end)-PDFx(1: ...
% $$$                                                   (end-1))));
% $$$         PDFy = PDFy/Total;
% $$$     plot(PDFx,PDFy,'Color',cc(TS,:),'LineWidth',3);
% $$$     hold on;
% $$$     Mean = sum((PDFy(1:(end-1)).*PDFx(1:(end-1))+PDFy(2:end).*PDFx(2:end))/2.*(PDFx(2:end)-PDFx(1: ...
% $$$                                                   (end-1))));
% $$$     plot([Mean Mean],[9e8 10e8],'Color',cc(TS,:),'LineWidth',2);
% $$$     end
% $$$ end
% $$$ %cb = colorbar;
% $$$ colorbar('off');
% $$$ %ylabel(cb,'Days to TIV entrance');
% $$$ caxis(caxs);
% $$$ colormap(cc);
% $$$ %title('basin TIV entrance Float evolution');
% $$$ xlabel('Potential Vorticity (s^{-3})');
% $$$ ylabel('Probability Density');
% $$$ xlim([-0.3 1]*1e-8);
% $$$ ylim([0 10]*1e8);
% $$$ %set(subplot(3,2,1),'Position',[0.05 0.74 0.4 0.2]);
% $$$ 
% $$$ %%%%%%PV_H:
% $$$ subplot(3,2,6);
% $$$ % $$$ PDF = FloatPDF([1 2 3 4],[10 33],250,1);
% $$$ % $$$ PDFtempSavedx = PDF(:,1);
% $$$ % $$$ PDFtempSaved = zeros(length(PDF(:,1)),LTSVEC);
% $$$ PDF = FloatPDF([-0.6:0.1:0.2]*1e-8,[-0.6 0.2]*1e-8,Bins,Filter);
% $$$ PDFx = PDF(:,1);
% $$$ PDFy = 0*PDF(:,2);
% $$$ for TS = 1:length(PDFtimeSHF)
% $$$     PDFy = 0*PDFy;
% $$$     cntPDFs = 0;
% $$$     for V = 1:NoV
% $$$         for t = 1:tLF
% $$$             timeshft = t+PDFtimeSHF(TS);
% $$$             if (NumENTf(V,t)>thFLTS & timeshft>0 & timeshft < tLF)
% $$$             
% $$$             Var = PVh(sub2ind(size(PVh),GROUPS{V,t},timeshft*ones(length(GROUPS{V,t}),1)));
% $$$             if (length(find(~isnan(Var)))>thFLTS)
% $$$             PDF = FloatPDF(Var(~isnan(Var)),[-0.6 0.2]*1e-8,Bins,Filter);
% $$$             PDFy = PDFy + PDF(:,2);
% $$$             cntPDFs = cntPDFs+1;
% $$$             end
% $$$             end
% $$$         end
% $$$     end
% $$$     if (cntPDFs > thPDFs)
% $$$         %Normalize:
% $$$         Total = sum((PDFy(1:(end-1))+PDFy(2:end))/2.*(PDFx(2:end)-PDFx(1: ...
% $$$                                                   (end-1))));
% $$$         PDFy = PDFy/Total;
% $$$     plot(PDFx,PDFy,'Color',cc(TS,:),'LineWidth',3);
% $$$     hold on;
% $$$     Mean = sum((PDFy(1:(end-1)).*PDFx(1:(end-1))+PDFy(2:end).*PDFx(2:end))/2.*(PDFx(2:end)-PDFx(1: ...
% $$$                                                   (end-1))));
% $$$     plot([Mean Mean],[54e8 60e8],'Color',cc(TS,:),'LineWidth',2);
% $$$     end
% $$$ end
% $$$ cb = colorbar;
% $$$ ylabel(cb,'Days to TIV entrance');
% $$$ caxis(caxs);
% $$$ colormap(cc);
% $$$ %title('basin TIV entrance Float evolution');
% $$$ xlabel('Baroclinic Potential Vorticity (s^{-3})');
% $$$ %ylabel('Probability Density');
% $$$ xlim([-0.3 0.1]*1e-8);
% $$$ ylim([0 60]*1e8);
% $$$ %set(subplot(3,2,1),'Position',[0.05 0.74 0.4 0.2]);
% $$$ 
% $$$ set(subplot(3,2,1),'Position',[0.05 0.74 0.4 0.2]);
% $$$ set(subplot(3,2,3),'Position',[0.05 0.44 0.4 0.2]);
% $$$ set(subplot(3,2,5),'Position',[0.05 0.13 0.4 0.2]);
% $$$ set(subplot(3,2,2),'Position',[0.5 0.74 0.4 0.2]);
% $$$ set(subplot(3,2,4),'Position',[0.5 0.44 0.4 0.2]);
% $$$ set(subplot(3,2,6),'Position',[0.5 0.13 0.4 0.2]);
% $$$ 
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ %%%%%%OTHER:
% $$$ 'Plotting...'
% $$$ fig1 = figure;
% $$$ set(fig1, 'Position', get(0,'Screensize')); % Maximize figure. 
% $$$ set(fig1,'Color',[1 1 1]);
% $$$ subplot(3,1,1);
% $$$ PDF = FloatPDF([-pi:1:pi],[-pi pi],Bins,Filter);
% $$$ PDFx = PDF(:,1);
% $$$ PDFy = 0*PDF(:,2);
% $$$ for TS = 1:length(PDFtimeSHF)
% $$$     PDFy = 0*PDFy;
% $$$     cntPDFs = 0;
% $$$     for V = 1:NoV
% $$$         for t = 1:tLF
% $$$             timeshft = t+PDFtimeSHF(TS);
% $$$             if (NumENTf(V,t)>thFLTS & timeshft>0 & timeshft < tLF)
% $$$             
% $$$             Var = SFa(sub2ind(size(SFa),GROUPS{V,t},timeshft*ones(length(GROUPS{V,t}),1)));
% $$$             if (length(find(~isnan(Var)))>thFLTS)
% $$$             PDF = FloatPDF(Var(~isnan(Var)),[-pi pi],Bins,Filter);
% $$$             PDFy = PDFy + PDF(:,2);
% $$$             cntPDFs = cntPDFs+1;
% $$$             end
% $$$             end
% $$$         end
% $$$     end
% $$$     if (cntPDFs > thPDFs)
% $$$         %Normalize:
% $$$         Total = sum((PDFy(1:(end-1))+PDFy(2:end))/2.*(PDFx(2:end)-PDFx(1: ...
% $$$                                                   (end-1))));
% $$$         PDFy = PDFy/Total;
% $$$     plot(PDFx,PDFy,'Color',cc(TS,:),'LineWidth',3);
% $$$     hold on;
% $$$     Mean = sum((PDFy(1:(end-1)).*PDFx(1:(end-1))+PDFy(2:end).*PDFx(2:end))/2.*(PDFx(2:end)-PDFx(1: ...
% $$$                                                   (end-1))));
% $$$     %    plot([Mean Mean],[63e8 70e8],'Color',cc(TS,:),'LineWidth',2);
% $$$     end
% $$$ end
% $$$ cb = colorbar;
% $$$ ylabel(cb,'Days to TIV entrance');
% $$$ caxis(caxs);
% $$$ colormap(cc);
% $$$ %title('basin TIV entrance Float evolution');
% $$$ xlabel('Angle (radians)');
% $$$ %ylabel('Probability Density');
% $$$ %xlim([-0.3 0.1]*1e-8);
% $$$ %ylim([0 70]*1e8);
% $$$ %set(subplot(3,2,1),'Position',[0.05 0.74 0.4 0.2]);
% $$$ 
% $$$ subplot(3,1,2);
% $$$ PDF = FloatPDF([-300:1:0],[-300 0],Bins,Filter);
% $$$ PDFx = PDF(:,1);
% $$$ PDFy = 0*PDF(:,2);
% $$$ for TS = 1:length(PDFtimeSHF)
% $$$     PDFy = 0*PDFy;
% $$$     cntPDFs = 0;
% $$$     for V = 1:NoV
% $$$         for t = 1:tLF
% $$$             timeshft = t+PDFtimeSHF(TS);
% $$$             if (NumENTf(V,t)>thFLTS & timeshft>0 & timeshft < tLF)
% $$$             
% $$$             Var = depth(sub2ind(size(depth),GROUPS{V,t},timeshft*ones(length(GROUPS{V,t}),1)));
% $$$             if (length(find(~isnan(Var)))>thFLTS)
% $$$             PDF = FloatPDF(Var(~isnan(Var)),[-300 0],Bins,Filter);
% $$$             PDFy = PDFy + PDF(:,2);
% $$$             cntPDFs = cntPDFs+1;
% $$$             end
% $$$             end
% $$$         end
% $$$     end
% $$$     if (cntPDFs > thPDFs)
% $$$         %Normalize:
% $$$         Total = sum((PDFy(1:(end-1))+PDFy(2:end))/2.*(PDFx(2:end)-PDFx(1: ...
% $$$                                                   (end-1))));
% $$$         PDFy = PDFy/Total;
% $$$     plot(PDFx,PDFy,'Color',cc(TS,:),'LineWidth',3);
% $$$     hold on;
% $$$     Mean = sum((PDFy(1:(end-1)).*PDFx(1:(end-1))+PDFy(2:end).*PDFx(2:end))/2.*(PDFx(2:end)-PDFx(1: ...
% $$$                                                   (end-1))));
% $$$     %    plot([Mean Mean],[63e8 70e8],'Color',cc(TS,:),'LineWidth',2);
% $$$     end
% $$$ end
% $$$ cb = colorbar;
% $$$ ylabel(cb,'Days to TIV entrance');
% $$$ caxis(caxs);
% $$$ colormap(cc);
% $$$ %title('basin TIV entrance Float evolution');
% $$$ xlabel('Depth (m)');
% $$$ %ylabel('Probability Density');
% $$$ %xlim([-0.3 0.1]*1e-8);
% $$$ %ylim([0 70]*1e8);
% $$$ %set(subplot(3,2,1),'Position',[0.05 0.74 0.4 0.2]);
% $$$ 
% $$$ subplot(3,1,3);
% $$$ PDF = FloatPDF([1e5:10:1e6],[1e2 5e6],Bins,Filter);
% $$$ PDFx = PDF(:,1);
% $$$ PDFy = 0*PDF(:,2);
% $$$ for TS = 1:length(PDFtimeSHF)
% $$$     PDFy = 0*PDFy;
% $$$     cntPDFs = 0;
% $$$     for V = 1:NoV
% $$$         for t = 1:tLF
% $$$             timeshft = t+PDFtimeSHF(TS);
% $$$             if (NumENTf(V,t)>thFLTS & timeshft>0 & timeshft < tLF)
% $$$             
% $$$             Var = SFr(sub2ind(size(SFr),GROUPS{V,t},timeshft*ones(length(GROUPS{V,t}),1)));
% $$$             if (length(find(~isnan(Var)))>thFLTS)
% $$$             PDF = FloatPDF(Var(~isnan(Var)),[1e2 5e6],Bins,Filter);
% $$$             PDFy = PDFy + PDF(:,2);
% $$$             cntPDFs = cntPDFs+1;
% $$$             end
% $$$             end
% $$$         end
% $$$     end
% $$$     if (cntPDFs > thPDFs)
% $$$         %Normalize:
% $$$         Total = sum((PDFy(1:(end-1))+PDFy(2:end))/2.*(PDFx(2:end)-PDFx(1: ...
% $$$                                                   (end-1))));
% $$$         PDFy = PDFy/Total;
% $$$     plot(PDFx,PDFy,'Color',cc(TS,:),'LineWidth',3);
% $$$     hold on;
% $$$     Mean = sum((PDFy(1:(end-1)).*PDFx(1:(end-1))+PDFy(2:end).*PDFx(2:end))/2.*(PDFx(2:end)-PDFx(1: ...
% $$$                                                   (end-1))));
% $$$     %    plot([Mean Mean],[63e8 70e8],'Color',cc(TS,:),'LineWidth',2);
% $$$     end
% $$$ end
% $$$ cb = colorbar;
% $$$ ylabel(cb,'Days to TIV entrance');
% $$$ caxis(caxs);
% $$$ colormap(cc);
% $$$ %title('basin TIV entrance Float evolution');
% $$$ xlabel('Radius (m)');
% $$$ %ylabel('Probability Density');
% $$$ %xlim([-0.3 0.1]*1e-8);
% $$$ %ylim([0 70]*1e8);
% $$$ %set(subplot(3,2,1),'Position',[0.05 0.74 0.4 0.2]);
% $$$ 
% $$$ %Richardson Number:
% $$$ subplot(3,1,3);
% $$$ PDF = FloatPDF([-3:1:3],[-3 3],Bins,Filter);
% $$$ PDFx = PDF(:,1);
% $$$ PDFy = 0*PDF(:,2);
% $$$ for TS = 1:length(PDFtimeSHF)
% $$$     PDFy = 0*PDFy;
% $$$     cntPDFs = 0;
% $$$     for V = 1:NoV
% $$$         for t = 1:tLF
% $$$             timeshft = t+PDFtimeSHF(TS);
% $$$             if (NumENTf(V,t)>thFLTS & timeshft>0 & timeshft < tLF)
% $$$             
% $$$             Var = Ri(sub2ind(size(Ri),GROUPS{V,t},timeshft* ...
% $$$                              ones(length(GROUPS{V,t}),1)));
% $$$             Var(Var<0) = 0;
% $$$             Var = log10(Var);
% $$$             if (length(find(~isnan(Var)))>thFLTS)
% $$$             PDF = FloatPDF(Var(~isnan(Var)),[-3 3],Bins,Filter);
% $$$             PDFy = PDFy + PDF(:,2);
% $$$             cntPDFs = cntPDFs+1;
% $$$             end
% $$$             end
% $$$         end
% $$$     end
% $$$     if (cntPDFs > thPDFs)
% $$$         %Normalize:
% $$$         Total = sum((PDFy(1:(end-1))+PDFy(2:end))/2.*(PDFx(2:end)-PDFx(1: ...
% $$$                                                   (end-1))));
% $$$         PDFy = PDFy/Total;
% $$$     plot(PDFx,PDFy,'Color',cc(TS,:),'LineWidth',3);
% $$$     hold on;
% $$$     Mean = sum((PDFy(1:(end-1)).*PDFx(1:(end-1))+PDFy(2:end).*PDFx(2:end))/2.*(PDFx(2:end)-PDFx(1: ...
% $$$                                                   (end-1))));
% $$$     %    plot([Mean Mean],[63e8 70e8],'Color',cc(TS,:),'LineWidth',2);
% $$$     end
% $$$ end
% $$$ cb = colorbar;
% $$$ ylabel(cb,'Days to TIV entrance');
% $$$ caxis(caxs);
% $$$ colormap(cc);
% $$$ %title('basin TIV entrance Float evolution');
% $$$ xlabel('Richardson Number');
% $$$ %ylabel('Probability Density');
% $$$ %xlim([-0.3 0.1]*1e-8);
% $$$ %ylim([0 70]*1e8);
% $$$ %set(subplot(3,2,1),'Position',[0.05 0.74 0.4 0.2]);


% $$$ %Frontogenesis function:
% $$$ %figure;
% $$$ PDF = FloatPDF([-26:1:-15],[-26 -15],Bins,Filter);
% $$$ PDFx = PDF(:,1);
% $$$ PDFy = 0*PDF(:,2);
% $$$ for TS = 1:length(PDFtimeSHF)
% $$$     PDFy = 0*PDFy;
% $$$     cntPDFs = 0;
% $$$     for V = 1:NoV
% $$$         for t = 1:tLF
% $$$             timeshft = t+PDFtimeSHF(TS);
% $$$             if (NumENTf(V,t)>thFLTS & timeshft>0 & timeshft < tLF)
% $$$             
% $$$             Var = FGF(sub2ind(size(FGF),GROUPS{V,t},timeshft* ...
% $$$                              ones(length(GROUPS{V,t}),1)));
% $$$             if (length(find(~isnan(Var)))>thFLTS)
% $$$             PDF = FloatPDF(Var(~isnan(Var)),[-26 -15],Bins,Filter);
% $$$             PDFy = PDFy + PDF(:,2);
% $$$             cntPDFs = cntPDFs+1;
% $$$             end
% $$$             end
% $$$         end
% $$$     end
% $$$     if (cntPDFs > thPDFs)
% $$$         %Normalize:
% $$$         Total = sum((PDFy(1:(end-1))+PDFy(2:end))/2.*(PDFx(2:end)-PDFx(1: ...
% $$$                                                   (end-1))));
% $$$         PDFy = PDFy/Total;
% $$$     plot(PDFx,PDFy,'Color',cc(TS,:),'LineWidth',3);
% $$$     hold on;
% $$$     Mean = sum((PDFy(1:(end-1)).*PDFx(1:(end-1))+PDFy(2:end).*PDFx(2:end))/2.*(PDFx(2:end)-PDFx(1: ...
% $$$                                                   (end-1))));
% $$$     plot([Mean Mean],[0.63 0.7],'Color',cc(TS,:),'LineWidth',2);
% $$$     end
% $$$ end
% $$$ cb = colorbar;
% $$$ ylabel(cb,'Days to TIV entrance');
% $$$ caxis(caxs);
% $$$ colormap(cc);
% $$$ %title('basin TIV entrance Float evolution');
% $$$ xlabel('Richardson Number');
% $$$ %ylabel('Probability Density');
% $$$ %xlim([-0.3 0.1]*1e-8);
% $$$ %ylim([0 70]*1e8);
% $$$ %set(subplot(3,2,1),'Position',[0.05 0.74 0.4 0.2]);

%Frontogenesis function:
%figure;
PDF = FloatPDF([-1e-8:1e-7:0.25e-6],[-1e-8 0.25e-6],Bins,Filter);
PDFx = PDF(:,1);
PDFy = 0*PDF(:,2);
MeanSave = zeros(length(PDFtimeSHF),1);
SDSave = zeros(length(PDFtimeSHF),1);
OM3Save = zeros(length(PDFtimeSHF),1);
for TS = 1:length(PDFtimeSHF)
    PDFy = 0*PDFy;
    cntPDFs = 0;
    for V = 1:NoV
        for t = 1:tLF
            timeshft = t+PDFtimeSHF(TS);
            if (NumENTf(V,t)>thFLTS & timeshft>0 & timeshft < tLF)
            Var = Gradhb(sub2ind(size(Gradhb),GROUPS{V,t},timeshft* ...
                             ones(length(GROUPS{V,t}),1)));
            if (length(find(~isnan(Var)))>thFLTS)
            PDF = FloatPDF(Var(~isnan(Var)),[-1e-8 0.25e-6],Bins,Filter);
            PDFy = PDFy + PDF(:,2);
            cntPDFs = cntPDFs+1;
            end
            end
        end
    end
    if (cntPDFs > thPDFs)
        %Normalize:
        Total = sum((PDFy(1:(end-1))+PDFy(2:end))/2.*(PDFx(2:end)-PDFx(1: ...
                                                  (end-1))));
        PDFy = PDFy/Total;
    plot(PDFx,PDFy,'Color',cc(TS,:),'LineWidth',3);
    hold on;
    Mean = sum((PDFy(1:(end-1)).*PDFx(1:(end-1))+PDFy(2:end).*PDFx(2:end))/2.*(PDFx(2:end)-PDFx(1: ...
                                                  (end-1))));
    SD = sqrt(sum((PDFy(1:(end-1)).*PDFx(1:(end-1)).^2+PDFy(2:end).*PDFx(2:end).^2)/2.*(PDFx(2:end)-PDFx(1: ...
                                                  (end-1)))));
    OM3 = (sum((PDFy(1:(end-1)).*PDFx(1:(end-1)).^3+PDFy(2:end).*PDFx(2:end).^3)/2.*(PDFx(2:end)-PDFx(1: ...
                                                  (end-1)))))^(1/3);
    plot([Mean Mean],[1.8 2]*1e7,'Color',cc(TS,:),'LineWidth',2);
    MeanSave(TS) = Mean;
    SDSave(TS) = SD;
    OM3Save(TS) = OM3;
    end
end
cb = colorbar;
ylabel(cb,'Days to TIV entrance');
caxis(caxs);
colormap(cc);
%title('basin TIV entrance Float evolution');
xlabel('$|\nabla_h b|$','Interpreter','latex','FontSize',25);
ylabel('Probability Density');
%xlim([-0.3 0.1]*1e-8);
%ylim([0 70]*1e8);
%set(subplot(3,2,1),'Position',[0.05 0.74 0.4 0.2]);
axis([-0.3e-8 2e-7 0 2e7]);
axes('Position',[0.4 0.55 0.35 0.25]);
plot(PDFtimeSHF,MeanSave,'-k','LineWidth',3);
hold on;
plot(PDFtimeSHF,MeanSave+SDSave,'--k','LineWidth',2);
plot(PDFtimeSHF,MeanSave-SDSave,'--k','LineWidth',2);
set(gca,'FontSize',15);
axis([-100 0 0.3e-7 0.6e-7]);
set(gca,'ytick',[0.3 0.4 0.5 0.6]*1e-7);
set(gca,'xtick',-100:20:0);
xlabel('Days to TIV Entrance');
ylabel('Mean');

figure;
subplot(3,1,1);
plot(PDFtimeSHF,MeanSave,'-k','LineWidth',3);
ylabel('Mean');
subplot(3,1,2);
plot(PDFtimeSHF,SDSave,'-k','LineWidth',3);
ylabel('SD');
subplot(3,1,3);
plot(PDFtimeSHF,OM3Save,'-k','LineWidth',3);
ylabel('3rd order');
