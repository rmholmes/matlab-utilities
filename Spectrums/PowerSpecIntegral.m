

%calculate a power spectrum integral over different frequency
%bands:

load filenames;
fname = BASF_his;

time = ncread(fname,'ocean_time')/86400;
lon = ncread(fname,'lon_rho');
lat = ncread(fname,'lat_rho');
[tmp loni] = min(abs(lon(:,1)+140));
[tmp lati] = min(abs(lat(1,:)));

signal = squeeze(GetVar(fname,0,{'v'},{loni*[1 1],lati*[1 1],[50 ...
                    50],0}));
tL = length(signal);
%Add some bad NaNs;
nNaNs = 0;
for i=1:nNaNs
    signal(round(rand(1)*tL)) = NaN;
end

%Replace with zeros?:
signal = Repl(signal,NaN,0);


signal = signal-nanmean(signal);

band1 = [0 1/75];
band2 = [1/75 1/12];
band3 = [1/12 1/2];
% $$$ band1 = [0 1/60];
% $$$ band2 = [1/60 1/15];
% $$$ band3 = [1/15 1/2];

%fake signal:
% $$$ signal = sin(2*pi*time/20)+0.5*sin(2*pi*time/10)+2*sin(2*pi* ...
% $$$                                                   time/80);
% $$$ signal = signal-mean(signal);
% $$$ band1 = [1/60 1/15];
% $$$ band2 = [0 1/60];
% $$$ band3 = [1/15 1/2];

%Compare AmpSpec and periodogram:
[freq,amp] = AmpSpec(time,signal);
% [power spec, freq] = periodogram(signal,window,nfft(use
% default),fs,'power')
[pxx,f] = periodogram(signal,hamming(length(signal)),[],1,'power');
ampp = sqrt(2*pxx);

subplot(2,1,1);
plot(time,signal);
subplot(2,1,2);
plot(f,ampp);
hold on;
plot(band1(1)*[1 1],[0 max(ampp)],'-b','LineWidth',2);
plot(band1(2)*[1 1],[0 max(ampp)],'-b','LineWidth',2);
plot(band2(1)*[1 1],[0 max(ampp)],'-b','LineWidth',2);
plot(band2(2)*[1 1],[0 max(ampp)],'-b','LineWidth',2);
plot(band3(1)*[1 1],[0 max(ampp)],'-b','LineWidth',2);
plot(band3(2)*[1 1],[0 max(ampp)],'-b','LineWidth',2);

%Do some bandpass filtering:
bandpower(signal) %This is the total variance in the signal.
std(signal,1)^2
norm(signal,2)^2/length(signal)

%Using the periodogram method:
[pxx,f] = periodogram(signal,hamming(length(signal)),[],1,'psd');
p = bandpower(pxx,f,'psd') %This is another way of calculating the
                           %total variance in the signal.
VAR1 = bandpower(pxx,f,band1,'psd');
VAR2 = bandpower(pxx,f,band2,'psd');
VAR3 = bandpower(pxx,f,band3,'psd');
TVAR = VAR1+VAR2+VAR3
TVARp = p

AMP1 = sqrt(2*VAR1)
AMP2 = sqrt(2*VAR2)
AMP3 = sqrt(2*VAR3)
TAMP = AMP1+AMP2+AMP3
TAMPp = sqrt(bandpower(signal)*2)

RMS1 = sqrt(VAR1)
RMS2 = sqrt(VAR2)
RMS3 = sqrt(VAR3)
TRMS = RMS1+RMS2+RMS3
TRMSp = std(signal,1)

%The fact that these aren't equal -> Is this covariance?

%Compare to filtering method:
subplot(2,1,1);
plot(time,signal);
hold on;
smooth = flowpass(time,signal,1/15,0);
plot(time,smooth,'-b');
plot(time,smooth-flowpass(time,smooth,1/60,0),'-r');
