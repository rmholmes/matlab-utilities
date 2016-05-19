function [freq,amplitude] = AmpSpec(time,signal)
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%
% This function returns the single sided amplitude spectrum of the
% signal on the time.
%
% signal = the signal field, assuming that time is the last
% dimension.
%
% time  = time vector in days.
%
% Outputs: 
%
% freq = frequency vector in cycles per day.
%
% amplitude = amplitude spectrum for each part in signal.
%
%---------------------------------------------------------------------
%
% Dependencies; none. 
%
% Ryans ROMS Matlab and netcdf Utilities 25/10/13
%
%---------------------------------------------------------------------
%---------------------------------------------------------------------

N = length(time); %Number of samples
Fs = 1/(time(2)-time(1));
T = N/Fs;
NFFT = 2^nextpow2(N); % Next power of 2 from length of y

%determine which dimension to do fourier along:
Size = size(signal);
if (length(Size) == 4)
    dim = 4;
    ham = repmat(permute(hamming(N),[2 3 4 1]),[Size(1) Size(2) ...
                        Size(3) 1]);
    rvec = [1 1 1 N];
elseif (length(Size) == 3)
    dim = 3;
    ham = repmat(permute(hamming(N),[2 3 1]),[Size(1) Size(2) 1]);
    rvec = [1 1 N];
elseif (length(Size) == 2 && Size(2)>1)
    dim = 2;
    ham = repmat(permute(hamming(N),[2 1]),[Size(1) 1]);
    rvec = [1 N];
else
    dim = 1;
    ham = hamming(N);
    rvec = [N 1];
end

amplitude = fft((signal-repmat(mean(signal,dim),rvec)).*ham,NFFT, ...
                dim)/N;
%Amplitude correction for hamming window (see
%http://www.diracdelta.co.uk/science/source/h/a/hamming%20window/source.html#.UtnHDvjTnCI)
amplitude = amplitude*1.855;

if (dim == 1)
    amplitude = 2*abs(amplitude(1:NFFT/2+1));
elseif (dim == 2)
    amplitude = 2*abs(amplitude(:,1:NFFT/2+1));
elseif (dim == 3)
    amplitude = 2*abs(amplitude(:,:,1:NFFT/2+1));
elseif (dim == 4)
    amplitude = 2*abs(amplitude(:,:,:,1:NFFT/2+1));
end

freq = Fs/2*linspace(0,1,NFFT/2+1);
end