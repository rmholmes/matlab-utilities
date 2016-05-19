function Ffilt = flowpass(T,F,cut,wind)
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%
% This function performs a low-pass filter of a vector/matrix
% etc. with cutoff frequency cut using a Hanning window or a square
% window in Fourier space.
%
% INPUTS:
%
% T = time x-axis.
% F = signal (with the time dimension being the last dimension with
% length > 1).
% cut = cut off frequency (1/day)
% wind = 0 (hanning) 
%      = 1 (square - sets fourier components above cut to zero)
%
%---------------------------------------------------------------------
%
% Dependencies; none. 
%
% Ryans ROMS Matlab and netcdf Utilities 17/7/13
%
%---------------------------------------------------------------------
%---------------------------------------------------------------------
[d1 d2 d3 d4] = size(F);

%Figure out which dimension is the time dimension, and set the
%permutation vectors:
if (d4>1)
    tdim = 4;
    tL = d4;
    Tvec = [d1 d2 d3 1];
    pvec = [2 3 4 1];
elseif (d3>1)
    tdim = 3;
    tL = d3;
    Tvec = [d1 d2 1 d4];
    pvec = [2 3 1 4];
elseif (d2>1)
    tdim = 2;
    tL = d2;
    Tvec = [d1 1 d3 d4];
    pvec = [2 1 3 4];
elseif (d1>1)
    tdim = 1;
    tL = d1;
    Tvec = [1 d2 d3 d4];
    pvec = [1 2 3 4];
else
   warning('Bad input to flowpass.m...');
end
tL2pow = 2^nextpow2(tL);
sampf = 1/(T(end)-T(end-1));

%obtain fourier space F:
Ffour = fft(F,tL2pow,tdim)/tL;
%Ffour = 2*abs(Ffour(1:tL2pow/2+1));
Tfour = zeros(tL2pow,1);
Tfour(1:(tL2pow/2+1)) = sampf/2*linspace(0,1,tL2pow/2+1);
Tfour((tL2pow/2+2):tL2pow) = flipud(Tfour(2:(tL2pow/2)));

%Filter:
if (wind == 1)
    Ffour(repmat(permute(Tfour,pvec)>cut,Tvec)) = 0;
elseif (wind == 0)
    Ffour = fftshift(Ffour,tdim);
    Tfour = fftshift(Tfour);
    Hn = zeros(size(Tfour));
    wL = length(find(Tfour<=cut));
    Hn(Tfour<=cut) = 0.5*(1-cos(2*pi* (0:1:(wL- 1))'/(wL- 1)));
    Ffour = Ffour.*repmat(permute(Hn,pvec),Tvec);
    Ffour = ifftshift(Ffour,tdim);
end
%Inverse transform:
Ffilt = ifft(Ffour,tL2pow,tdim)*tL;
if (tdim==1)
Ffilt = real(Ffilt(1:tL,:,:,:));
elseif (tdim==2)
Ffilt = real(Ffilt(:,1:tL,:,:));
elseif (tdim==3)
Ffilt = real(Ffilt(:,:,1:tL,:));
elseif (tdim==4)
Ffilt = real(Ffilt(:,:,:,1:tL));
end
end