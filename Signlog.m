function VarOUT = Signlog(Var,minMag,maxMag)
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%
% This function maps the values Var which take both +ve and -ve values
% to a log scale between -1 and 1. the negative values are the -ve
% values of Var with -1 corresponding to the maximum magnitude and 0
% corresponding to the minimum magnitude. The positive values are
% equivalent with 1 corresponding to the maximum magnitude
% For now all exponents assumed negative.
%
% INPUTS:
%
% Var = field (any size)
%
% minMag = lowest magnitude (power of 10) being mapped to 0.
%
% maxMag = highest magnitude (power of 10) mapped to +-1
%
%
%---------------------------------------------------------------------
%
% Dependencies; none. 
%
% Ryans ROMS Matlab and netcdf Utilities 17/7/13
%
%---------------------------------------------------------------------
%---------------------------------------------------------------------

%Take signed log:
VarOUTt = sign(-Var).*log10(abs(Var));

%Max/min out out of range values:
VarOUTt((VarOUTt<0)&(VarOUTt>maxMag)) = maxMag;
VarOUTt(VarOUTt<minMag) = minMag;
VarOUTt((VarOUTt>0)&(VarOUTt<(-maxMag))) = -maxMag;
VarOUTt(VarOUTt>(-minMag)) = -minMag;

%Map negative values:
VarOUT = VarOUTt;
VarOUT(VarOUTt<0) = (-VarOUT(VarOUTt<0)+minMag)/(maxMag-minMag);
%Map positive values:
VarOUT(VarOUTt>0) = (-VarOUT(VarOUTt>0)-minMag)/(maxMag-minMag);
end
