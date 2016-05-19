
function [IN] = Repl(IN,toreplace,replacement)
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%
% This function replaces elements of a matrix that are = a with b.
%
% INPUTS
%
% IN = in field (any size)
%
% toreplace = real number to replace
%
% replacement = number to replace toreplace with.
%
% OUTPUTS:
%
% New replaced IN.
%
%---------------------------------------------------------------------
%
% Dependencies; none. 
%
% Ryans ROMS Matlab and netcdf Utilities 17/7/13
%
%---------------------------------------------------------------------
%---------------------------------------------------------------------

if (isnan(toreplace))
    IN(isnan(IN)) = replacement;
else
    IN(IN == toreplace) = replacement;
end

end