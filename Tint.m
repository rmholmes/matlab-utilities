function integrated_field = Tint(tvec,field);
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%
% This function time integrates the time series/series' in field (time
% must be last dimension) and returns a vector of the time integrated
% values at each time. tvec is in days.
%
% INPUTS:
%
% tvec = vector of time values in days
%
% field = some field where time is the last dimension.
%
% OUTPUTS:
%
% integrated_field = time integrated field (same size as field). 
%
%---------------------------------------------------------------------
%
% Dependencies; none. 
%
% Ryans ROMS Matlab and netcdf Utilities 6/8/13
%
%---------------------------------------------------------------------
%---------------------------------------------------------------------

sizef = size(field);
tL = length(tvec);
integrated_field = zeros(sizef);
if (sizef(2) == 1)
for t = 2:tL
    integrated_field(t) = integrated_field(t-1)+86400*...
        (tvec(t)-tvec(t-1))*(field(t)+field(t-1))/2;
end
elseif (length(sizef) == 2)
for t = 2:tL
    integrated_field(:,t) = integrated_field(:,t-1)+86400*...
        (tvec(t)-tvec(t-1))*(field(:,t)+field(:,t-1))/2;
end
elseif (length(sizef) == 3)
for t = 2:tL
    integrated_field(:,:,t) = integrated_field(:,:,t-1)+86400*...
        (tvec(t)-tvec(t-1))*(field(:,:,t)+field(:,:,t-1))/2;
end
elseif (length(sizef) == 4)
for t = 2:tL
    integrated_field(:,:,:,t) = integrated_field(:,:,:,t-1)+86400*...
        (tvec(t)-tvec(t-1))*(field(:,:,:,t)+field(:,:,:,t-1))/2;
end    
end
end
