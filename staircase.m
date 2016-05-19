function [outx,outy] = staircase(xb,y)
% This function constructs a staircase data set. xb is the position of
% the grid cell transitions, and y is the value at the grid cells
% (thus length(xb) = length(y)+1). The output is a staircase data set
% with transitions at the half-way points between the grids.

Nt = length(xb);
outx = [];
outy = [];

for ii=1:length(y)
    outx = [outx xb(ii) xb(ii+1)];
    outy = [outy y(ii) y(ii)];
end
end