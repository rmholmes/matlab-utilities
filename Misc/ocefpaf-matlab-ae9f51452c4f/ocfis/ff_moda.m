function function [moda] = ff_moda(data)
% ff_moda.m -> compute the mode and plot a histogram of modes...
%
%   call:  [moda] = ff_moda(data);
%
%  input:  data -> data to compute the mode
%
% output:  moda -> the mode computed
%
% example: [moda] = ff_moda(data);
%

% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     22-Feb-2009
% modified: 21-Nov-2009
%
% obs: 
%

[value,i,j] = unique(data);
h = hist(j,length(value));
m=max(h); moda=value(h==m);

hist(moda)
    ylabel('frequency','fontsize',20)
    xlabel('modes','fontsize',20)
    set(gca,'fontsize',20)