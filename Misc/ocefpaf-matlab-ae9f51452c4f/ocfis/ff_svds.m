function function [ec,eof,var] = ff_svds(data,dim,modes)
% ff_svds.m -> Compute EOF via svds
%
%   call:  [ec, eof, var] = ff_svds(data,dim,modes);
%
%  input:  data  -> must be demeand and [n,p]=size(data) n->time, p->station
%          dim   -> dimension for reshaping
%          modes -> number of modes for decomposition
%
% output:  ec   -> expansion coefficients time series
%          eof  -> eigvectors
%          var  -> variance explained
%
% example: [ec, eof, var] = ff_svds(data,dim,modes);
%

% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     02-Nov-08
% modified: 21-Nov-2009
%
% obs: Use matlab internal
%      var = diag(S*S')/trace(S*S') % if full svd is used
%

%% EOF by svds (also try svd and eig)
[U,S,V]=svds(data,modes);

eign = diag(S).^2;             % eignvalues
ec   = U*S;                    % expansion coeff
var  = eign/trace(data'*data); % variance explained
eof  = reshape(V,[dim modes]); % reshape eof back to i,j or lon,lat