function [f95] = ff_montecarlo(M,edof,nsamples);
% ff_montecarlo.m - compute the confidence limits for significant eigenvalues in an EOF
%
% use:  [f95] = ff_montecarlo(M,dof,nsamples);;
% input:
%            M        -> spatial dimension
%            edof     -> temporal dimension N or the edof
%            nsamples -> realizations for Monte Carlo, the higher the better
%
% output:
%            f95 -> percent variance in each mode occuring by chance
%
% other m-files required: none
%
% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     18-May-2009
% modified: 18-May-2009
%
% obs: based on eofsig_montecarlo.m from John Wilkin
%      http://marine.rutgers.edu/dmcs/ms615/wilkin/matlab/
%
%      Create a set of random data (dimension M x N)
%      and finds the percent variance (f95) in each mode
%      for which there would be only a 5% probability of occuring by chance.
%
%      Jolliffe, I., 2002, Principal Component Analysis, 2nd ed., 502pp, Springer. chapter 6

if nargin < 3
  nsamples = 1000;
end

for i=1:nsamples
  if mod(i,round(nsamples/10)) == 0
    disp(['realization ' int2str(i)])
  end
  % create random matrix of Normal(0,1) distributed values
  Dmc = randn(M,edof);
  [Ur,Sr,Vr] = svd(Dmc,0);
  % fraction of variance
  fraction = diag(Sr*Sr')/trace(Sr*Sr');
  % preallocate for speed
  if i==1
    F = ones([length(fraction) nsamples]);
  end
  % save realization
  F(:,i) = fraction;
end

% find the 95-percentile value from the Monte Carlo set.
F   = sort(F,2);
f95 = F(:,round(95/100*nsamples));