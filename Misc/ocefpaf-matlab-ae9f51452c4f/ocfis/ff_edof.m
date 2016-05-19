function [dof,Tscale,xc,lags] = ff_edof(D);
% ff_edof.m - Estimate degrees of freedom and auto-correlation scale
%             for the ensemble of time series in rows of data matrix D.
%
% use:  [dof,Tscale,xc,lags] = ff_edof(D);
%
% input:
%            rows of D are a set of time series of
%            length N = size(D,2) at M = size(D,1) locations.
%
% output:
%             dof    -> degrees of freedom
%             Tscale -> non-dimensional temporal scale
%             xc:    -> normalized auto-correlation
%             lags:  -> non-dimensional lags
%
% other m-files required: matlab's xcov.m, nanmean.m
%
% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     18-May-2009
% modified: 18-May-2009
%
% obs: based on degreesfreedom.m from John Wilkin
%      http://marine.rutgers.edu/dmcs/ms615/wilkin/matlab/
%      time series must be recorded at uniform sampling intervals

% compute the covariance of each time series
for m=1:size(D,1)
   % 'coeff' for normalized autocovariance (xc(0)=1)
   [xc(m,:),lags] = xcov(D(m,:),'coeff');
end

% average over all M time series to form
% a ensemble autocorelation. nanmean in case of masked data!
xc = nanmean(xc,1);

% retain half of auto-correlation, it is symmetric
xc = xc(find(lags==0):end); lags = lags(find(lags==0):end);

% integrate to first zero crossing, at t=t*
% find auto-correlation scale: Tscale = integral_0^t* xc(t') dt'
% Tscale will be in the normalized units of the sampling interval.
tstar = min(find(xc<=0))-1; Tscale = mean(xc(1:tstar))*tstar;

% effective degrees of freedom is estimated as the number
% of correlation time scales in length of the data record
dof = round(size(D,2)/Tscale);

%  compare with the example below
% function [tau,df] = time_scale(xx)
% computes the autocovariance of xx
% then computes the no. of degrees of freedom in xx
%
%  function [tau,df] = time_scale(xx)
% compute autocovariance
%  mlag=round(length(xx)/2);
%  rho=xcov(xx,mlag,'unbiased');
% keep only one side of the symmetric autocovariance
%  corr=rho(mlag+1:2*mlag+1)/max(rho);
%  lags=0:mlag;

% integral time scale (approximate)
% integrate only while positive
% USE TRAPEZOID RULE
%  tau=0.; i=1;
%  while corr(i)>0
%    tau=tau+(corr(i)+corr(i+1))/2;
%    i=i+1;
%  end

% find degrees of freedom
%  df=floor((length(xx)/tau));