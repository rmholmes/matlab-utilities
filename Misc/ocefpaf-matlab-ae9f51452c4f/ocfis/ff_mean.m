function [ci,des,med] = ff_mean(x,dof,conf);
% ff_mean.m - computes mean, std and confidence interval for large data sets
%
% use:  [ci,des,med] = ff_mean(x,dof,conf);
% input:
%    input - data
%            degree of freedom
%            confidence
%
% output:
%    output -  ci: confidence interval
%             std: standard deviation
%            mean: mean
%
% example: [ci,std,mean] = ff_mean(x,length(x),0.95);
%
%
% other m-files required: mean, std
%
% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     16-May-2009
% modified: 16-May-2009
%
% obs:
%

med = mean(x);
des = std(x);
se = des / sqrt(length(x)));
ci  = abs(icdf('normal',(1-conf)/2,0,1)) * se;

%N = 1000;             % sample size
%X = exprnd(3,N,1);    % sample from a non-normal distribution
%mu = mean(X);         % sample mean (normally distributed)
%sig = std(X)/sqrt(N); % sample standard deviation of the mean
%alphao2 = .05/2;      % alpha over 2
%CI = [mu + norminv(alphao2)*sig ,...
      %mu - norminv(alphao2)*sig  ]

%CI =

%2.9369    3.312
