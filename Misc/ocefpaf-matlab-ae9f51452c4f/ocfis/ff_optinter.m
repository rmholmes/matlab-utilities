function [res,err, Theta] = ff_optinter(R,lambda,X,Y,data,xdata,ydata)
% ff_optinter.m -> % code based on OI exercise at:
%                  http://marine.rutgers.edu/dmcs/ms615/wilkin/matlab/oi_mercator.m
%                  The code uses optimal interpolation to map irregular
%                  spaced observations onto a regular grid.
%
%   call:  [res,err, Theta] = ff_optinter(R,lambda,X,Y,data,xdata,ydata);
%
%  input:  R            -> square root of the de-correlation length scale in units of sqrt(deg)
%          lambda       -> error squared to signal squared or E
%          X,Y          -> grid of the locations for Theta
%          data         -> observations
%          xdata, ydata -> observations location
%
% output:  Theta -> optimal interpolated data
%          err   -> estimated optimun error (add coment on perc of true)
%          res   -> residue fit
%
% example: [res,err, Theta] = ff_optinter(R,lambda,X,Y,data,xdata,ydata)
%
% subfunctions: ff_function.m

% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     21-Nov-2009
% modified: 21-Nov-2009
%
% obs: The covariance matrices are computed using a Gaussian shaped function
%      need to add ff_function to add more options
%

data = data(:);
disp(['length scale ' num2str(R) ' lon/lat'])
disp(['and error to signal variance ratio lambda of ' num2str(lambda)])
disp('Computing, please wait....')

%% Ars
xr = repmat(xdata(:) ,[1 length(xdata)]);
yr = repmat(ydata(:)',[length(xdata) 1]);

% compute the distance of each data point:
Rdist = sqrt((xr-xr').^2+(yr-yr').^2);

% covariance function.
Cdd0 = exp(-Rdist.^2/R.^2);

%% noise
noise = 1-lambda; % filipe

% Add error/signal variance to the diagonal:
% Final Data covariance Matrix between data points.
Cdd = noise*Cdd0 + lambda*eye(size(Cdd0)); % added noise (filipe)

%% Cxr
xr = repmat(xdata(:),[1 length(X(:))]);
yr = repmat(ydata(:),[1 length(X(:))]);
xg = repmat(X(:),  [1 length(xdata)]);
yg = repmat(Y(:),  [1 length(xdata)]);

% distance between observation r and grid g
Rmd = sqrt( (xr-xg').^2 + (yr-yg').^2 );

% again plug into covariance function
Cmd = noise*exp(-Rmd.^2/R.^2); Cmd=Cmd'; % added noise (filipe)

res = mean(data) + Cdd0*inv(Cdd)*(data-mean(data));
res = data(:)-res(:);

err = diag(1 - ( Cmd*inv(Cdd)*Cmd' )./noise ); % normalized by the error variance (filipe)
err = reshape(err,size(X))*100; % error in percentages

%% The one line estimation (the data were demeaned!)
% mean removed (filipe)
Theta = Cmd*(Cdd\data);
Theta = reshape(Theta,size(X));