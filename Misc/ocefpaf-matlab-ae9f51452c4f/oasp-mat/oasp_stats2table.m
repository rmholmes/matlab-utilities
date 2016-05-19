function [stats] = oasp_stats2table(data,tabfile)
%
% oasp_stats2table.m - save a table of basic statistics
%
% use:  [stats] = oasp_stats2table(data);
% input:
%    input - data, for different datasets use a cell array
%            tabfilele
%
% output:
%    output - mmean   -> mean
%             mmedian -> media
%             mmode   -> mode
%             vvar    -> variance
%             sstd    -> standard deviation
%             bvar    -> biased variance
%             bstd    -> biased standard deviation
%
% example:
%    [stats] = oasp_stats2table([{dataset1} {dataset2}]);
%
%    other m-files required: saveascii.m

% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     07-Oct-2009
% modified: 07-Oct-2009
%
% obs: use all "ignore nans" version
%

if ~iscell(data)
  f = find(isnan(data));
  if isempty(f);
    mmean   = mean(data);
    mmedian = median(data);
    mmode   = mode(data);
    vvar    = var(data);
    sstd    = std(data);
%      bvar    = sum((data-mmean).^2) / length(data);
%      bstd    = sqrt(bvar);
  else
    mmean   = nanmean(data);
    mmedian = nanmedian(data);
    mmode   = mode(data);
    vvar    = nanvar(data);
    sstd    = nanstd(data);
    data(f) = [];
%      bvar    = sum((data-mmean).^2) / length(data);
%      bstd    = sqrt(bvar);
  end
stats=[mmean,mmedian,mmode,vvar,sstd,bvar,bstd];
else
  for k=1:length(data)
    f = find(isnan(data{k}));
    if isempty(f);
      mmean   = mean(data{k});
      mmedian = median(data{k});
      mmode   = mode(data{k});
      vvar    = var(data{k});
      sstd    = std(data{k});
%        bvar    = sum((data{k}-mmean).^2) / length(data{k});
%        bstd    = sqrt(bvar);
    else
      mmean   = nanmean(data{k});
      mmedian = nanmedian(data{k});
      mmode   = mode(data{k});
      vvar    = nanvar(data{k});
      sstd    = nanstd(data{k});
      data(f) = [];
%        bvar    = sum((data{k}-mmean).^2) / length(data{k});
%        bstd    = sqrt(bvar);
    end
%    stats(k,:)=[mmean,mmedian,mmode,vvar,sstd,bvar,bstd];
  stats(k,:)=[mmean,mmedian,mmode,vvar,sstd];
  end
end

if nargin==2
  % make latex table using saveascii
%    saveascii('Mean Median Mode Variance Standard Deviation Biased-var Biased-std',tabfile)
%    saveascii(stats,tabfile,{'%g', '%g', '%g', '%g' , '%g', '%g', '%g'},{' & ', ' & ', ' & ', ' & ', ' & ', ' & '},'a')
  saveascii('Mean Median Mode Variance Standard Deviation',tabfile)
  saveascii(stats,tabfile,{'%g', '%g', '%g', '%g' , '%g'},{' & ', ' & ', ' & ', ' & '},'a')
  eval(['!sed -i ''s/$/\\\\/'' ', tabfile,''])
end

return