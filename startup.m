%---------------------------------------------------------------------
%---------------------------------------------------------------------
%
% Startup script for plot settings in MATLAB 
%
%---------------------------------------------------------------------
%
% Dependencies; none. 
%
% Ryans ROMS Matlab and netcdf Utilities 22/2/16
%
%---------------------------------------------------------------------
%---------------------------------------------------------------------

%Add to path without using sudo:
% $$$ startdir = pwd;
% $$$ tmp = char(userpath);
% $$$ pathdir = tmp(1:(end-1));
% $$$ cd(pathdir)
% $$$ path(pathdef)
% $$$ cd(startdir)

set(0,'DefaultAxesColorOrder',[0 0 0], ...
      'DefaultAxesLineStyleOrder','-|--|:|-.')
set(0,'defaulttextfontsize',20);
set(0,'defaultaxesfontsize',20);
set(0,'DefaultFigureColor',[1 1 1]);
set(0,'defaultFigureRenderer','painters');
set(0,'DefaultTextInterpreter', 'latex')
set(0,'defaultlinelinewidth',1)
set(0,'defaultaxeslinewidth',1)
set(0,'defaultpatchlinewidth',1)
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');
set(0,'defaultColorbarTickLabelInterpreter','latex');

if (exist('/home/z3500785/Research/Data_Analysis'))
    cd('/home/z3500785/Research/Data_Analysis');
end
% $$$ cd 'C:\Users\RYANHO~1\Research\Data_Analysis\';
% $$$ addpath('C:\Users\RYANHO~1\.emacs.d\matlab-emacs\toolbox','-begin');
% $$$ rehash;
% $$$ emacsinit('C:\PROGRA~1\Emacs\bin\emacsclient.exe -n');

