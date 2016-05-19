%% filipe is just a place holder and is chande to $USER via bashrc alias that call matlab%%

%% add my toolboxes path at the end of the path
addpath(genpath('/home/filipe/svn-tools/mymatlab'),'-end')
addpath(genpath('/home/filipe/svn-tools/ocefpaf-matlab'),'-end')
addpath(genpath('/home/filipe/svn-tools/mat-svn-toolboxes'),'-end')

cls

%% snctools with java
javaaddpath('/home/filipe/svn-tools/mymatlab/netcdfAll-4.2.jar','-end');
setpref ( 'SNCTOOLS', 'USE_JAVA', true );

%% Netcdf Toolbox global options (turn autoscale & autonan on)
javaaddpath('/home/filipe/svn-tools/mymatlab/njTools-2.0.10_jre1.6.jar','-end');
global nctbx_options
nctbx_options.theAutoNaN=1;
nctbx_options.theAutoscale=1;

%% figure defaults
set(0,'DefaultAxesColor',      'white',...          % axes color white
      'DefaultFigureColor',    'white', ...         % figure color white
      'DefaultAxesTickLength',  [0.01, 0.01], ...   % tick length
      'DefaultAxesFontSize',    16, ...             % axes font size
      'DefaultTextFontSize',    16, ...             % text font size
      'DefaultTextVerticalAlignment', 'middle', ... % vertical font alignment
      'DefaultTextHorizontalAlignment', 'left', ... % horizontal font alignment
      'DefaultFigurePaperPositionMode','auto')      % WYSIWYG
