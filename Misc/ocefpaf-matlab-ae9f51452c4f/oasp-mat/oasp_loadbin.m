function [data,time]=oasp_loadbin(filein,bo)
%
% oasp_loadbin.m - load  OASP binary file
%
% use:  [data,time]=oasp_loadbin(filein,bo)
% input:
%    input  - filein: binary file
%           - bo: byte order -> 'be' for Big-endian 
%                               'le' for Little-endian%
% output:
%    output - data: binary data and date in matlab format
%           - time: matlab serial date format
%
% example:
%    [data,time]=oasp_loadbin('unix_ec011.tn','be')
%
% other m-files required: oasp_jhour2matlab.m

% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     25-Jun-2009
% modified: 25-Jun-2009
%
% obs: this version just loads data into matlab and convert
%      OASP date format to matlab serial date
%

%% check endianness, here is where the trick is!
  if bo == 'be'
    byteorder = 'ieee-be';
  else
    byteorder = 'ieee-le';
  end

%% open input file for reading
  [fidr,msg] = fopen(filein,'rb', byteorder);

%% check if file was successfully opened
  if fidr < 0, error(['cannot read file: ' filein]); end

%% get header
  % ascii part
  entries = 38; % is it 38 always ?! need to check !
  header = fread(fidr, entries, 'char'); headertxt = char(header');
  % numeric part
  [rec, count]  = fread(fidr, 4, 'double', byteorder); % start date
  [num, count]  = fread(fidr, 3, 'int32',  byteorder); % record length
  [inc, count]  = fread(fidr, 2, 'double', byteorder); % increment
%% data -> contain some kind of line break beteween them need to check later (FORTRAN?)
  [data, count] = fread(fidr, [3 num(3)], 'real*4' , byteorder);
                                            % real*4, and all other precision info, are hard coded.
                                            % if binary change this might be a source of error.

  data = data(3,:); % strip the the FORTRAN record index
  time = oasp_jhour2matlab(rec(2)):inc(2)/24:oasp_jhour2matlab(rec(4)); % time vector
