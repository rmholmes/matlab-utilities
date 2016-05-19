function oasp_bin2asc(filein,fileout,bo)
% oasp_bin2asc.m - convert OASP binary file to ASCII
%
% use:  oasp_bin2asc(filein,fileout,bo)
%
% input:
%           - filein: binary file
%           - fileout: ascii file
%           - bo: byte order -> 'be' for Big-endian 
%                               'le' for Little-endian
%
% example:
%    oasp_bin2asc('unix_ec011.tn','unix_ec011.asc','be')
%

% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     25-Jun-2009
% modified: 25-Jun-2009
%
% obs: need to test with a larger sample
%      some records spacing are unexplained
%

% check endianness, here is where the trick is!
if bo == 'be'
  byteorder = 'ieee-be';
else
  byteorder = 'ieee-le';
end

% open input file for reading
[fidr,msg] = fopen(filein,'rb', byteorder);

% check if file was successfully opened
  if fidr < 0, error(['cannot read file: ' filein]); end

%% get header
  % ascii part
  entries = 38; % there are some spaces and line breaks that need to be checked !
                % is it 38 always ?! need to check !
  header = fread(fidr, entries, 'char'); headertxt = char(header');

  % numeric part
  [rec, count]   = fread(fidr, 4, 'double', byteorder); % start date
  [num, count]   = fread(fidr, 3, 'int32', byteorder);  % record length
  [inc, count]   = fread(fidr, 2, 'double', byteorder); % increment
                                                        % they all contain some kind of line break beteween them
                                                        % need to check later

%% data
  [data, count]   = fread(fidr, [3 num(3)], 'real*4' , byteorder);
                                            % real*4, and all other precision info, are hard coded.
                                            % if binary change this might be a source of error.

  % write ascii output
  [fidw,msg] = fopen(fileout,'w');
  % check if file was created
  if fidw < 0, error(['cannot write file: ' fileout]); end

  % write header
  fprintf(fidw, '%s\n',['% filename ' headertxt(5:13)]);
  fprintf(fidw, '%s\n',['% date ' headertxt(23:30)]);
  fprintf(fidw, '%s\n',['% time ' headertxt(31:38)]);
  fprintf(fidw, '%s %6.8f\n','% start',rec(2));
  fprintf(fidw, '%s %6.8f\n','% finish',rec(4));
  fprintf(fidw, '%s\n',['% num ' num2str(num(3))]);
  fprintf(fidw, '%s %E\n','% inc',inc(2));

  % write data
  fprintf(fidw, '%g \n',data(3,:)');
% here the 3 by length (num(3)) is a "magic" number that needs an explanation !!!
% probably only by reading the full original fortran code and/or asking the guy that
% wrote the binary file... I'm throwing away the first two colouns and keeping onlyt the data

%% close ASCII file
  status = fclose(fidw);
  % check if file was successfully closed
  if status < 0, error(['cannot write file: ' fileout]); end

return
