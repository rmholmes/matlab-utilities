function oasp_bin2dat(filein,fileout,bo)
%
% oasp_bin2dat.m - convert OASP binary file to ASCII
%
% use:  fileout = oasp_bin2dat(filein,fileout,bo)
% input:
%    input  - filein: binary file
%           - fileout: ascii file
%           - bo: byte order -> 'be' for Big-endian 
%                               'le' for Little-endian

% example:
%    oasp_bin2dat('unix_ec011.tn','unix_ec011.dat','be')
%
% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     25-Jun-2009
% modified: 25-Jun-2009
%
% obs: this version is different from oasp_bin2asc it also
%      saves the time vector in matlab serial date format
%

% check endianness, here is where the trick is!
if bo == 'be'
  byteorder = 'ieee-be';
else
  byteorder = 'ieee-le';
end

% open inoput file for reading
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
  [num, count]   = fread(fidr, 3, 'int32',  byteorder); % record length
  [inc, count]   = fread(fidr, 2, 'double', byteorder); % increment
                                                        % they all contain some kind of line break beteween them
                                                        % need to check later

%% data
  [data, count]   = fread(fidr, [3 num(3)], 'real*4' , byteorder);
                                            % real*4, and all other precision info, are hard coded.
                                            % if binary change this might be a source of error.

%% here the 3 by length (num(3)) is because FORTRAN add a record dimension every data
  data=data(3,:);
%% time vector
  time=oasp_jhour2matlab(rec(2)):inc(2)/24:oasp_jhour2matlab(rec(4));

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
  fprintf(fidw, '%g \t %e\n',[data;time]);

%% close ASCII file
  status = fclose(fidw);
  % check if file was successfully closed
  if status < 0, error(['cannot write file: ' fileout]); end

return