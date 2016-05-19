function [sst,xlons,ylats] = read_goes(fn,varargin)
%
% FILENAME: read_goes.m
%
% USAGE: To run this program, use the following command:
%
%     [sst,xlons,ylats] = read_goes(FILENAME)
% or  [sst,xlons,ylats] = read_goes(FILENAME,CLOUD CONTAMINATION PERCENTAGE)')
%
% DESCRIPTION: This file contains one (1) Matlab program to read the GOES-11
%                 and 12 combined data set.  READ_GOES reads in sea surface
%                 temperature (sst) scaled values from GOES geostationary
%                 satellites equipped with the GOES IMAGER.  It's output is a
%                 sea surface temperature matrix with longitudes and latitudes
%                 (see below).  Optionally it produces figures of cloud density
%                 (Bayesian only), scaled sst values and sea surface
%                 temperature(SST).
%
%              The newer Bayesian data sets contain a lower case 'b' in the file
%                 name.  These contain a cloud possibility percentage for every
%                 pixel.  The default is to remove any SST value that has a
%                 cloud contamination possibility(pcloud) above 2%.  This
%                 threshold can be overwritten by the user to any value by the
%                 user.  The scales are quasi-logarithmic to offer more
%                 differentiation for various cloud probabilities.
% 
%
% NOTES:
%
%      1. This read software was created using Matlab version 7.4.0.287.
%
%      2. This program reads the GOES binary data as big-endian.
%
%      3. Please email all comments and questions concerning these routines
%           to podaac@podaac.jpl.nasa.gov.
%
% Created:
%         10 January 2008: G. Foti
%
% INPUT
% FN       (string) containing GOES data file name including path if necessary.
%             Use real name as it contains needed data
% VARARGIN (Optional). (number)  Can contain cloud contamination control
%             percentage.
%             If pre-Bayesian (no lower case 'b' in name) this variable will be
%                ignored.
%
%             The newer Bayesian cloud mask has a second set of values in the 
%                data files which scale to a probable percentage of cloud
%                contamination.  These scaled integer values (0-255) represent a
%                quasi-logarithmic scale to offer more differentiation at low
%                cloud probabilities. 
%                The cloud probability percentage maps to the scale using the
%                following equation:
%
%                'Scaled Value' =
%
%                2 - 40.546121 * ln( 0.002 + 'cloud probability percentage'/100)
%
%                                where "ln" is the natural logarithm.
%
%                The default contamination percentage is 2.0.  This is the same
%                percentage used to calculate the cloud contamination flag in
%                the previous data product (pre-Bayesian).  If the user wishes
%                to override the cloud probability percentage threshold enter the
%                percentage in this input variable.  For example if you wanted
%                to filter out SSTs with higher than 5.5% cloud contamination
%                probablity you would run:
%
%                    [sst,xlons,ylats] = read_goes(FILENAME,5.5);
%                
%                Running
%                   [sst,xlons,ylats] = read_goes(FILENAME,2.0) or
%                   [sst,xlons,ylats] = read_goes(FILENAME);
%                would produce the same results.
%                
% OUTPUT
% SST   An M x N array of Sea Surface Temperatures in degrees Celsius.  The
%          rectangle is parallel with lines on longitude and latitude. 
% XLONS A vector of length M representing longitude.
% YLATS A vector of length N representing latitude.
%
%======================================================================
% Copyright (c) 2008, California Institute of Technology
%======================================================================
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIPLAY ACKNOWLEDGEMENT INFO
%
disp(' ')
disp('If this data product is used for publication please include the')
disp('following acknowledgement:')
disp(' ')
disp('''The GOES Sea Surface Temperature data was obtained from the JPL')
disp('Physical Oceanography Distributed Active Archive Center (PODAAC).''')
disp(' ')
disp('and if appropriate:')
disp(' ')
disp('''Sea surface temperature image courtesy of Physical Oceanography')
disp('Distributed Active Archive Center (PODAAC).''')
disp(' ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK ARGUMENT IN COUNT
%
if nargin < 1 | nargin > 2
   disp('USAGE: [sst,xlons,ylats] = read_goes(filename)')
   disp('   or  [sst,xlons,ylats] = read_goes(filename,quality integer)')
   disp('function ended.  Debug or enter ''dbquit'' to exit program')
   disp('or enter ''return'' to continue')
   %   keyboard
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display(['Reading file ' fn]); 
dsname = 'SST'; %data set
instr  = 'GOES IMAGER'; %instrument name
sst_flag_val = 7; %Ignore all temperatures under this value as they are flags
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARSE INFO FROM FILENAME
%
i = findstr(fn,'sst'); %indicator of where the filename w/o directory begins
fn2 = fn(i(end):end); %strip off directory part of name
i = findstr(fn2, '_') ; %find where 2 or 3 dashes in name are

if length(i)==3 | length(i)==2
   yyyy=str2num(fn2(i(1)+1:i(1)+4)); %year
   yearday=str2num(fn2(i(2)+1:i(2)+3)); %day of year
   ver=fn2(i(1)-1); %if file is new Bayesian the value will be 'b'
   if length(i)==3 
      freq=fn2(4); %how many hours are averaged together to create this dataset.
      hr = [' hour ' fn2(i(3)+1:i(3)+2)]; %hour of averaged data
   else 
      freq= '24';
      hr = '';
   end
else
   disp('file name should have 2 or 3 dashes')
   disp('function ended.  Debug or enter ''dbquit'' to exit program')
   disp('or enter ''return'' to continue')
   %   keyboard
end

% text describing this data set
datestring=datestr(datenum(yyyy-1,12,31) + yearday); %string of data date
title_txt = [dsname ' from ' instr ': ' freq ' hour average on ' datestring hr];
disp(' ')
disp(['Processing ' title_txt])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET QUALITY FLAG VALUES WITH INPUT OR DEFAULTS
%
if strcmp(ver,'b')
   if nargin == 1
      %Ignore all SSTs with cloud contamination possibilities less than this
      %  value.  The default is 2.0%
      cloud_flag_pct = 2.0;
   else
      cloud_flag_pct = cell2mat(varargin); %user overwrites value
      if cloud_flag_pct > 100 | cloud_flag_pct < 0 
         disp('2nd optional input parameter if included must be a percentage')
         disp('value between 0 and 100.')
         disp('function ended.  Debug or enter ''dbquit'' to exit program')
         disp('or enter ''return'' to continue')
         %         keyboard
      end
   end
   %convert cloud probability percentage to scale of flag.
   cloud_flag_val = 2 - 40.546121 * log( 0.002 + cloud_flag_pct/100);
   disp(' ')
   disp(['Filtering SSTs with greater than ' num2str(cloud_flag_pct)...
         '% chance of cloud contamination']);
else
   if nargin == 2
      disp('file name does not contain lower case b and is therefore not')
      disp(' Bayesian.  The second input value will be ignored.')
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS THE DATA
%
fid=fopen(fn,'r','b'); % open file
if fid < 0
   disp(['Could not open file ' fn '.'])
   disp('Make sure path is correct and program exists.')
   disp('function ended.  Debug or enter ''dbquit'' to exit program')
   disp('or enter ''return'' to continue')
end
   
fseek(fid,0,-1); %make sure pointer is at beginning

%read in metadata    
if strcmp(ver,'b') %it is Bayesian data set with some meta data.
   header    = fgetl(fid); %get header info
   head_dat  = strread(header,'%s','delimiter',' '); %separate header meta data
   xlon_dim  = str2num(char(head_dat(2))); %length of longitude dimension
   ylat_dim  = str2num(char(head_dat(4))); %length of latitude dimension
   xlon1     = str2num(char(head_dat(5))); %first longitude indicator
   xlon_step = str2num(char(head_dat(6))); %width of longitudinal pixel
   ylat1     = str2num(char(head_dat(7))); %first longitude indicator
   ylat_step = str2num(char(head_dat(8))); %width of latitudinal pixel
   fseek(fid,xlon_dim,-1);%go to end of header line
else
   %pre 2008 data sets did not contain metadata, these are constant 
   xlon_dim  = 3000;
   ylat_dim  = 2100;
   xlon1     = 180.025;
   xlon_step = 0.0500;
   ylat1     = 60.025;
   ylat_step = 0.0500;
end

%calculate the lons & lats
xlons = -xlon1 + xlon_step*(1:xlon_dim);
ylats =  ylat1 - ylat_step*(1:ylat_dim);

%read the sst weighted values and cloud contamination values if data set is
% Bayesian
if strcmp(ver,'b') %get Bayesian cloud data

   % Get 2 sets of data.  The first ylat_dim rows are sst scaled values.  The
   %    next ylat_dim rows contain cloud contamination scaled values. 
   ds_data=fread(fid,[xlon_dim, ylat_dim * 2]);

   %seperate the data from the percentage of cloud contamination
   ds_pcloud = ds_data(:, ylat_dim+1:end); %cloud contamination matrix
   ds_data(:,ylat_dim+1:end)=[]; %leave sst scaled values in ds_data
else
   %pre Bayesian data sets were simply flat files containing sst scaled values
   ds_data=fread(fid,[xlon_dim,ylat_dim]); 
end

fid=fclose(fid); %close file

sst = ds_data; %copy scaled values to variable which will eventually hold sst
               %values in degrees C

%filter data for quality
if strcmp(ver,'b') % Bayesion
   %use pcloud values to filter data
   cloud_index=find(ds_pcloud < cloud_flag_val); %index of "too high" cloud
                                                  %probabilities
   sst(cloud_index)=NaN; %set cloud contaminated values to Not A Number.
else
   flag_index=find(ds_data < sst_flag_val); %index of sst scaled values not used as
                                         % they are flags
   sst(flag_index)=NaN; %set flagged values to Not A Number.
end

% convert scaled data to degrees centigrade
sst = sst * 0.15 - 3.15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT DATA SETS IF DESIRED by setting conditional to true to plot.
%
%set values to false for any plot not wanted
% plot cloud flag values
if false & strcmp(ver,'b') %change to false to skip plot
   figure(1)
   pcolor(xlons,ylats,ds_pcloud') 
   shading flat
   title(['pcloud:' title_txt])
   colorbar
end

% plot scaled sst values
if false %change to false to skip plot
   figure(2)
   pcolor(xlons,ylats,ds_data')
   shading flat
   title(['no filter no conv' title_txt])
   colorbar
end

% plot filtered sst values converted to centigrade
if false %change to false to skip plot
   figure(3)
   pcolor(xlons,ylats,sst')
   shading flat
   title([title_txt])
   colorbar
end
% keyboard
