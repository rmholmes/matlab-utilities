function digimg(input)
% digimg.m -> georeference an image
%
%   call:  digimg(input);
%
%  input:  input -> image
%
%
% example: digimg('map.png');
%

% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     15-Mar-2009
% modified: 28-Dec-2009
%
% obs: 
%

% read image file
fmt      = input(end-2:end);
[im,map] = imread(input,fmt);
imshow(im); % show image

% image limits menu
prompt   = {'south lat','north lat','east lon',' west lon'};
def      = {'-22.7','-19.7','-41','-38'};
dlgTitle = 'figure limites in degrees decimal';
lineNo   = 1; answer=inputdlg(prompt,dlgTitle,lineNo,def);

lat0 = str2num(answer{1}); latf = str2num(answer{2});
lon0 = str2num(answer{3}); lonf = str2num(answer{4});

% click to get positions
title('upper left than lower right','FontSize',20)
posic = ginput(2);
xmin  =  min(posic(:,1));
ymin  =  min(posic(:,2));
larg  = diff(posic(:,1));
altu  = diff(posic(:,2));

% crop image
im2 = imcrop(im,[xmin ymin larg altu]); close

% show crop referenced figure
iptsetpref('ImshowAxesVisible','on')
imshow(im2,map,'XData',[lon0 lonf],'YData',[latf lat0]); axis xy; hold on
iptsetpref('ImshowAxesVisible','off')