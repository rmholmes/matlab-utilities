function pausa(time);
% 
% pausa.m
% purpose:  wrapper for pause with messages
% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     11-Sep-2009
% modified: 11-Sep-2009
%
% obs: 
%

% force plotting of any figure
drawnow
if nargin == 1
         disp('Pausing...')
         disp(['for ',num2str(time),' seconds...'])
         pause(time)
         disp('...continuing')
else
         disp('...press key to continue...')
         pause
         disp('...continuing')
end