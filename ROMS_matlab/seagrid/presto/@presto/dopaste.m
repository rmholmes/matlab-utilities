function theResult = dopaste(self, varargin)

% presto/dopaste -- Stub for "paste" event.
%  dopaste(self) returns self.

% svn $Id: dopaste.m 361 2009-07-02 15:43:20Z arango $
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 02-Nov-1999 23:18:31.
% Updated    02-Nov-1999 23:18:31.

if nargout > 0, theResult = []; end
if nargin < 1, help(mfilename), return, end

if nargout > 0
	theResult = self;
end