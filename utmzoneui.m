function [zone, varargout] = utmzoneui(arg)
%UTMZONEUI Choose or identify UTM zone by clicking on map
%
%  ZONE = UTMZONEUI will create an interface for choosing a UTM zone on a world
%  display map.  It allows for clicking on an area for its appropriate zone, or 
%  entering a valid zone to identify the zone on the map.
%
%  ZONE = UTMZONEUI(InitZone) will initialize the displayed zone to the zone 
%  string given in INITZONE.

% Copyright 1996-2013 The MathWorks, Inc.
% Written by A. Kim

msg = '';

if nargout > 1
    warnObsoleteMSGSyntax(mfilename)
    varargout{1} = msg;
end

if nargin == 0
    initzone = [];
    action = 'initialize';
elseif nargin == 1
    if length(arg)>3	% action string
        initzone = [];
        action = arg;
    else	% initial zone
        [~,~] = utmzone(upper(arg));
        initzone = upper(arg);
        action = 'initialize';
    end
end

zone = '';

if isempty(action),	action = 'initialize'; end

if strcmp(action,'initialize')

	hndl = initgui(initzone);
	set(gcf,'UserData',hndl)

	set(gca,'ButtonDownFcn','utmzoneui(''zoneclick'');');	%  Disable default button down function
	hch = get(gca,'Children');
	set(hch,'ButtonDownFcn','utmzoneui(''zoneclick'');');	%  Disable default button down function

	while 1

		uiwait(hndl.fig);
        if ~ishghandle(hndl.fig)
            return
        end

		btn = get(hndl.fig,'CurrentObject');

		zone = upper(get(hndl.zoneedit,'String'));
        try
		    [latlim,lonlim] = utmzone(zone);
        catch e
            latlim = [];
            lonlim = [];
            msg = e.message;
        end

		if isempty(btn)
			disp('???')
		elseif btn==hndl.zoneedit

			hndl = get(gcf,'UserData');
			hboxold = get(hndl.zoneedit,'UserData');

			if ~isempty(msg)   %  Error condition
				uiwait(errordlg(msg,'Invalid Zone','modal'))
			else         %  Valid zone
				if ~isempty(hboxold),  delete(hboxold),  end
				ltbox = [latlim(1) latlim(2) latlim(2) latlim(1) latlim(1)];
				lnbox = [lonlim(1) lonlim(1) lonlim(2) lonlim(2) lonlim(1)];
				hbox = patchm(ltbox,lnbox,-1,[1 .2 .2]);

				set(hndl.zoneedit,'String',zone)
				set(hndl.zoneedit,'UserData',hbox)
			end

		elseif btn==hndl.apply

			hndl = get(gcf,'UserData');

			if isempty(msg)
				zone = upper(get(hndl.zoneedit,'String'));
			else
				zone = [];
			end

 			delete(hndl.fig)
 			break

		else	% cancel button pushed

			zone = initzone;
			delete(hndl.fig)
			break

		end
	
	end

elseif strcmp(action,'zoneclick')

	hndl = get(gcf,'UserData');
	hboxold = get(hndl.zoneedit,'UserData');

	mat = gcpmap; lt = mat(1,1); ln = mat(1,2);
	mat2 = get(gca,'CurrentPoint');

    try
	    zone = utmzone(lt,ln);
    catch e
        zone = [];
        msg = e.message;
    end
	
	if mat2(1,1)>pi || mat2(1,1)<-pi
		msg = 'Coordinates not within UTM zone limits';
	end

	if ~isempty(msg)   %  Error condition
		zone = [];
		uiwait(errordlg(msg,'Invalid Zone','modal'))
	else         %  Valid zone
		[latlim,lonlim] = utmzone(zone);
		if ~isempty(hboxold),  delete(hboxold),  end
		ltbox = [latlim(1) latlim(2) latlim(2) latlim(1) latlim(1)];
		lnbox = [lonlim(1) lonlim(1) lonlim(2) lonlim(2) lonlim(1)];
		hbox = patchm(ltbox,lnbox,-1,[1 .2 .2]);

		set(hndl.zoneedit,'String',zone)
		set(hndl.zoneedit,'UserData',hbox)
	end

	set(gcf,'UserData',hndl)
end

%**************************************************************************
%**************************************************************************
%**************************************************************************

function h = initgui(initzone)

PixelFactor = guifactm('pixels');
FontScaling =  guifactm('fonts');

h.fig = figure('NumberTitle','off', 'Name','Pick UTM Zone', ...
       'Units','Points', 'Position',PixelFactor*[96 128 894 616], ...
       'Resize','off', 'WindowStyle','modal', 'Visible','off');
%
%adjust window position if corners are offscreen
%
shiftwin(h.fig)

colordef(h.fig,'white')
figclr = get(h.fig,'Color');
frameclr = brighten(figclr,0.5);

%  Display map of world with utm zone designations
lts = [-80:8:72 84]';
lns = (-180:6:180)';

axesm('miller','maplatlim',[-80 84],'maplonlim',[-180 180],...
	  'mlinelocation',lns,'plinelocation',lts,...
	  'mlabellocation',-180:24:180,'plabellocation',lts)
framem; gridm; mlabel; plabel
set(gca,'position',[.02 .12 .96 .87],'xlim',[-3.5 3.3],'ylim',[-2 2.2])

load('coast');
hold on
geoshow(lat, long, 'Color','k');
clear lat long

if ~isempty(initzone)
	[latlim,lonlim] = utmzone(initzone);
	ltbox = [latlim(1) latlim(2) latlim(2) latlim(1) latlim(1)];
	lnbox = [lonlim(1) lonlim(1) lonlim(2) lonlim(2) lonlim(1)];
	hbox = patchm(ltbox,lnbox,-1,[1 .2 .2]);
else
	hbox = [];
end

%  Set up gui display below map

uicontrol(h.fig, 'Style','frame', ...
	             'Units','normalized', 'Position',[0.2 0.02 0.15 0.08], ...
	             'FontSize',FontScaling*12, 'FontWeight','bold', ...
	             'BackgroundColor',frameclr, 'ForegroundColor','black');

uicontrol(h.fig, 'Style','text', 'String','Zone:', ...
	             'Units','normalized', 'Position',[0.21 0.03 0.065 0.05], ...
	             'FontSize',FontScaling*12, 'FontWeight','bold', ...
	             'HorizontalAlignment','left',  ...
	             'BackgroundColor',figclr, 'ForegroundColor','black');

h.zoneedit = uicontrol(h.fig, 'Style','edit', 'String',initzone, ...
	             'Units','normalized', 'Position',[0.28 .035 0.06 0.05], ... 
	             'FontSize',FontScaling*12,    'FontWeight','bold', ...
	             'HorizontalAlignment','center', 'UserData',hbox, ...
	             'BackgroundColor',figclr, 'ForegroundColor','red', ...
				 'CallBack','uiresume');

h.apply = uicontrol(h.fig, 'Style','push', 'String','Accept', ...
	        'Units','normalized', 'Position',[0.46 0.02 0.1 0.08], ...
	        'ForegroundColor','black', 'BackgroundColor',figclr, ...
	        'FontName','Helvetica', 'FontSize',FontScaling*12, ...
			'FontWeight','bold', 'CallBack','uiresume');

h.help = uicontrol(h.fig, 'Style','push', 'String','Help', ...
	        'Units','normalized', 'Position',[0.58 0.02 0.1 0.08], ... 
	        'ForegroundColor','black', 'BackgroundColor',figclr, ...
	        'FontName','Helvetica', 'FontSize',FontScaling*12, ...
			'FontWeight','bold',  'Interruptible','on', ...
			'Callback', @(~,~) doc('utmzoneui'));

h.cancel = uicontrol(h.fig, 'Style','push', 'String','Cancel', ...
	        'Units','normalized', 'Position',[0.7 0.02 0.1 0.08], ...  
	        'ForegroundColor','black', 'BackgroundColor',figclr, ...
	        'FontName','Helvetica', 'FontSize',FontScaling*12, ...
			'FontWeight','bold', 'CallBack','uiresume');

set(h.fig,'Visible','on')
