function varargout = matslise(varargin)
% MATSLISE M-file for matslise.fig


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @matslise_OpeningFcn, ...
                   'gui_OutputFcn',  @matslise_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before matslise is made visible.
function matslise_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to matslise (see VARARGIN)

% Choose default command line output for matslise
handles.output = hObject;

set(handles.xmin_edit,'String','0');
set(handles.xmax_edit,'String','inf');
set(handles.a0_edit,'String','1');
set(handles.a1_edit,'String','1');
set(handles.b0_edit,'String','0');
set(handles.b1_edit,'String','0');

startup;

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = matslise_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;





% --- Executes on button press in compute_button.
function compute_button_Callback(hObject, ed, handles)
try
    set(handles.info_edit,'ForeGroundColor','k');
    %disable components while computation:
    hs=findobj(gcf,'Style','pushbutton','Enable','on','-or',...
        'style','text','-or','style','edit','-or','style','listbox',...
        '-or','style','popupmenu','-or','style','listbox');
    set(hs,'Enable','off')
    set(handles.info_edit,'String','Busy...','Enable','on')
    drawnow
    %str=get(handles.info_edit,'String');
    try
        [p,q,w,xmin,xmax,a0,b0,a1,b1]=readInputFields(handles);
    catch
        %errordlg('The problem specification fields need to be filled out correctly.','Error','modal');
        %set(handles.info_edit,'String',str)
        set(handles.info_edit,'String','The problem specification fields need to be filled out correctly.');
        set(handles.info_edit,'ForeGroundColor','r');
        set(hs,'Enable','on')
        return;
    end
    %create slp-object:
    handles.slp=slp(p,q,w,xmin,xmax,a0,b0,a1,b1);
    %remove old output:
    set(handles.eigenvalues_listbox,'String','');
    set(handles.eigenfunctionsData_listbox,'String','');
    set([handles.eigenfunctionWorkspace_button],'Enable','off')
    drawnow
    %display classification information:
    set(handles.classifInfo_edit,'String',classificationInfo(handles.slp))
    %start computation of eigenvalues
    pmin=eval(get(handles.pmin_edit,'String'));
    pmax=eval(get(handles.pmax_edit,'String'));
    tol=eval(get(handles.tol_edit,'String'));
    indices=get(handles.indices_menu,'Value');
    if indices>1
        indices=false;
    else
        indices=true;
    end
    warning('off');
    [E,handles.meshData] = computeEigenvalues(handles.slp,pmin,pmax,tol,indices);
    warning('on');
    if ~E.success && length(E.msg)>1
        error(E.msg)
    end
    data = [E.indices' E.eigenvalues' E.errors']; 
    %display eigenvalue results:
    format=['%' num2str(max(ceil(log10(E.indices(end)))+1,4)) '.0f'];
    s1 = sprintf(format,data(1,1));
    tmp=max(abs(E.eigenvalues(end)),abs(E.eigenvalues(1)));
    format2=['%22.' num2str(17-max(ceil(log10(tmp)),1)) 'f      '];
    s2 = sprintf(format2,data(1,2));
    s3 = sprintf('%5.2e  ',data(1,3));
    resultsStr = {[s1,s2,s3]};
    for i=2:length(E.indices)
            s1 = sprintf(format,data(i,1));
            s2 = sprintf(format2,data(i,2));
            s3 = sprintf('%5.2e  ',data(i,3));
            resultsStr = [resultsStr;...
                {[s1,s2,s3]}] ;
    end
    set(handles.eigenvalues_listbox,'Value',1)
    set(handles.eigenvalues_listbox,'String',resultsStr);
    handles.eigenvalueData=data;
    set([handles.plot_button,handles.workspace_button, handles.eigenfunctionCompute_button],'Enable','on')
    set(hs,'Enable','on')
    resultsStr = {['Eigenfunction ',num2str(E.indices(1)) ]};
    for i=2:length(E.eigenvalues)
        if E.indices(i)>99
            resultsStr = [resultsStr;{['Eigenf. ',num2str(E.indices(i))]}] ;
        else
            resultsStr = [resultsStr;{['Eigenfunction ',num2str(E.indices(i))]}] ;
        end
    end
    set(handles.eigenfunctions_popupmenu, 'Value', 1);
    set(handles.eigenfunctions_popupmenu,'String',resultsStr);
    if any(E.status)
        tmp=find(E.status);
        s=sprintf('WARNING: DIFFICULTIES ENCOUNTERED IN THE COMPUTATION OF EIGENVALUES WITH INDEX ');
        for i=tmp
            s=[s num2str(E.indices(i)) ', '];
        end
        s=[s 'THE RESULT MAY BE UNACCURATE.'];
        set(handles.info_edit,'String',s)
        set(handles.info_edit,'ForeGroundColor','r');
    else
        if length(E.eigenvalues)>1
           set(handles.info_edit,'String',[num2str(length(E.eigenvalues)) ' eigenvalues computed.'])
        else
           set(handles.info_edit,'String','1 eigenvalue computed.') 
        end
    end
    %check close eigenvalues
    if any(abs(diff(E.eigenvalues))<max(tol,1e-10))
        s1= get(handles.info_edit,'String');
        s2=sprintf('\nClose eigenvalues detected. The computation of the eigenfunctions of the clustered eigenvalues is very ill-conditioned, eigenfunction results may be inaccurate.');
        s=strcat(s1,s2);
        set(handles.info_edit,'String',s);
        set(handles.info_edit,'ForeGroundColor','r');
    end
    guidata(hObject, handles);
catch ME
    set(hs,'Enable','on')
    %set(handles.info_edit,'String','')
    %errordlg(lasterr,'Error','modal');
    set(handles.info_edit,'ForeGroundColor','r');
    set(handles.info_edit,'String',ME.message)
end

function [p,q,w,xmin,xmax,a0,b0,a1,b1]=readInputFields(handles)
%reads the data from the input fields
[paramNames,paramValues] = processParameters(get(handles.parameterEdit,'String'));
q=get(handles.q_edit,'String');
p=get(handles.p_edit,'String');
w=get(handles.w_edit,'String');
xmin=get(handles.xmin_edit,'String');
xmax=get(handles.xmax_edit,'String');
a0=get(handles.a0_edit,'String');
a1=get(handles.a1_edit,'String');
b0=get(handles.b0_edit,'String');
b1=get(handles.b1_edit,'String');
for i=1:length(paramNames)
    q=strreplace(q,paramNames{i},paramValues{i});
    p=strreplace(p,paramNames{i},paramValues{i});
    w=strreplace(w,paramNames{i},paramValues{i});
    xmin=strreplace(xmin,paramNames{i},paramValues{i});
    xmax=strreplace(xmax,paramNames{i},paramValues{i});
    a0=strreplace(a0,paramNames{i},paramValues{i});
    a1=strreplace(a1,paramNames{i},paramValues{i});
    b0=strreplace(b0,paramNames{i},paramValues{i});
    b1=strreplace(b1,paramNames{i},paramValues{i});
end
xmin=eval(xmin);
xmax=eval(xmax);
a0=eval(a0);
a1=eval(a1);
b0=eval(b0);
b1=eval(b1);
p=str2func(['@(x)' p]);
q=str2func(['@(x)' q]);
w=str2func(['@(x)' w]);


function s=strreplace(s,old,new)
    %strrep such that e.g. strrep('n*sin(x)','n','5') returns '5*sin(x)'
    %and not '5*si5(x)'
    inds=strfind(s,old);
    while ~isempty(inds)
        ind=inds(1);
        if (ind>1 && isletter(s(ind-1))) || (ind(end)+length(old)<=length(s) && isletter(s(ind+length(old))))
            ind=ind+1;
        else
            s=[s(1:ind-1) '(' new  ')' s(ind+length(old):end)];
            ind=ind+length(new);
        end
        inds=strfind(s(ind:end),old)+ind-1;
    end
 
%-------------------------------------------------------------------
 function [names,values] = processParameters(s)
 try
 remainder = s;
 i=1;
 names={}; values={};
 while any(remainder)
    [chopped,remainder] = strtok(remainder,',');
    ind=find(chopped=='=');
    names{i}=strtrim(chopped(1:ind-1));
    values{i}=chopped(ind+1:end) ;
    i=i+1;
 end
 catch
    error('Error in the evaluation of the parameters.');
 end


function File_menu_Callback(h, ed, handles)



% --------------------------------------------------------------------
function Open_menu_Callback(hObject, ed, handles)
[filename, pathname] = uigetfile( ...
{'*.mat','MAT-files (*.mat)'}, ...
    'Select the problem to open');
if isequal([filename,pathname],[0,0])
    return;
end
file=fullfile(pathname,filename);
if exist(file) == 2
    S = load(file);
else
    errordlg('file doesn''t exist','Error','modal');
    return;
end
clearFields(handles)
set(handles.p_edit,'String',S.p);
set(handles.q_edit,'String',S.q);
set(handles.w_edit,'String',S.w);
set(handles.tol_edit,'String',S.tol);
set(handles.a0_edit,'String',S.a0);
set(handles.a1_edit,'String',S.a1);
set(handles.b0_edit,'String',S.b0);
set(handles.b1_edit,'String',S.b1);
set(handles.xmin_edit,'String',S.a);
set(handles.xmax_edit,'String',S.b);
set(handles.parameterEdit,'String',S.param);



function clearFields(handles)
hs=findobj(handles.SLP_panel,'style','edit');
set(hs,'String','')
set(handles.info_edit,'String','');
set(handles.classifInfo_edit,'String','');
set(handles.eigenvalues_listbox,'String','');
set(handles.eigenfunctionsData_listbox,'String','');
set([handles.plot_button,handles.workspace_button,handles.eigenfunctionWorkspace_button,...
    handles.eigenfunctionCompute_button],'Enable','off')



% --------------------------------------------------------------------
function Save_menu_Callback(hObject, ed, handles)
% hObject    handle to Save_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
[FileName,PathName,FilterIndex]=uiputfile({'*.mat','MAT-files (*.mat)'},'Save problem as ...');
if strcmp(FileName,'0') || strcmp(PathName,'0') || FilterIndex == 0
    return;
end
p = get(handles.p_edit,'String');
q = get(handles.q_edit,'String');
w = get(handles.w_edit,'String');
a = get(handles.xmin_edit,'String');
b = get(handles.xmax_edit,'String');
if isempty(w) || isempty(q) || isempty(p) || isempty(a) || isempty(b)
    error('You should fill in all the fields');
end
a0 = get(handles.a0_edit,'String');
b0 = get(handles.b0_edit,'String');
a1 = get(handles.a1_edit,'String');
b1 = get(handles.b1_edit,'String');
if isempty(a0) || isempty(a1) || isempty(b0) || isempty(b1) 
    error('You must enter a numeric value for a0, b0 ,a1 ,b1');
end
tol = get(handles.tol_edit,'String');
if isempty(tol) 
    error('You must enter a positive numeric value for tol');
end
param = get(handles.parameterEdit,'String');
file=fullfile(PathName,FileName);
save(file,'p','q','w','a','b','a0','b0','a1','b1','tol','param');
catch
    errordlg(lasterr,'Error','modal');
end


% --- Executes on button press in plot_button.
function plot_button_Callback(hObject, eventdata, handles)
    %plot of eigenvalue results
fhandle = figure('IntegerHandle','off','NumberTitle','off',...
    'Visible','off','Name','Eigenvalues plot');
movegui(fhandle,'west')
ahandle = axes('Parent',fhandle);
plot(handles.eigenvalueData(:,1),handles.eigenvalueData(:,2),'bo',...
     'MarkerFaceColor','b','MarkerSize',5,'Parent',ahandle);
set(ahandle,'XGrid','on');
set(ahandle,'YGrid','on');
if length(handles.eigenvalueData(:,1)) <50
   set(ahandle,'xtick',handles.eigenvalueData(:,1))
   set(ahandle, 'XTickMode', 'manual');
end
set(get(ahandle,'XLabel'),'String','eigenvalue index','FontSize',11);
legend(ahandle,'eigenvalues','Location','Best');
set(fhandle,'Visible','on');
set(fhandle,'HandleVisibility','off');
if handles.slp.classification.LNF %show potential
    fhandle = figure('IntegerHandle','off','NumberTitle','off',...
    'Visible','off','Name','Eigenvalues + Potential');
    movegui(fhandle,'center')
    ahandle = axes('Parent',fhandle);
    a = handles.meshData.Infa;
    b = handles.meshData.Infb;
    step = min(abs(b-a)/1500,3);
    w = a:step:b;
    VV = arrayfun(handles.slp.q,w);
    plot(w,VV,':k','Parent',ahandle);
    indices=handles.eigenvalueData(:,1);
    eigenvalues=handles.eigenvalueData(:,2);
    colormap = jet(length(indices));
    legendStr{1} = 'Potential function';
    hold on;
    %to make the colors in the legend match with the color of the eigenvalues_window:
    for j = length(indices):-1:1
       plot(b+2,eigenvalues(j),'Color',colormap(j,:),'Parent',ahandle);
       hold on;
       legendStr = {legendStr{:} ['E_{', num2str(indices(j)),'}']}; 
    end
    savedLines=[];
    for j = 1:length(indices)
      lines = []; %linesegments representing energy levels
      eigv = eigenvalues(j); %current eigenvalue
      start = a; %start of linesegment
      for i=1:length(VV)-1
        if VV(i) >= eigv && VV(i+1) < eigv
            start = (w(i) + w(i+1))/2; 
        end
        if VV(i) < eigv && VV(i+1) >= eigv
            lines = [lines; start (w(i) + w(i+1))/2];
        end
      end
      if VV(end) < eigv
          lines = [lines; start b];
      end
     %draw these linesegments
     for s=1:size(lines,1)
       plot([lines(s,1) lines(s,2)],[eigv eigv],'Color',colormap(j,:),'Parent',ahandle);
       hold on;
     end
     if isempty(lines)
         plot([a b],[eigv eigv],'Color',colormap(j,:),'Parent',ahandle);
     end
      if j==length(indices)
         savedLines=lines;
     end
    end
    %find good xlimits and ylimits for the axis: 
    if isempty(savedLines)
        xmax = b;
        xmin = a;
    else
      xmax = min(max([max(savedLines) a])+1,b); 
      xmin = max(min([min(savedLines) b]-1),a);
    end
    hold off;
    ymin = min(min(eigenvalues),min(VV(15:end)));
    ymax = max(eigenvalues);
    ymin=ymin-(ymax-ymin)/30;
    if length(indices) > 1
     [s1,s2]=size(eigenvalues);
     if s1<s2
       eigvdiff = [eigenvalues 0]-[0 eigenvalues];
     else
       eigvdiff = [eigenvalues; 0]-[0; eigenvalues];  
     end
     eigvdiff = min(eigvdiff(2:length(eigvdiff)-1));
     if (eigvdiff/(ymax-ymin)) < 0.005
      ymin = min(eigenvalues) - abs(ymax-min(eigenvalues))/6;
     end
    end
    ymax = ymax + abs(ymax-ymin)/6;
    if ymax <= ymin
        ymax = ymin + 1;
        ymin = ymin - 1;
    end
    axis(ahandle,[xmin xmax ymin ymax]);
    leg=legend(ahandle,legendStr{:},'Location','BestOutside');
    set(leg,'FontSize',9);
    set(get(ahandle,'XLabel'),'String','x','FontSize',11);
    set(fhandle,'Visible','on');
    set(fhandle,'HandleVisibility','off');
end


% --------------------------------------------------------------------
function Help_menu_Callback(hObject, eventdata, handles)
% hObject    handle to Help_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in workspace_button.
function workspace_button_Callback(hObject, eventdata, handles)
prompt={'Enter name for the indices:','Enter name for the eigenvalues:',...
        'Enter name for the errors:'};
title = 'Enter variable names';
lines = 1;
def = {'indices','eigenvalues','errors'};
answer = inputdlg(prompt,title,lines,def);
if ~isempty(answer)
  assignin('base',answer{1},handles.eigenvalueData(:,1));
  assignin('base',answer{2},handles.eigenvalueData(:,2));
  assignin('base',answer{3},handles.eigenvalueData(:,3));
end




% --- Executes on button press in eigenfunctionCompute_button.
function eigenfunctionCompute_button_Callback(hObject, eventdata, handles)
try
    hs=findobj(gcf,'Style','pushbutton','Enable','on','-or',...
        'style','text','-or','style','edit','-or','style','listbox',...
        '-or','style','popupmenu','-or','style','listbox');
    set(hs,'Enable','off')
    set(handles.info_edit,'ForeGroundColor','k');
    set(handles.info_edit,'String','Busy...','Enable','on')
    drawnow
    %which eigenfunction?:
    ind = get(handles.eigenfunctions_popupmenu,'Value');
    %which evaluation points?:
    ind2 = get(handles.xvalues_popupmenu,'Value');
    if ind2==1 %only in mesh points
       [handles.x,handles.y,handles.yp,nrm] = computeEigenfunction(handles.slp,handles.eigenvalueData(ind,1:2),handles.meshData);
       handles.evpoints=0;
    else
        prompt={'Enter the name of the workspace variable containing the eigenfunction evaluation points:'};
        title = 'Enter variable name';
        lines = 1;
        try
           answer = inputdlg(prompt,title,lines);
           x = evalin('base', answer{1});
        catch
           set(hs,'Enable','on') 
           set(handles.info_edit,'String','','Enable','on')
           errordlg('No valid variable name was entered') 
        end
        try
           [handles.x,handles.y,handles.yp,nrm] = computeEigenfunction(handles.slp,handles.eigenvalueData(ind,1:2),handles.meshData,x);
        catch
           set(hs,'Enable','on') 
           set(handles.info_edit,'String','','Enable','on')
           errordlg(lasterr,'Error','modal');  
           return;
        end
        handles.evpoints=1;
    end
    %display results:
    nrs= ceil(log10(max(abs(handles.x(end)),abs(handles.x(1)))));
    format=['%18.' num2str(16-nrs) 'f'];
    s1 = sprintf(format,handles.x(1));
    s2 = sprintf('%19.14f',handles.y(1));
    if abs(handles.yp(1))>10000
      s3 = sprintf('  %-23.10f',handles.yp(1));
    else
      s3 = sprintf('%19.14f',handles.yp(1));
    end
    resultsStr = {[s1,' ',s2,' ',s3]};
    for i=2:length(handles.x)
            s1 = sprintf(format,handles.x(i));
            s2 = sprintf('%19.14f',handles.y(i));
            if abs(handles.yp(i))>10000
                s3 = sprintf('  %-23.10f',handles.yp(i));
            else
                s3 = sprintf('%19.14f',handles.yp(i));
            end
            resultsStr = [resultsStr;...
                {[s1,' ',s2,' ',s3]}] ;
    end
    set(handles.eigenfunctionsData_listbox,'String',resultsStr);
    set(handles.eigenfunctionsData_listbox,'Value',[]);
    set(handles.info_edit,'String',['Computed Norm: ' num2str(nrm,6)]);

    %generate smooth plot:
    ind = get(handles.eigenfunctions_popupmenu,'Value');
    a = eval(get(handles.xmin_edit,'String'));
    b = eval(get(handles.xmax_edit,'String'));
    if isinf(a)
        a=handles.meshData.Infa;
    else
        a=handles.x(1);
    end
    if isinf(b)
        b=handles.meshData.Infb;
    else
        b=handles.x(end);
    end
    npts=150;
    x=linspace(a,b,npts);
    %add evaluation points in oscillatory region, i.e. the region where
    %(q-Ew)/p<0
    [p,q,w]=readInputFields(handles);
    E=handles.eigenvalueData(ind,2);
    pf=arrayfun(p,x); qf=arrayfun(q,x); wf=arrayfun(w,x);
    Z=(qf-E*wf)./pf;
    inds=find(Z<0);
    if any(diff(Z))
        while length(inds)<2 
            npts=npts*2;
            x=linspace(a,b,npts);
            pf=arrayfun(p,x);
            qf=arrayfun(q,x);
            wf=arrayfun(w,x);
            Z=(qf-E*wf)./pf;
            inds=find(Z<0);
        end
    else
        inds=[0 length(x)];
    end
    ind1=max(0,inds(1)-10);
    ind2=min(inds(end)+10,length(x)+1);
    index=handles.eigenvalueData(ind,1);
    if x(min(ind2,length(x)))-x(max(1,ind1))> (x(end)-x(1))/2
        %heuristic to find appropriate number of points
        npts=max(npts*3,floor(399*100./(1+exp(-index/400))-19900));
    else
        npts=min((handles.eigenvalueData(ind,1)+1)*100,2e4);
    end
    x=[x(1:ind1) linspace(x(ind1+1),x(ind2-1),npts) x(ind2:end) handles.x];
    [x,y,yp] = computeEigenfunction(handles.slp,handles.eigenvalueData(ind,1:2),handles.meshData,x);
    %should the eigenfunction be multiplied by -1?
    [m,ix]=max(handles.y);
    %find corresponding value in y:
    xm=handles.x(ix);
    inds=find(x>xm);
    if isempty(inds)
        inds=length(x);
    end
    if y(inds(1))*m<0
        y=-y; yp=-yp;
    end
    %create figure eigenfunction derivative
    fhandle = figure('IntegerHandle','off','NumberTitle','off',...
        'Visible','off','Name','Eigenfunction Derivative');
    pos=get(fhandle,'Position');
    set(fhandle,'Position',[pos(1:2) 650 250])
    movegui(fhandle,'north')
    pos=get(fhandle,'Position');
    ahandle = axes('Parent',fhandle);
    phandles=plot(handles.x,handles.yp,'Parent',ahandle);
    if length(handles.x)<25
        set(phandles,'Marker','o','LineStyle','none','Color','r','MarkerSize',5,'MarkerFaceColor','r')
    else
        set(phandles,'Marker','.','LineStyle','none','Color','r')
    end
    hold on
    phandles=plot(x,yp,'LineWidth',1,'Color','b','Parent',ahandle);
    hold off
    axis(ahandle,'tight');
    set(get(ahandle,'XLabel'),'String','x','FontSize',11);
    if handles.evpoints
        set(get(ahandle,'Title'),'String',['p(x)y''_{' num2str(handles.eigenvalueData(ind,1)) '}(x) and its evaluation in some user-specified evaluation points'])
    else
        set(get(ahandle,'Title'),'String',['p(x)y''_{' num2str(handles.eigenvalueData(ind,1)) '}(x) and its evaluation in the mesh points'])
    end
    set(fhandle,'Visible','on')

    %create figure eigenfunction 
    fhandle = figure('IntegerHandle','off','NumberTitle','off',...
        'Visible','off','Name','Eigenfunction');
    set(fhandle,'Position',[pos(1) pos(2)-340 650 250])
    %movegui(fhandle,'center')
    ahandle = axes('Parent',fhandle);
    phandles=plot(handles.x,handles.y,'Parent',ahandle);    
    if length(handles.x)<25
        set(phandles,'Marker','o','LineStyle','none','Color','r','MarkerSize',5,'MarkerFaceColor','r')
    else
        set(phandles,'Marker','.','LineStyle','none','Color','r')
    end
    hold on
    phandles2=plot(x,y,'LineWidth',1,'Parent',ahandle);
    axis(ahandle,'tight');
    set(get(ahandle,'XLabel'),'String','x','FontSize',11);
    if handles.evpoints
        set(get(ahandle,'Title'),'String',['Eigenfunction y_{' num2str(handles.eigenvalueData(ind,1)) '} and its evaluation in some user-specified evaluation points'])
    else
        set(get(ahandle,'Title'),'String',['Eigenfunction y_{' num2str(handles.eigenvalueData(ind,1)) '} and its evaluation in the mesh points'])
    end
    set(fhandle,'Visible','on')

    set([handles.eigenfunctionWorkspace_button],'Enable','on')
    set(hs,'Enable','on');
    guidata(hObject, handles);
catch
  set(hs,'Enable','on') 
  set(handles.info_edit,'String','','Enable','on')
  errordlg(lasterr,'Error','modal'); 
  return;
end




% --- Executes on button press in coeffPlot_button.
function coeffPlot_button_Callback(hObject, eventdata, handles)
try
    [p,q,w,xmin,xmax,a0,b0,a1,b1]=readInputFields(handles);
catch
    errordlg('The problem specification fields need to be filled out correctly.','Bad Input','modal');
    return;
end
try
h_figure = figure('Visible','off','Position',[360,500,530,450],'IntegerHandle','off','NumberTitle','off');
% create structure of handles
myhandles = guihandles(h_figure); 
myhandles.q = q; 
myhandles.p = p; 
myhandles.w = w; 
hpanel=uipanel('Title','x-range','Position',[.1,.02,.82,.12],'ForegroundColor','k');
hplot = uicontrol('Parent',hpanel,'Style','pushbutton','String','Plot',...
          'Position',[280,10,70,25],...
          'Callback',{@plotbutton_Callback});
hrefine = uicontrol('Parent',hpanel,'Style','pushbutton','String','Refine',...
          'Position',[355,10,70,25],...
          'Callback',{@refinebutton_Callback});
htext1=uicontrol('Parent',hpanel,'Style','text','String','xmin: ',...
          'Position',[8,15,60,15]);
htext2=uicontrol('Parent',hpanel,'Style','text','String','xmax: ',...
          'Position',[130,15,60,15]);
b=xmax;
if isinf(xmin)
    %compute a cutoff value
      start=min(0,xmax);
      vmin=q(start)/w(start);
      while isinf(vmin)
          start=start-0.0001;
          vmin=q(start)/w(start);
      end
      sums=0;
      xmin=start;
      h=0.1;
      it=0;
      vprev=vmin;
      while sums<35 && it<70
            it=it+1;
            h=h*1.01;
            xmin=xmin-h;
            va=q(xmin)/w(xmin);
            if vprev>va
                sums=0;
            end
            if va >= vmin 
               sums=sums+sqrt(va-vmin)*h;
            else
               sums=0;
               vmin=va;
            end
            vprev=va;
      end 
      if abs((q(xmin*2)/w(xmin*2)-q(xmin)/w(xmin))/xmin) < 3
          xmin=xmin*2;
      end
end     
hedit1=uicontrol('Parent',hpanel,'Style','edit','String',num2str(xmin),...
          'Position',[50,10,80,25],'HorizontalAlignment','right');   
if isinf(xmax)
      start=max(0,xmin);
      vmin=q(start)/w(start);
      while isinf(vmin)
          start=start+0.0001;
          vmin=q(start)/w(start);
      end
      sums=0;
      xmax=start;
      h=0.1;
      it=0;
      vbprev=vmin;
      while sums<35 && it<70
            it=it+1;
            h=h*1.01;
            xmax=xmax+h;
            vb=q(xmax)/w(xmax);
            if vbprev>vb
                sums=0;
                vmin=vb;
            end
            if vb >= vmin 
               sums=sums+sqrt(vb-vmin)*h;
            else
               sums=0;
               vmin=vb;
            end
            vbprev=vb;
      end 
      if q(xmax*2)/w(xmax*2)<(q(xmin)/w(xmin)) || isnan(q(xmin)/w(xmin))
          xmax=xmax*2;
      end
      if abs((q(xmax*2)/w(xmax*2)-q(xmax)/w(xmax))/xmax) < 3
          xmax=xmax*2;
      end
end        
hedit2=uicontrol('Parent',hpanel,'Style','edit','String',num2str(xmax),...
          'Position',[178,10,80,25],'HorizontalAlignment','right'); 
c=get(htext1,'BackgroundColor');
set(h_figure,'Color',c);
if ispc && isequal(get(hedit1,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hedit1,'BackgroundColor','white');
    set(hedit2,'BackgroundColor','white');
end      
ha = axes('Units','Pixels','Position',[55,115,430,305]); 
myhandles.ha=ha;
set([h_figure,ha,hplot,htext1,htext2,hedit1,hedit2],...
   'Units','normalized');
myhandles.npts=2000;
x=linspace(xmin,xmax,2001);
qf=arrayfun(q,x);pf=arrayfun(p,x);wf=arrayfun(w,x);
starti=1;
endi=length(x);
if all(pf) && all(wf) %schrodinger form
    [mi,ix]=min(qf);
    m=max(qf(ix:ix+floor((length(qf)-ix)/2)));% 
    tmp=max(abs(mi),abs(m));
    while starti<endi/2 && (isnan(qf(starti)) || isinf(qf(starti)) || (abs(qf(starti))-tmp)>10^4) 
        starti = starti+1;
    end
    while max(qf(starti+1:endi))+tmp < qf(starti) && min(qf(starti:endi)) ~= max(qf(starti:endi)) && ~(xmin==-xmax && ~isinf(b)) %morse1
        starti = starti+1;
    end
    while endi>1 && (isnan(qf(endi)) || isinf(qf(endi)) || (abs(qf(endi))-tmp)>10^4)
        endi = endi-1;
    end
    if starti>endi
        starti=1;
        endi=length(x);
    end
end
plot(x(starti:endi),pf(starti:endi),'b');
hold on
plot(x(starti:endi),qf(starti:endi),'r');
plot(x(starti:endi),wf(starti:endi),'g');
set(hedit2,'String',num2str(x(endi)));
set(hedit1,'String',num2str(x(starti)));
hold off
axis(ha,'tight');
leg=legend(ha,{'p(x)','q(x)','w(x)'},'Location','Best');
set(leg,'FontSize',9);
set(get(ha,'XLabel'),'String','x','FontSize',11);
set(h_figure,'Name','SLP Coefficient Functions')
movegui(h_figure,'center')
set(h_figure,'Visible','on');
guidata(h_figure,myhandles) 
catch EM
    errordlg(EM.message,'Error','modal');
end


function plotbutton_Callback(hObject,ed)
hs=findobj(get(hObject,'Parent'),'Style','edit');
try
  x1=eval(get(hs(1),'String'));
  x2=eval(get(hs(2),'String'));
catch
  errordlg('You must enter a numeric value','Bad Input','modal')  
end
myhandles=guidata(gcbf);
xmin=min(x1,x2); xmax=max(x1,x2);
h=(xmax-xmin)/myhandles.npts;
x=xmin:h:xmax;
plot(myhandles.ha,x,arrayfun(myhandles.p,x),'b');
hold on
plot(myhandles.ha,x,arrayfun(myhandles.q,x),'r');
plot(myhandles.ha,x,arrayfun(myhandles.w,x),'g');
hold off
axis(myhandles.ha,'tight');
leg=legend(myhandles.ha,{'p(x)','q(x)','w(x)'},'Location','Best');
set(leg,'FontSize',9);
set(get(myhandles.ha,'XLabel'),'String','x','FontSize',11);


function refinebutton_Callback(hObject,ed)
hs=findobj(get(hObject,'Parent'),'Style','edit');
try
  x1=eval(get(hs(1),'String'));
  x2=eval(get(hs(2),'String'));
catch
  errordlg('You must enter a numeric value','Bad Input','modal')  
end
myhandles=guidata(gcbf);
xmin=min(x1,x2); xmax=max(x1,x2);
myhandles.npts=myhandles.npts*2;
h=(xmax-xmin)/myhandles.npts;
x=xmin:h:xmax;
plot(myhandles.ha,x,arrayfun(myhandles.p,x),'b');
hold on
plot(myhandles.ha,x,arrayfun(myhandles.q,x),'r');
plot(myhandles.ha,x,arrayfun(myhandles.w,x),'g');
hold off
axis(myhandles.ha,'tight');
leg=legend(myhandles.ha,{'p(x)','q(x)','w(x)'},'Location','Best');
set(leg,'FontSize',9);
set(get(myhandles.ha,'XLabel'),'String','x','FontSize',11);
guidata(gcbf,myhandles) 


% --- Executes on button press in eigenfunctionWorkspace_button.
function eigenfunctionWorkspace_button_Callback(hObject, eventdata, handles)
prompt={'Enter name for the x-values:','Enter name for the y-values:',...
        'Enter name for the y''-values:'};
title = 'Enter variable names';
lines = 1;
def = {'x','y','yp'};
answer = inputdlg(prompt,title,lines,def);
if ~isempty(answer)
  assignin('base',answer{1},handles.x);
  assignin('base',answer{2},handles.y);
  assignin('base',answer{3},handles.yp);
end

function clearbutton_Callback(hObject, eventdata, handles)
clearFields(handles)


% --------------------------------------------------------------------
function helpMenuItem_Callback(hObject, eventdata, handles)
HelpPath = which('matslise.m');
[pathstr,name,ext] =fileparts(HelpPath);
HelpPath=fullfile(pathstr,'Helpfiles','matslise.html');
web(HelpPath);

% --------------------------------------------------------------------
function aboutMenuItem_Callback(hObject, eventdata, handles)
s1=['MATSLISE(2) is a MATLAB toolbox that provides a user-friendly interface to Sturm-Liouville computations. It features highly efficient algorithms to compute eigenvalues and eigenfunctions of regular and'...
' singular self-adjoint Sturm-Liouville problems. The numerical methods used are based on the so-called Constant (reference potential) Pertubation technique.'];
s2='This software package does not require any specific Matlab toolbox.';
s2b=[char(169) ' Veerle.Ledoux@UGent.be, Marnix.VanDaele@UGent.be'];
s3='DISCLAIMER';
s4='MATSLISE is freely available for non-commercial use. In no circumstances can the authors be held responsible for any deficiency, fault or other mishappening with regard to the use or performance of MATSLISE.';
s5='Last Update: June 2013';
msgbox(sprintf('%s\n\n%s\n\n\n%s\n%s\n\n\n%s\n%s',s1,s2,s3,s4,s2b,s5),'About MATSLISE')


function p_edit_Callback(hObject, eventdata, handles)

function p_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function q_edit_Callback(hObject, eventdata, handles)

function q_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function w_edit_Callback(hObject, eventdata, handles)

function w_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function xmin_edit_Callback(hObject, eventdata, handles)

function xmin_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function xmax_edit_Callback(hObject, eventdata, handles)

function xmax_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tol_edit_Callback(hObject, eventdata, handles)

function tol_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function a0_edit_Callback(hObject, eventdata, handles)

function a0_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function a1_edit_Callback(hObject, eventdata, handles)

function a1_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function b1_edit_Callback(hObject, eventdata, handles)

function b1_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function b0_edit_Callback(hObject, eventdata, handles)

function b0_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pmin_edit_Callback(hObject, eventdata, handles)

function info_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pmin_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pmax_edit_Callback(hObject, eventdata, handles)

function pmax_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function indices_menu_Callback(hObject, eventdata, handles)

function indices_menu_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');

end

% --- Executes on selection change in xvalues_popupmenu.
function xvalues_popupmenu_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function xvalues_popupmenu_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function eigenfunctions_popupmenu_Callback(hObject, eventdata, handles)


function eigenfunctions_popupmenu_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function eigenfunctionsData_listbox_Callback(hObject, eventdata, handles)


function eigenfunctionsData_listbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function eigenvalues_listbox_Callback(hObject, eventdata, handles)


function eigenvalues_listbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function classifInfo_edit_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function classifInfo_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function parameterEdit_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function parameterEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
