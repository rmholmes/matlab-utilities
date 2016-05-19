function shade_anomaly(x,y,hcolor,lcolor,offset)
% shade_anomaly(x,y,hcolor,lcolor,offset)
%
% Version 1.0
%
% This function will shade the area above and below the zero mark different
% colors for anomaly plots. You can also create an offset that will shade
% only above a certain value. This script is based off symmetry, but can be
% edited to allow for different thresholds above and below the x-axis.
% --- INPUT ---
% x: x-values (If you only have y-values, just input those in place of x.)
% y: y-values (Maybe only input)
% hcolor/lcolor: Color of shading above and below x-axis (i.e. 'r')
%                (Default colors - red for upper, blue for lower)
% offset: The offset for shading above a certain value only, default = 0
%
% For an example of the plot, do not enter any inputs.
%
% Created by Alec
% asb05j@fsu.edu with any questions.

% Create a new figure window
% $$$ figure

% Based on the number of input arguments, make assumptions.
if nargin==0
    y=[rand(1,10),-1*rand(1,10),rand(1,10),-1*rand(1,10),rand(1,10),-1*rand(1,10),rand(1,10),-1*rand(1,10)];
    y=y(:);
    x=1:length(y);
    hcolor='r';
    lcolor='b';
    offset=.3;
elseif nargin==1
    y=x;
    x=1:length(y);
    hcolor='r';
    lcolor='b';
    offset=0;
elseif nargin==2
    hcolor='r';
    lcolor='b';
    offset=0;
elseif nargin==3
    error('You need to enter a lower color as well.')
elseif nargin==4
    offset=0;
elseif nargin==5
    disp('Nothing to assume.')
else
    error('You entered too many arguments')
end

% Make sure the offset is a positive number
offset=abs(offset);

%% TOP HALF

% Look at the top half of the data by setting y-values below zero to zero.
top=y;
top=top-offset;
top(top<=0)=0;
% Set up an x-vector for the top half
xx=x;

% If a value falls above the threshold set z to infinite otherwis 0.
z=top;
z(z>0)=Inf;
z(z<=0)=0;

% Find where values change from 0 to Inf. This finds points where the the
% plot crosses the threshold value. This makes the shading line up with the
% plot correctly. 
z=find(abs(diff(z))==Inf);

% Find the x-value when the y-value is equal to zero (threshold value). 
for i=1:length(z)
    m=(y(z(i)+1)-y(z(i)))/(x(z(i)+1)-x(z(i)));
    addon=((offset-y(z(i)))/m)+x(z(i));
    xx(length(xx)+1)=addon;
    top(length(top)+1)=0;
end

% Resort the data so it's in ascending order. 
q=sortrows([xx(:) top],1);

% Reset the x and y vectors, and remove values at the ends that are bad.
xx=q(:,1);
xx=xx(:);
xx(xx>max(x)+length(z))=NaN;
xx(xx<min(x))=NaN;
top(isnan(xx))=NaN;
top=q(:,2);
top=top(:);

% Create the offset array. This is because the area command works by
% stacking areas on one another.
up=offset.*ones(length(xx),1);

% Clear z to reuse.
clear z

%% BOTTOM HALF

% Look at the bottom half of the data by setting y-values above zero to zero.
bottom=y;
bottom=bottom+offset;
bottom(bottom>=0)=0;
% Set up an x-vector for the top half
xxx=x;

% If a value falls below the threshold set z to infinite otherwis 0.
z=bottom;
z(z<0)=-Inf;
z(z>=0)=0;

% Find where values change from 0 to -Inf. This finds points where the the
% plot crosses the threshold value. This makes the shading line up with the
% plot correctly. 
z=find(abs(diff(z))==Inf);

% Find the x-value when the y-value is equal to zero (threshold value). 
for i=1:length(z)
    m=(y(z(i)+1)-y(z(i)))/(x(z(i)+1)-x(z(i)));
    addon=((-offset-y(z(i)))/m)+x(z(i));
    xxx(length(xxx)+1)=addon;
    bottom(length(bottom)+1)=0;
end

% Resort the data so it's in ascending order. 
q=sortrows([xxx(:) bottom],1);

% Reset the x and y vectors, and remove values at the ends that are bad.
xxx=q(:,1);
xxx=xxx(:);
xxx(xxx>max(x)+length(z))=NaN;
xxx(xxx<min(x))=NaN;
bottom(isnan(xxx))=NaN;
bottom=q(:,2);
bottom=bottom(:);

% Create the offset array. This is because the area command works by
% stacking areas on one another.
down=-offset.*ones(length(xxx),1);

%% PLOT
% Plot the upper area.
h=area(xx,[up, top]);
set(h(1),'FaceColor','none','LineStyle','none');
set(h(2),'FaceColor',hcolor,'LineStyle','none');
hold on
% Plot the lower area.
g=area(xxx,[down, bottom]);
set(g(1),'FaceColor','none','LineStyle','none');
set(g(2),'FaceColor',lcolor,'LineStyle','none');

% Plot the actual data.
plot(x,y,'k','LineWidth',2)

% Turn grid on and make solid line at 0. 
grid on
set(gca,'Layer','Top')
hline(0,'k')