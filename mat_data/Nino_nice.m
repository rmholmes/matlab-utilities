load('Nino3p4new.mat');
x = double(Nino3p4YRDEC);
y = double(Nino3p4F);

ENx = x(y>0.5);
ENy = y(y>0.5);

gaps = diff(ENx)>0.2;
ginds = find(gaps);

for ii=1:length(ginds)
    ENx((ginds(ii)+4):(end+3)) = ENx((ginds(ii)+1):end);
    ENy((ginds(ii)+4):(end+3)) = ENy((ginds(ii)+1):end);
    ENy(ginds(ii)+1) = 0.5;
    ENy(ginds(ii)+2) = NaN;
    ENy(ginds(ii)+3) = 0.5;
    dx = diff(x);
    dx = dx(1);
    ENx(ginds(ii)+1) = ENx(ginds(ii))+dx/2;
    ENx(ginds(ii)+2) = NaN;
    ENx(ginds(ii)+3) = ENx(ginds(ii)+4)-dx/2;
    ginds((ii+1):end) = ginds((ii+1):end)+3;
end
plot(ENx,ENy,'-r','LineWidth',2);

ENx = x(y<-0.5);
ENy = y(y<-0.5);

gaps = diff(ENx)>0.2;
ginds = find(gaps);

for ii=1:length(ginds)
    ENx((ginds(ii)+4):(end+3)) = ENx((ginds(ii)+1):end);
    ENy((ginds(ii)+4):(end+3)) = ENy((ginds(ii)+1):end);
    ENy(ginds(ii)+1) = -0.5;
    ENy(ginds(ii)+2) = NaN;
    ENy(ginds(ii)+3) = -0.5;
    dx = diff(x);
    dx = dx(1);
    ENx(ginds(ii)+1) = ENx(ginds(ii))+dx/2;
    ENx(ginds(ii)+2) = NaN;
    ENx(ginds(ii)+3) = ENx(ginds(ii)+4)-dx/2;
    ginds((ii+1):end) = ginds((ii+1):end)+3;
end
plot(ENx,ENy,'-b','LineWidth',2);