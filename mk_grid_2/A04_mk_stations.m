clear all
x=load('lon_curv.txt');
y=load('lat_curv.txt');
dep1=load('dep_curv_meas_data.txt');


eval(['cd ' '/Users/fengyanshi/WORK/work/ROCKY_BEACH/bathy_data/elevation_products/']);

figure(1)
wid=8;
len=7.;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid+1 len+1],'paperposition',[0 0 wid len]);

clf
contour(x,y,dep1,[0 0],'r','LineWidth',2)
hold on

%chinarock
load zmsl_ChinaRock_gridded.mat
lat=bathymetry.latitude;
lon=bathymetry.longitude;

x1=lon(1,1);
x2=lon(1,end);
x3=lon(end,end);
x4=lon(end,1);
y1=lat(1,1);
y2=lat(1,end);
y3=lat(end,end);
y4=lat(end,1);
plot([x1 x2 x3 x4 x1],[y1 y2 y3 y4 y1])

%asilomar

load zmsl_Asilomar_gridded.mat

lat=bathymetry.latitude;
lon=bathymetry.longitude;

x1=lon(1,1);
x2=lon(1,end);
x3=lon(end,end);
x4=lon(end,1);
y1=lat(1,1);
y2=lat(1,end);
y3=lat(end,end);
y4=lat(end,1);
plot([x1 x2 x3 x4 x1],[y1 y2 y3 y4 y1])

% Hopkin
load zmsl_Hopkins_gridded.mat

lat=bathymetry.latitude;
lon=bathymetry.longitude;

x1=lon(1,1);
x2=lon(1,end);
x3=lon(end,end);
x4=lon(end,1);
y1=lat(1,1);
y2=lat(1,end);
y3=lat(end,end);
y4=lat(end,1);
plot([x1 x2 x3 x4 x1],[y1 y2 y3 y4 y1])

x1=x(1,1);
x2=x(1,end);
x3=x(end,end);
x4=x(end,1);
y1=y(1,1);
y2=y(1,end);
y3=y(end,end);
y4=y(end,1);
plot([x1 x2 x3 x4 x1],[y1 y2 y3 y4 y1])

eval(['cd ' '/Users/fengyanshi/WORK/work/ROCKY_BEACH/NearCoM/make_grid_bathy']);

x_sta=[40 45 50 55 60 65 70 75  40  45  50  55  60  65  70  75    121 120 119 118 116 115];
y_sta=[66 68 69 71 73 74 76 78  116 117 118 119 121 122 123 125   160 165 170 175 180 185];

for k=1:length(x_sta)
plot(x(y_sta(k),x_sta(k)),y(y_sta(k),x_sta(k)),'ko')
end

writeout(:,1)=x_sta;
writeout(:,2)=y_sta;

fid = fopen('stations.txt', 'wt');
  fprintf(fid, ['%5d','%5d', '\n'],writeout');
fclose(fid);

print -djpeg100 stations.jpg



