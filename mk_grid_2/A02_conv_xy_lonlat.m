clear all

lon2=load('x_curv.txt');
lat2=load('y_curv.txt');
dep2=load('dep_curv.txt');

% different from the Norfork work, the grid generation provides a up-side-down array so correct it

lon1=lon2';
lat1=lat2';
dep1=dep2';

lon=flipud(lon1);
lat=flipud(lat1);
dep=flipud(dep1);

Lon=lon(1:168,1:256);
Lat=lat(1:168,1:256);

save -ASCII lon_curv_correct.txt Lon
save -ASCII lat_curv_correct.txt Lat


lon_m=min(min(lon));
lat_m=min(min(lat));

R=6371000.0;

x=(lon-lon_m)*R*pi/180.0*cos(pi/180.0*30.075);
y=(lat-lat_m)*R*pi/180.0;

save -ASCII xx_curv.txt x
save -ASCII yy_curv.txt y
save -ASCII dep_curv_correct.txt dep

[dxx dxy]=gradient(x);
[dyx dyy]=gradient(y);

dx=sqrt(dxx.^2+dyx.^2);
dy=sqrt(dxy.^2+dyy.^2);

ds=sqrt(dx.^2+dy.^2);

dx_min=min(min(dx));
dy_min=min(min(dy));

figure(1)
wid=8;
len=12;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid+1 len+1],'paperposition',[0 0 wid len]);
clf
subplot(211)
pcolor(lon,lat,dx),shading interp;
ch=colorbar;

xlabel('Lon (Deg)')
ylabel('Lat (Deg)')
grid
title(['grid size dx (m), Min:' num2str(dx_min)])

subplot(212)
pcolor(lon,lat,dy),shading interp;
ch=colorbar;

title(['grid size dy (m), Min:' num2str(dy_min)])

print -djpeg100 gridsize.jpg
return

figure(2)
clf
line(X,Y)
line(X',Y')
grid
xlabel('EAST (m)')
ylabel('NORTH (m)')
print -djpeg100 grid_in_xy.jpg

