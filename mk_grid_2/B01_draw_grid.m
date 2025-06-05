clear all

lon_curv1=load('lon_curv.txt');
lat_curv1=load('lat_curv.txt');
lon_curv=lon_curv1(1:end-1,4:end);
lat_curv=lat_curv1(1:end-1,4:end);
bathy_curv=-load('dep_circ.txt');

save -ASCII lon_circ.txt lon_curv
save -ASCII lat_circ.txt lat_curv

ax1=[-76.345 -76.25 36.90 36.98];

figure(1)
wid=8;
len=6;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid+1 len+1],'paperposition',[0 0 wid len]);
clf
pcolor(lon_curv,lat_curv,bathy_curv),shading interp;
hold on
sk=1;
lonc=lon_curv;
latc=lat_curv;

lonc(bathy_curv>1.0)=NaN;

line(lonc(1:sk:end,1:sk:end),latc(1:sk:end,1:sk:end),'Color','k')
line(lonc(1:sk:end,1:sk:end)',latc(1:sk:end,1:sk:end)','Color','k')
xlabel('lon (deg)')
ylabel('lat (deg)')
grid
%plot([ax1(1) ax1(2) ax1(2) ax1(1) ax1(1)],[ax1(3) ax1(3) ax1(4) ax1(4) ax1(3)],'w-','LineWidth',1)
demcmap(bathy_curv);
print -djpeg100 grid_bathy_latlon_B.jpg
print -depsc2 grid_bathy_latlon_B.eps

% zoom
wid=8;
len=6;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid+1 len+1],'paperposition',[0 0 wid len]);
axis([-85.66 -85.54 30.02 30.14])
print -djpeg100 grid_bathy_latlon_zoom.jpg
print -depsc2 grid_bathy_latlon_zoom.eps



