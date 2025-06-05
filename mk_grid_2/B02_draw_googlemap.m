clear all

lon_curv=load('lon_circ.txt');
lat_curv=load('lat_circ.txt');
bathy_curv=load('dep_circ.txt');

bathy=-bathy_curv;
bathy(bathy>0)=NaN;

bd=load('boundary.txt');

x0=-85.98898;
y0=29.8334293090301;

figure(1)
clf
wid=8;
len=6;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid+1 len+1],'paperposition',[0 0 wid len]);
colormap jet
plot([x0 x0+0.5],[y0 y0+0.5],'.r','MarkerSize',1)
plot_google_map('maptype','satellite','APIKey','AIzaSyBeu2oRBtLClpcm4i2VDIXltuzMAOY5yX4')
hold on

xlabel('Lon(deg) ','fontsize',12,'fontweight','bold');
ylabel('Lat(deg) ','fontsize',12,'fontweight','bold');

b=pcolor(lon_curv,lat_curv,bathy),shading interp;
set(b,'FaceAlpha',0.8);
set(b,'AlphaData',~isnan(bathy))
colorbar
caxis([-30 10])
axis([-85.95 -85.4 29.85 30.3])
grid
for k=1:length(bd)
plot(bd(k,2),bd(k,1),'r*')
text(bd(k,2),bd(k,1),num2str(k),'Color','w')
end

print -djpeg100 dep_googlemap.jpg
