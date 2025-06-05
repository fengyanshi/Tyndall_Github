clear all
fdir='/Users/fyshi/OUTSIDE_Google/GITHUB_M3/Tyndall/DEM/';
corr=load('../grid_2/grid_171x257_xy.txt');
m=171;  
n=257;  
     for i=1:m;                                 
      for j=1:n;                                          
        x(j,i)=corr((j-1)*m+i,1)*1.0;               
        y(j,i)=corr((j-1)*m+i,2)*1.0;               
      end                                       
   end  

bathy=load([fdir 'dep_8158x5633.txt']);
bathy=-bathy;

lon=load([fdir 'lon_8158x5633.txt']);
lat=load([fdir 'lat_8158x5633.txt']);
[X Y]=meshgrid(lon,lat);
[Lon Lat]=meshgrid(lon,lat);


% interp
bathy_curv=griddata(X,Y,bathy,x,y);

%bathy_curv=load('dep_curv.txt');

figure(1)
wid=8;
len=6;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid+1 len+1],'paperposition',[0 0 wid len]);
clf
pcolor(x,y,bathy_curv),shading interp;
hold on
sk=1;
xx=x;
yy=y;

xx(bathy_curv>1.0)=NaN;

line(xx(1:sk:end,1:sk:end),yy(1:sk:end,1:sk:end),'Color','k')
line(xx(1:sk:end,1:sk:end)',yy(1:sk:end,1:sk:end)','Color','k')
xlabel('Lon (Deg)')
ylabel('Lat (Deg)')
grid
demcmap(bathy_curv);
print -djpeg100 grid_bathy_xy.jpg

save -ASCII dep_curv.txt bathy_curv
save -ASCII x_curv.txt x
save -ASCII y_curv.txt y

% latlon

lat_curv=griddata(X,Y,Lat,x,y);
lon_curv=griddata(X,Y,Lon,x,y);

%lat_curv=load('lat_curv.txt');
%lon_curv=load('lon_curv.txt');

figure(2)
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
demcmap(bathy_curv);
print -djpeg100 grid_bathy_latlon.jpg

save -ASCII lon_curv.txt lon_curv
save -ASCII lat_curv.txt lat_curv


