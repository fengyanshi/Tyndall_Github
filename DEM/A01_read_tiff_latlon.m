%% Read image file using readgeoraster

clear all

fname = "exportImage.tiff";
[A,R] = readgeoraster(fname);
  %          LatitudeLimits: [29.8405580908311 30.3621319091689]
   %         LongitudeLimits: [-86.10157 -85.3462]
    %             RasterSize: [5633 8158]

dep=-flipud(A);

x0=-86.10157;
y0=29.8405580908311;

dx=1.0/3600.0/3.0;

[n m]=size(dep);

x=x0+[0:m-1]*dx;
y=y0+[0:n-1]*dx;

% from NAVD88(0.179) to MTL(0.214) +0.035

depp=double(dep)+0.035;

%save -ASCII dep_8158x5633.txt depp 
%save -ASCII lon_8158x5633.txt x
%save -ASCII lat_8158x5633.txt y

bd=load('boundary.txt');

% contours
[C,h]=contour(x,y,depp,[0 0]);
lat=C(2,1:end);
lon=C(1,1:end);
lat(lat>30.35)=NaN;
lat(lat<29.85)=NaN;
lon(lon<-86.1)=NaN;
lon(lon>-85.3)=NaN;

writeout(:,1)=lon;
writeout(:,2)=lat;

save -ASCII coast_xy.txt writeout
figure
plot(lon,lat)
hold on
for k=1:length(bd)
plot(bd(k,2),bd(k,1),'r*')
text(bd(k,2),bd(k,1),num2str(k),'Color','w')
end
xlabel('Lon(deg) ','fontsize',12,'fontweight','bold');
ylabel('Lat(deg) ','fontsize',12,'fontweight','bold');
print('-djpeg', 'coast_data.jpg')

figure
pcolor(x(1:10:end),y(1:10:end),-dep(1:10:end,1:10:end)),shading flat
demcmap(dep)
colorbar
hold on
for k=1:length(bd)
plot(bd(k,2),bd(k,1),'r*')
text(bd(k,2),bd(k,1),num2str(k),'Color','w')
end
xlabel('Lon(deg) ','fontsize',12,'fontweight','bold');
ylabel('Lat(deg) ','fontsize',12,'fontweight','bold');
print('-djpeg', 'map_data.jpg')



%caxis([-15 15])
