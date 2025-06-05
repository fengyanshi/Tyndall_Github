
clear all
close all

case1='Baseline_CHM';
fdir_hydro=['/Users/fyshi/OUTSIDE_Google/Users/XBeach_data/' case1 '/Xbeach_hydro_' case1 '/'];
fdir_waves=['/Users/fyshi/OUTSIDE_Google/Users/XBeach_data/' case1 '/bnd_windwaves_' case1 '/'];

bd=load('boundary.txt');


myVideo = VideoWriter('videoETA.mp4','MPEG-4');
myVideo.FrameRate = 25;  
myVideo.Quality = 100;
%vidHeight = 576; %this is the value in which it should reproduce
%vidWidth = 1024; %this is the value in which it should reproduce
open(myVideo);


for k=1:length(bd)
lat=sprintf('%.4f',bd(k,1));
lon=sprintf('%.4f',bd(k,2));
fname=[fdir_hydro 'wl_XB' num2str(k) '_' lat '_' lon '_' case1 '.csv'];
Tbl=readtable(fname);
Ele(:,k)=Tbl.wl;
Time=Tbl.date_time;
end

% interp
fdir='../mk_grid_2/';
x=load([fdir 'lon_curv_correct.txt']);
y=load([fdir 'lat_curv_correct.txt']);

sk=10;
xx=x(1:sk:end,1:sk:end);
yy=y(1:sk:end,1:sk:end);

xb=bd(:,2);
yb=bd(:,1);

for k=1:size(Ele,1) % time
ele=Ele(k,:)';
zz(:,:,k)=griddata(xb,yb,ele,xx,yy,'linear');
end


fig=figure(1);
wid=10;
len=12;
set(fig,'units','inches','paperunits','inches','papersize', [wid len],'position',[2 2 wid len],'paperposition',[0 0 wid len]);
colormap jet

% movie

time=[0:size(Ele,1)];

for k=1:10:size(Ele,1)

clf

subplot(211)
x0=-86.0;
y0=29.75;
plot([x0 x0+0.6],[y0 y0+0.6],'.r','MarkerSize',1)
plot_google_map('maptype','satellite','APIKey','AIzaSyBeu2oRBtLClpcm4i2VDIXltuzMAOY5yX4')
hold on

pcolor(xx,yy,zz(:,:,k)),shading flat
caxis([0 2.5])
colorbar

xlabel('Lon(deg) ','fontsize',12,'fontweight','bold');
ylabel('Lat(deg) ','fontsize',12,'fontweight','bold');

axis([-85.70 -85.48 29.93 30.1])

title(string(Time(k)))

subplot(212)
plot(time(1:k),squeeze(zz(5,5,1:k)))
axis([0 length(zz) 0 2.5])
grid
xlabel('time (frame) ')
ylabel('elevation(m)')

pause(0.1)

F = print('-RGBImage','-r300');
%J = imresize(F,[vidHeight vidWidth]);
mov(k).cdata = F;
writeVideo(myVideo,mov(k).cdata);

end

close(myVideo)



