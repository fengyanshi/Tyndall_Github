
clear all
close all

case1='Baseline_HM';
fdir_hydro=['/Users/fyshi/OUTSIDE_Google/Users/XBeach_data/' case1 '/Xbeach_hydro_' case1 '/'];
fdir_waves=['/Users/fyshi/OUTSIDE_Google/Users/XBeach_data/' case1 '/bnd_windwaves_' case1 '/'];

bd=load('boundary.txt');


myVideo = VideoWriter('videoETA.mp4','MPEG-4');
myVideo.FrameRate = 25;  
myVideo.Quality = 100;
%vidHeight = 576; %this is the value in which it should reproduce
%vidWidth = 1024; %this is the value in which it should reproduce
open(myVideo);


fig=figure(1);
wid=10;
len=6;
set(fig,'units','inches','paperunits','inches','papersize', [wid len],'position',[2 2 wid len],'paperposition',[0 0 wid len]);
colormap jet


for k=1:length(bd)
lat=sprintf('%.4f',bd(k,1));
lon=sprintf('%.4f',bd(k,2));
fname=[fdir_hydro 'wl_XB' num2str(k) '_' lat '_' lon '_' case1 '.csv'];
Tbl=readtable(fname);
Ele=Tbl.wl;


x0=-86.0;
y0=29.75;

clf

plot([x0 x0+0.6],[y0 y0+0.6],'.r','MarkerSize',1)
plot_google_map('maptype','satellite','APIKey','AIzaSyBeu2oRBtLClpcm4i2VDIXltuzMAOY5yX4')
hold on

plot(bd(:,2),bd(:,1),'r*')
text(bd(:,2),bd(:,1),num2str(Ele(k)))

xlabel('Lon(deg) ','fontsize',12,'fontweight','bold');
ylabel('Lat(deg) ','fontsize',12,'fontweight','bold');

pause(0.1)

F = print('-RGBImage','-r300');
%J = imresize(F,[vidHeight vidWidth]);
mov(k).cdata = F;
writeVideo(myVideo,mov(k).cdata);

end
close(myVideo)



