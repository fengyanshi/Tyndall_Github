
clear all
close all

bd=load('boundary.txt');

%set(0,'DefaultFigureColormap',feval('jet'));
figure(1)
wid=10;
len=6;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[2 2 wid len],'paperposition',[0 0 wid len]);

colormap parula

x0=-85.98898;
y0=29.8334293090301;

plot([x0 x0+0.5],[y0 y0+0.5],'.r','MarkerSize',1)
plot_google_map('maptype','satellite','APIKey','AIzaSyBeu2oRBtLClpcm4i2VDIXltuzMAOY5yX4')
hold on

for k=1:length(bd)
plot(bd(k,2),bd(k,1),'r*')
text(bd(k,2),bd(k,1),num2str(k),'Color','w')
end

xlabel('Lon(deg) ','fontsize',12,'fontweight','bold');
ylabel('Lat(deg) ','fontsize',12,'fontweight','bold');





