
clear all
close all

case1='Baseline_CHM';
fdir_hydro=['/Users/fyshi/OUTSIDE_Google/Users/XBeach_data/' case1 '/Xbeach_hydro_' case1 '/'];
fdir_waves=['/Users/fyshi/OUTSIDE_Google/Users/XBeach_data/' case1 '/bnd_windwaves_' case1 '/'];

foutput=['/Users/fengyanshi/OUTSIDE_Google/GITHUB/Norfolk/Tyndall/Coupling_data/' case1 '/'];

%eval(['mkdir ' foutput case1]);

bd=load('boundary.txt');

fdir='../mk_grid_2/';
x=load([fdir 'lon_curv_correct.txt']);
y=load([fdir 'lat_curv_correct.txt']);


for k=1:length(bd)
lat=sprintf('%.4f',bd(k,1));
lon=sprintf('%.4f',bd(k,2));
fname=[fdir_hydro 'wl_XB' num2str(k) '_' lat '_' lon '_' case1 '.csv'];
Tbl=readtable(fname);
Ele(:,k)=Tbl.wl;
Time=Tbl.date_time;
end

% interp


sk=1;
xx=x(1,1:sk:end);
yy=y(1,1:sk:end);

xb=bd(:,2);
yb=bd(:,1);

for k=1:size(Ele,1) % time
ele=Ele(k,:)';
zz(:,:,k)=griddata(xb,yb,ele,xx,yy,'linear');
end

% write --------------
time2=[0:length(Time)-1]*20*60; % second, the data has a 20 min step
Eta2=squeeze(zz(1,:,:));

% -------------
%coupling_filename={'coupling.txt'};
logfile='log.txt';

nstart=175;
nend=250;

time=([nstart:nend]-nstart)*1200.0;
southpoints=[1:256];
northpoints=[1:256];

for sp=1:length(southpoints)
for ti=1:length(time)
southfine(sp,1,ti)=0.0;  % u
southfine(sp,2,ti)=0.0;  % v
southfine(sp,3,ti)=Eta2(sp,ti+nstart-1);  %eta
end
end

%southfine=[];
northfine=[];
westfine=[];
eastfine=[];


filename=['coupline_file.txt'];
FIN = fopen(filename,'w');           

Finfo = fopen(logfile,'w');  
% log file
fprintf(Finfo,'coupling data\nboundary info: num of points, start point');

if isempty(eastfine)==0
fprintf(Finfo,'\nEAST\n\t%d\t\t%d',size(eastfine,1),1);
else
fprintf(Finfo,'\nEAST\n\t%d\t\t%d',-1,1);
end

if isempty(westfine)==0
fprintf(Finfo,'\nWEST\n\t%d\t\t%d',size(westfine,1),1);
else
fprintf(Finfo,'\nWEST\n\t%d\t\t%d',-1,1);
end

if isempty(southfine)==0
fprintf(Finfo,'\nSOUTH\n\t%d\t\t%d',size(southfine,1),1);
else
fprintf(Finfo,'\nSOUTH\n\t%d\t\t%d',-1,1);
end

if isempty(northfine)==0
fprintf(Finfo,'\nNORTH\n\t%d\t\t%d',size(northfine,1),1);
else
fprintf(Finfo,'\nNORTH\n\t%d\t\t%d',-1,1);
end

% end log file

%%
fprintf(FIN,'coupling data\nboundary info: num of points, start point');

if isempty(eastfine)==0
fprintf(FIN,'\nEAST\n\t%d\t\t%d',size(eastfine,1),1);
else
fprintf(FIN,'\nEAST\n\t%d\t\t%d',-1,1);
end

if isempty(westfine)==0
fprintf(FIN,'\nWEST\n\t%d\t\t%d',size(westfine,1),1);
else
fprintf(FIN,'\nWEST\n\t%d\t\t%d',-1,1);
end

if isempty(southfine)==0
fprintf(FIN,'\nSOUTH\n\t%d\t\t%d',size(southfine,1),1);
else
fprintf(FIN,'\nSOUTH\n\t%d\t\t%d',-1,1);
end

if isempty(northfine)==0
fprintf(FIN,'\nNORTH\n\t%d\t\t%d',size(northfine,1),1);
else
fprintf(FIN,'\nNORTH\n\t%d\t\t%d',-1,1);
end

%%


%fprintf(FIN,'coupling data\nboundary info: num of points, start point');
%fprintf(FIN,'\nEAST\n\t%d\t\t%d',size(eastfine,1),1);
%fprintf(FIN,'\nWEST\n\t%d\t\t%d',size(westfine,1),1);
%fprintf(FIN,'\nSOUTH\n\t%d\t\t%d',-1,1);
%fprintf(FIN,'\nNORTH\n\t%d\t\t%d',-1,1);

fprintf(FIN,'\nTIME SERIES');
for t = 1:length(time)
    disp(sprintf('Writing Time Step No. %d    of   %d',t,length(time) ))
    fprintf(FIN,'\n\t%f',time(t));
    printside(FIN,'EAST',eastfine,t)
    printside(FIN,'WEST',westfine,t)
    printside(FIN,'SOUTH',southfine,t)
    printside(FIN,'NORTH',northfine,t)
end
fclose(FIN);
disp('Finished!')
%disp(sprintf('NOTE: This coupling file starts at time = %d sec',sample(start,1)))
%disp(sprintf('      ends at time = %d sec',sample(stop,1)))

%fprintf(Finfo,['\n' sprintf('NOTE: This coupling file starts at time = %d sec',sample(start,1))]);
%fprintf(Finfo,['\n' sprintf('      ends at time = %d sec',sample(stop,1))]);
fprintf(Finfo,['\n' 'Total time steps: ' num2str(length(time))]);
fprintf(Finfo,['\n' 'Time interval is AROUND: ' num2str(time(2)-time(1))]);

%clear eta_fine u_fine v_fine eastfine westfine southfine northfine

fclose(Finfo);


% -------------






