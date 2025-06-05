
clear all
close all

case1='Baseline_CHM';
fdir_hydro=['/Users/fyshi/OUTSIDE_Google/Users/XBeach_data/' case1 '/Xbeach_hydro_' case1 '/'];
fdir_waves=['/Users/fyshi/OUTSIDE_Google/Users/XBeach_data/' case1 '/bnd_windwaves_' case1 '/'];

bd=load('boundary.txt');


for k=1:length(bd)
lat=sprintf('%.4f',bd(k,1));
lon=sprintf('%.4f',bd(k,2));
fname=[fdir_waves 'windwaves_XB' num2str(k) '_' lat '_' lon '_' case1 '.csv'];
Tbl=readtable(fname);
Hs1(:,k)=Tbl.Hs;
Tp1(:,k)=Tbl.tp;
Dm1(:,k)=Tbl.mwd;
Ds1(:,k)=Tbl.ws;
Time1=Tbl.date_time;
end

% -------------
nstart=175;
nend=250;
% ------------

k_bc=20;

Time=Time1(nstart:nend);
Hs = Hs1(nstart:nend,k_bc) ;
Tp= Tp1(nstart:nend,k_bc);
Dm=Dm1(nstart:nend,k_bc);  
Ds=Ds1(nstart:nend,k_bc); 

set(0,'DefaultFigureColormap',feval('jet'));
figure(1)
subplot(211)
plot(Time1,Hs1)
%text(170,Hs1(169),fpoints{k})
hold on
plot(Time,Hs,'r','LineWidth',2)

% write
for k=1:length(Time)
   nyear=num2str(year(Time(k)),'%.4d');
   nmonth=num2str(month(Time(k)),'%.2d');
   ndate=num2str(day(Time(k)),'%.2d');
   nhour=num2str(hour(Time(k)),'%.2d');
   nminute=num2str(minute(Time(k)),'%.2d');
   yeartime=[nyear nmonth ndate '.' nhour nminute '00'];
   hh=num2str(Hs(k),'%4.1f');
   wp=num2str(Tp(k),'%4.1f');
   wa=num2str(Dm(k),'%4.1f');
   ws=num2str(Ds(k),'%4.1f');
   tot=[yeartime ' ' hh ' ' wp ' ' wa ' ' ws];
wavenew{k}=tot;
end

%eval(['mkdir ' 'data/']);
%eval(['mkdir ' 'data/' fcase]);

fid=fopen(['wave.txt'],'w','n');
fprintf(fid,'%s\n','TPAR');
for k=1:length(wavenew)
%fprintf(fid, '%s\n', strtrim(wavenew(k,:)));
fprintf(fid, '%s\n', wavenew{k});
end
fclose(fid);


