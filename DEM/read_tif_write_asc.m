


clear all
fdir='/Users/fyshi/OUTSIDE_Google/GITHUB_M3/Tyndall/DEM/';
fcur='/Users/fyshi/OUTSIDE_Google/GITHUB_M3/Tyndall/DEM/';
eval(['cd ' fdir])
uiopen('exportImage.tiff');

% break

dep=-flipud(exportImage.tiff);

[n m]=size(dep);

%            LatitudeLimits: [29.8334293090301 30.3395406909699]
%            LongitudeLimits: [-85.98898 -85.29685]
%                 RasterSize: [5466 7475]

x0=-78.0-25./60.;
y0=34.0;

dx=1.0/3600.0/9.0;

x=x0+[0:m-1]*dx;
y=(y0-(n-1)*dx)+[0:n-1]*dx;

% 1/3 s
sk=3;

dep1=dep(1:sk:end,1:sk:end);
x1=x(1:sk:end);
y1=y(1:sk:end);
[n1 m1]=size(dep1);
pcolor(x1,y1,-dep1),shading flat

% header
ncol=m1;
mrow=n1;
xll=x1(1);
yll=y1(1);
cells=1./3600./3.;

eval(['cd ' fcur])

writeout=-flipud(dep1);
filename=['wilmingtonnc.asc'];
fid = fopen(filename,'wt');
fprintf(fid, '%s\t%d\n', 'ncols',ncol);
fprintf(fid, '%s\t%d\n', 'nrows',mrow);
fprintf(fid, '%s\t%16.8f\n', 'xllcorner',xll);
fprintf(fid, '%s\t%16.8f\n', 'yllcorner',yll);
fprintf(fid, '%s\t%16.8f\n', 'cellsize',cells);
fprintf(fid, '%s\t%d\n', 'NODATA_value',-99999);
fclose(fid);
dlmwrite(filename,writeout,'delimiter',' ','precision',['%10.',num2str(4),'f'],'-append');


