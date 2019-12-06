% **** Read SEAMAP Zooplankton Data ****
% Written by: Taylor Shropshire
% Date: 12/5/18
% *************************************** %
clear;clc;close all

% Paths 
file='/nexsan/people/taylor/projects/mit_gcm/zooplankton_data/SEAMAP_zoo_displacement_vol_1982_2017.xls';
load('/nexsan/people/taylor/projects/NPZ/data/matlab/eastcoastmask.mat'); clear lon lat			% coastline 

% Original file for reference (not used in code)
% 6: Lat                 % 7: Lon
% 8: Date                % 12: bottom depth 
% 13: volume filtered    % 14: volume displaced
% 15: min depth          % 16: max depth  
% 25: mesh size

% Data File 
% 4: lat                 % 5: lon 
% 6: date                % 7: gmt 
% 10: bottom depth       % 11: volume filtered 
% 12: volume displaced   % 13: min depth
% 14: max depth

% Read Data
display('Reading Data');tic 
data_all=importdata([file]);
text_data=data_all.textdata;
data=data_all.data;
ind_data=[4, 5, 6, 7, 10, 11, 12, 13, 14];
data=data(2:end,ind_data);toc

lat_all=data(:,1);
lon_all=data(:,2);
date_all=data(:,3);
time_all=data(:,4);
bot_depth_all=data(:,5);
vol_filt_all=data(:,6);
vol_disp_all=data(:,7);
depth_lim_all(:,1)=data(:,8);
depth_lim_all(:,2)=data(:,9);
mesh_all=str2num(cell2mat(text_data(2:end,end)));
zoo_all=vol_disp_all./vol_filt_all;

% Scale date
date_sample=datenum(1982,6,16);
sf_date=date_sample-date_all(1);
date_all=date_all+sf_date;

% ------- Quality Control ------ % 
% A few data points with weird lat lon format
lat_all2=str2num(cell2mat(text_data(2:end,6)));
lon_all2=str2num(cell2mat(text_data(2:end,7)));
lon_all2=lon_all2/100*-1;
lat_all2=lat_all2/100;

lon_all3=[]; lat_all3=[];
for ii=1:length(lon_all2)
n1=num2str(lon_all2(ii));
ind1=find(n1=='.');
if ~isempty(ind1)
dec=str2num(n1(ind1:end));
dec=dec/0.6;
dec=num2str(dec);
ind2=find(dec=='.');
n2=[n1(1:ind1),dec(ind2+1:end)];
lon_all3(ii)=str2num(n2);
end

n1=num2str(lat_all2(ii));
ind1=find(n1=='.');
if ~isempty(ind1)
dec=str2num(n1(ind1:end));
dec=dec/0.6;
dec=num2str(dec);
ind2=find(dec=='.');
n2=[n1(1:ind1),dec(ind2+1:end)];
lat_all3(ii)=str2num(n2);
end
end 
lat_all3=lat_all3';
lon_all3=lon_all3';

ind=find(isnan(lon_all));
lon_all(ind)=lon_all3;
lat_all(ind)=lat_all3;
clear lat_all2; lon_all2;  lat_all3; lon_all3;

% Strange (time) hour data 
ind=find(isnan(time_all));
text_data2=text_data(2:end,:);
text_data2=text_data2(ind,:);
text_data2=text_data2(:,9);
temp=str2num(cell2mat(text_data2));
temp=temp/2400;
time_all(ind)=temp;

% Add time to date
date_all=date_all+time_all;

% There are 13 points where volume filtered is not available 
ind=find(isnan(zoo_all));
zoo_all(ind)=[];
lat_all(ind)=[];lon_all(ind)=[];date_all(ind)=[];
time_all(ind)=[];bot_depth_all(ind)=[];vol_filt_all(ind)=[];
vol_disp_all(ind)=[];depth_lim_all(ind,:)=[];mesh_all(ind)=[];

% Missing time (hour) data
ind=find(isnan(time_all));
time_all(ind)=0.5; % Set to Noon for an average (not important).

% Missing depth max data
temp=mean(depth_lim_all,2);
ind=find(isnan(temp));
zoo_all(ind)=[];temp(ind)=[];
lat_all(ind)=[];lon_all(ind)=[];date_all(ind)=[];
time_all(ind)=[];bot_depth_all(ind)=[];vol_filt_all(ind)=[];
vol_disp_all(ind)=[];depth_lim_all(ind,:)=[];mesh_all(ind)=[];

% Check all data for nan
data2=[];
data2(1,:)=lon_all;data2(2,:)=lat_all;
data2(3,:)=date_all;data2(4,:)=time_all;
data2(5,:)=temp;data2(6,:)=bot_depth_all;
data2(7,:)=mesh_all;data2(8,:)=vol_filt_all;
data2(9,:)=vol_disp_all;data2(10,:)=zoo_all;
ind1=find(isnan(data2));
ind2=find(data2==-9);
[indx,indy]=find(data2==0);
% length(ind3) is equal to 4... 3 of the points are time (representing mid-night) 
% one of the points is volume 
unique(indx);

% A handful of points with lat, lon, or volume disp = 0;
ind_remove=[];
ind=find(lon_all==0);
ind_remove=cat(1,ind_remove,ind);
ind=find(lat_all==0);
ind_remove=cat(1,ind_remove,ind);
%ind=find(time_all==0); % Time is ok
ind=find(vol_filt_all==0);
ind_remove=cat(1,ind_remove,ind);
ind=ind_remove;
zoo_all(ind)=[];temp(ind)=[];
lat_all(ind)=[];lon_all(ind)=[];date_all(ind)=[];
time_all(ind)=[];bot_depth_all(ind)=[];vol_filt_all(ind)=[];
vol_disp_all(ind)=[];depth_lim_all(ind,:)=[];mesh_all(ind)=[];

% Check all data (AGAIN)
data2=[];
data2(1,:)=lon_all;data2(2,:)=lat_all;
data2(3,:)=date_all;data2(4,:)=time_all;
data2(5,:)=temp;data2(6,:)=bot_depth_all;
data2(7,:)=mesh_all;data2(8,:)=vol_filt_all;
data2(9,:)=vol_disp_all;data2(10,:)=zoo_all;
ind1=find(isnan(data2));
ind2=find(data2==-9);
ind3=find(data2==0);
% lenght(ind3)==7 all from time variable - this is ok

% Remove 202 samples (333/202 ratio computed in previous script, also this dataset is missing some 202 samples)
% mesh 12 = 202 sample
% mesh 3 = 333 sample
ind=find(mesh_all==12);
zoo_all(ind)=[];temp(ind)=[];
lat_all(ind)=[];lon_all(ind)=[];date_all(ind)=[];
time_all(ind)=[];bot_depth_all(ind)=[];vol_filt_all(ind)=[];
vol_disp_all(ind)=[];depth_lim_all(ind,:)=[];mesh_all(ind)=[];

% Simple Filter 
UB=mean(zoo_all)+2.5*std(zoo_all);
LB=mean(zoo_all)-2.5*std(zoo_all);
ind=find(zoo_all<LB | zoo_all>UB);
zoo_all(ind)=[];temp(ind)=[];
lat_all(ind)=[];lon_all(ind)=[];date_all(ind)=[];
time_all(ind)=[];bot_depth_all(ind)=[];vol_filt_all(ind)=[];
vol_disp_all(ind)=[];depth_lim_all(ind,:)=[];mesh_all(ind)=[];
% removed 212 outliers

% ******* Testing for Replicas *********
data2=[];
data2(:,1)=lon_all;data2(:,2)=lat_all;
data2(:,3)=date_all;data2(:,4)=time_all;
data2(:,5)=depth_lim_all(:,1); data2(:,6)=depth_lim_all(:,2);
data2(:,7)=bot_depth_all;
data2(:,8)=mesh_all;data2(:,9)=vol_filt_all;
data2(:,10)=vol_disp_all;data2(:,11)=zoo_all;
pos=unique(data2,'rows');
size(data2,1)-size(pos,1);
% 5 replicas? 

data2=[];
data2(:,1)=lon_all;data2(:,2)=lat_all;
data2(:,3)=date_all;data2(:,4)=time_all;
data2(:,5)=depth_lim_all(:,1); data2(:,6)=depth_lim_all(:,2);
data2(:,7)=mesh_all;data2(:,8)=vol_filt_all;
pos=unique(data2,'rows');
size(data2,1)-size(pos,1);
% 35 replicas 

data2=[];
data2(:,1)=lon_all;data2(:,2)=lat_all;
data2(:,3)=date_all;data2(:,4)=time_all;
data2(:,5)=depth_lim_all(:,1); data2(:,6)=depth_lim_all(:,2);
data2(:,7)=mesh_all;
pos=unique(data2,'rows');
size(data2,1)-size(pos,1);
% 199 replicas 

data2=[];
data2(:,1)=lon_all;data2(:,2)=lat_all;
data2(:,3)=date_all;data2(:,4)=time_all;
pos=unique(data2,'rows');
size(data2,1)-size(pos,1);
% 214 replicas

data2=[];
data2(:,1)=lon_all;data2(:,2)=lat_all;
pos=unique(data2,'rows');
size(data2,1)-size(pos,1);
% 1387 replicas
% ****************************** % 

% ################################ % 
% Need to be able to identify observations based on time and location
% so I need to fix the 169 issue
data2=[];
date_all_new=round(date_all);
lat_all_new=round(lat_all*1e5)/1e5; 	% Same precision that I read the model output coordinates 
lon_all_new=round(lon_all*1e5)/1e5;	% Same precision that I read the model output coordinates
data2(:,1)=lon_all_new;data2(:,2)=lat_all_new; data2(:,3)=date_all_new;
pos=unique(data2,'rows');
size(data2,1)-size(pos,1);
% 413 replicas

lon_all2=[]; lat_all2=[]; date_all2=[]; time_all2=[]; depth_lim_all2=[]; 
bot_depth_all2=[]; vol_filt_all2=[]; vol_disp_all2=[]; mesh_all2=[]; zoo_all2=[];
for ii=1:size(pos,1)
ind=find(lon_all_new==pos(ii,1) & lat_all_new==pos(ii,2) & date_all_new==pos(ii,3));

ind=ind(1);	% just grab the first one

lon_all2(ii,1)=lon_all_new(ind);
lat_all2(ii,1)=lat_all_new(ind);
time_all2(ii,1)=time_all(ind);
date_all2(ii,1)=date_all_new(ind);
bot_depth_all2(ii,1)=bot_depth_all(ind);
depth_lim_all2(ii,:)=depth_lim_all(ind,:);
vol_filt_all2(ii,1)=vol_filt_all(ind);
vol_disp_all2(ii,1)=vol_disp_all(ind);
mesh_all2(ii,1)=mesh_all(ind);
zoo_all2(ii,1)=zoo_all(ind);
end 
lon_all=lon_all2; clear lon_all2;
lat_all=lat_all2; clear lat_all2;
time_all=time_all2; clear time_all2;
date_all=date_all2; clear date_all2;
bot_depth_all=bot_depth_all2; clear bot_depth_all2;
depth_lim_all=depth_lim_all2; clear depth_lim_all2;
vol_filt_all=vol_filt_all2; clear vol_filt_all2;
vol_disp_all=vol_disp_all2; clear vol_disp_all2;
mesh_all=mesh_all2; clear mesh_all2;
zoo_all=zoo_all2; clear zoo_all2;
% ####################################

% ---------------------------- %
% Get month, day, year from date 
temp=datestr(date_all,'mm/dd/yyyy');
day_all=str2num(temp(:,4:5));
month_all=str2num(temp(:,1:2));
year_all=str2num(temp(:,end-3:end));

% **************************************** %
% ************ Conversion **************** % 
% **************************************** %

display('Test 110519')
pause 

% Convert zooplankton data from displacement volume (ml/m3) to carbon mass (mg C/m3)
zoo_all=(log10(zoo_all)+1.434)./0.820;  % log mg C/m3
zoo_all=10.^zoo_all;                    % mg C/m3
% Equation: log CM = (log DV + 1.434)/0.820
% Note: found equation in (Moriarty and O'Brien, 2013), originally published by (Wiebe, 1988)
% emailed Todd O'Brien (todd.obrien@noaa.gov), the equation is log base 10.
% log rules: b^y=x, log base b (x) =y 

% Convert from mg C/m3 to mmol N/ m3 (for model comparison)
zoo_all=zoo_all/1000;          % g C/m3
zoo_all=zoo_all/12.011;        % mol C/m3
zoo_all=zoo_all*1000;          % mmol C/m3
zoo_all=zoo_all*16/106;        % mmol N/m3

% Convert to 202 (so we can compare with (LZ+PZ) in the model
sf_zoo=0.50731;		% computed in a previous script 
zoo_all=zoo_all/sf_zoo;

% **************************************** %

% ----- Grid SEAMAP data ------- % 
res=0.5;
lon=[min(lon_all):res:max(lon_all)];
lat=[min(lat_all):res:max(lat_all)];
[xc,yc]=meshgrid(lon,lat);

xc2=[];yc2=[];zoo_grid=[];num_grid=[];year_grid=[];
for xx=1:size(xc,2)-1
for yy=1:size(xc,1)-1
ind1=find(lon_all>=xc(yy,xx) & lon_all<xc(yy,xx+1) & lat_all>=yc(yy,xx) & lat_all<=yc(yy+1,xx));
z_box=zoo_all(ind1);
b1=nanmean(z_box)-2.5*nanstd(z_box);
b2=nanmean(z_box)+2.5*nanstd(z_box);
ind2=find(z_box<=b1 | z_box>=b2);
z_box(ind2)=[];
ind1(ind2)=[];

if ~isempty(ind1)
zoo_grid(yy,xx)=nanmean(z_box);
num_grid(yy,xx)=length(ind1);
year_grid(yy,xx)=length(unique(year_all(ind1)));
else
zoo_grid(yy,xx)=NaN;
num_grid(yy,xx)=NaN;
year_grid(yy,xx)=NaN;
end

xc2(yy,xx)=(xc(yy,xx)+xc(yy,xx+1))/2;
yc2(yy,xx)=(yc(yy,xx)+yc(yy+1,xx))/2;
end
end
% --------------------------- %

display('Done Before Figures');pause 

% ######################################## %
% ############## Figures ################# % 
% ######################################## %

%display('Pause Before Figures');pause 

% Zooplankton Abundance
figure(1)
subplot('position',[0.1 0.3 0.8 0.4])
%subplot('position',[0.1 0.55 0.8 0.4])
plot(cc(:,1),cc(:,2),'k'); hold on;
pcolor(xc2,yc2,zoo_grid); shading flat;colorbar
%axis([-99 -79 24 33]); set(gca,'xticklabel',[]); set(gca,'yticklabel',[])
caxis([0 0.35]); box on;axis equal; axis([-98.1 -77.3 24 31.2]);
%title('Zooplankton Biomass - >200 um (mmol N/m3)  ')

figure(2)
%subplot('position',[0.1 0.1 0.8 0.4])
subplot('position',[0.1 0.3 0.8 0.4])
field=log10(zoo_grid); x=[-2:0.25:0];
plot(cc(:,1),cc(:,2),'k'); hold on; contourf(xc2,yc2,field,x);
%axis([-99 -79 24 33]);
caxis([-2 0]); x=colorbar;
set(x,'YTick',-2:0.5:0); set(gca,'Ytick',[24:2:32])
xlabel(['Longiutude (','\circ',')']); ylabel(['Latitude (','\circ',')'])
axis equal;axis([-98.1 -77.3 24 31.2]);box on
%title('Log_1_0 Zoo Biomass')

% ####################################################
% Sample Density
figure(2)
subplot('position',[0.1 0.3 0.8 0.4])
plot(cc(:,1),cc(:,2),'k'); hold on
pcolor(xc2,yc2,num_grid);shading flat;
%y=[0:10:100];
%contourf(xc2,yc2,num_grid,y);shading flat;
x=colorbar;
caxis([0 100])
%set(gca,'xticklabel',[]); set(gca,'yticklabel',[])
axis equal;axis([-98.1 -80 24 31.2]);box on;
%title('Total Number of Samples')
set(x,'position',[0.85 0.315 0.03 0.371])
set(gca,'ytick',[26,28,30])
set(gca,'xtick',[-96:2:-82])

% Sample Density Year
figure(3)
subplot('position',[0.1 0.3 0.8 0.4])
plot(cc(:,1),cc(:,2),'k'); hold on
pcolor(xc2,yc2,year_grid);shading flat;
x=colorbar;
caxis([0 35])
%set(gca,'xticklabel',[]); set(gca,'yticklabel',[])
axis equal;axis([-98.1 -80 24 31.2]);box on;
%title('Total Number of Samples')
set(x,'position',[0.85 0.315 0.03 0.371])
set(gca,'ytick',[26,28,30])
set(gca,'xtick',[-96:2:-82])
set(x,'ytick',[0:5:35])
% ##########################################

subplot('position',[0.1 0.1 0.8 0.4])
plot(cc(:,1),cc(:,2),'k'); hold on
pcolor(xc2,yc2,year_grid);shading flat;colorbar
caxis([0 35]); axis([-98.5 -79 24 33]); set(gca,'Ytick',[24:2:32])
xlabel(['Longiutude (','\circ',')']); ylabel(['Latitude (','\circ',')'])
axis equal;axis([-98.1 -80 24 31.2]);box on
title('Years Sampled')

% Locations of Samples
figure(3)
plot(cc(:,1),cc(:,2),'k'); hold on; axis([-99 -79 16 33])
scatter(lon_all,lat_all,15,month_all,'filled'); colorbar;caxis([1 12])
title('Sample Locations and Month')

% ----- Time of samples ----- %
y=[min(year_all):max(year_all)];
year_count=[];
for ii=1:length(y)
ind=find(year_all==y(ii));
year_count(ii)=length(ind);
end
m=[1:12];
month_count=[];
for ii=1:length(m)
ind=find(month_all==m(ii));
month_count(ii)=length(ind);
end

% Timing of Samples
figure(4)
subplot('position',[0.1 0.575 0.8 0.4])
bar(y,year_count,'facecolor',[0.75 0.75 0.75])
xlabel('Year'); ylabel('Number of Samples')
axis([1980 2020 0 600])
subplot('position',[0.1 0.075 0.8 0.4])
bar(m,month_count,'facecolor',[0.75 0.75 0.75])
xlabel('Month'); ylabel('Number of Samples')
axis([0 13 0 4000])

% Depth Sampled
figure(5)
subplot('position',[0.1 0.575 0.8 0.4])
hist(depth_lim_all(:,2),20)
axis([0 450 0 4000])
xlabel('Tow Depth'); ylabel('Number of Samples')
subplot('position',[0.1 0.075 0.8 0.4])
hist(bot_depth_all,20)
xlabel('Bottom Depth'); ylabel('Number of Samples')

display('Done');pause

% Check
temp=[lon_all,lat_all,date_all];
pos=unique(temp,'rows');
if size(pos,1)~=length(lon_all)
display('Error');break;
end

% ---- Save data ----- %
%lat_all=round(lat_all*1e5)/1e5;
%lon_all=round(lon_all*1e5)/1e5;

out='/Net/gleam/taylor/SEAMAP/';
SM_lat=lat_all;
SM_lon=lon_all;
SM_depth_lim=depth_lim_all;
SM_bot_depth=bot_depth_all;
SM_date=date_all;
SM_zoo=zoo_all;
SM_zoo_grid=zoo_grid;
SM_xc_grid=xc2;
SM_yc_grid=yc2;

whos SM_lat SM_lon SM_depth_lim SM_bot_Depth SM_date SM_zoo SM_zoo_grid SM_xc_grid SM_yc_grid
data2=[];
data2(:,1)=SM_lat;
data2(:,2)=SM_lon;
data2(:,3)=SM_date;
data2(:,4)=SM_depth_lim(:,1);
data2(:,5)=SM_depth_lim(:,2);
data2(:,6)=SM_bot_depth;
data2(:,7)=SM_zoo;
ind=find(isnan(data2));
pos1=unique(data2,'rows');
pos2=unique(data2(:,1:3),'rows');
% It is ok that pos2 has less because I am using rounded lat lon positions - when I use this
% data for sampling the model I will need to repeat the code written above where I find each 
% unique position and just select the first observation. In the code I will need to round the time 
% as well for preparation for model input.

% I went ahead an did that portion in this code - this is not highest accuracy SEAMAP data but 
% these rounded positions are needed for the model 

display('pause before save');pause 
save([out,'SEAMAP_data_072419.mat'],'SM_lat','SM_lon','SM_depth_lim','SM_bot_depth','SM_date','SM_zoo', 'SM_zoo_grid','SM_xc_grid','SM_yc_grid')
% ---------------------- % 







