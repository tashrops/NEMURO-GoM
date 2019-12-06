% ******************************* %
% Read SEAMAP CTD Flouresence Data 
% Data source: Glenn Zapfre (glenn.zapfe@noaa.gov)
% Written by: Taylor Shropshire
% Date: 06/24/18
% ******************************* %
clear;clc;close all

% Parameters 
create_mat_flag=0;
first_date=datenum(2008,02,07,23,21,44);

% ------- Read in and create Matlab File (done once) ------ % 
% Paths 
file1='/nexsan/people/taylor/projects/mit_gcm/zooplankton_data/CTD/WINTER_CTD_CASTS.xls';
file2='/nexsan/people/taylor/projects/mit_gcm/zooplankton_data/CTD/SPRING_CTD_CASTS.xls';
file3='/nexsan/people/taylor/projects/mit_gcm/zooplankton_data/CTD/FALL_CTD_CASTS.xls';
file4='/nexsan/people/taylor/projects/mit_gcm/zooplankton_data/CTD/SPRING2009_CTD_CASTS.xls';
file5='/nexsan/people/taylor/projects/mit_gcm/zooplankton_data/CTD/SPRING2011_CTD_CASTS.xls';
out='/nexsan/people/taylor/projects/mit_gcm/zooplankton_data/';

coast_file='/nexsan/people/taylor/projects/NPZ/data/matlab/eastcoastmask.mat';
load(coast_file);clear lat lon;

file4
tic
data_all=importdata(file4);		% < ------------- just change file1 to file1-5
data=data_all.data;

%data=cell2mat(struct2cell(data));

sheet_name=fieldnames(data);
data2=[];
for ii=1:length(sheet_name)
name=cell2mat(sheet_name(ii));
sheet=getfield(data,name);
%sheet=sheet(:,1:9);
if size(sheet,2)<11
sheet(:,end:11)=Inf;	% just a place holder
end
data2=cat(1,data2,sheet);
end 
data=data2;clear data2;

% --------- Inital quality control ------- % 
% Remove headers 
ind=find(isnan(data(:,1)));
data(ind,:)=[];

% Remove text columns
ind=[3,10];
data(:,ind)=[];

% Remove other data not needed 
ind=[2,8,9];
data(:,ind)=[];

% Bad data
[row,col]=find(isnan(data));
row=unique(row);
data(row,:)=[];

% Check
ind=find(isnan(data));
if ~isempty(ind);
display('Still NaN Values'); pause 
end

% Some zeros exist for fluoresence
[row,col]=find(data(:,6)==0);
row=unique(row);
data(row,:)=[];
% --------------------------- %

% Remove bad depths 
[row,col]=find(data(:,5)>500);
row=unique(row);

if ~isempty(row)	% 44 bad casts exist for spring file
b_lat=data(row,2);
b_lon=data(row,3);
b_time=data(row,1);
temp=[b_time,b_lat,b_lon];
pos_temp=unique(temp,'rows');
for ii=1:size(pos_temp,1);
ind=find(data(:,1)==pos_temp(ii,1) & data(:,2)==pos_temp(ii,2) & data(:,3)==pos_temp(ii,3));
if isempty(ind)
display('Error');pause
end
data(ind,:)=NaN;
end 
[row,col]=find(isnan(data));
row=unique(row);
data(row,:)=[];
end 

% Define variables 
time_all=data(:,1);
%time_scale=first_date-time_all(1);
time_scale=693960;
time_all=time_all+time_scale;

lat_all=data(:,2);
lon_all=data(:,3);
bot_depth_all=data(:,4);
depth_all=data(:,5);
flo_all=data(:,6);

display(['Samples from = ', datestr(min(time_all)), ' ----> ', datestr(max(time_all))]);

% Find unique cast postions
temp=[time_all,lat_all,lon_all, bot_depth_all];
pos=unique(temp,'rows');
time=pos(:,1);
lat=pos(:,2);
lon=pos(:,3);
bot_depth=pos(:,4);

figure(1)
plot(cc(:,1),cc(:,2),'k');axis([-98 -78 17 32]); hold on
plot(lon,lat,'.r')

% Find max number of measurements in a cast
num_obs=[];
for ii=1:length(time)
ind=find(time_all==time(ii) & lat_all==lat(ii) & lon_all==lon(ii));
temp=flo_all(ind);
num_obs(ii)=length(ind);
end

% Check 
if min(num_obs)==0
display('Profile Error');pause 
end 

% Organize profiles
flo_prof=zeros(max(num_obs),length(time))*NaN; depth_prof=flo_prof;
count_neg=0; count_prof=1; count_bad=1; bad_prof=[];
for ii=1:length(time)
ind=find(time_all==time(ii) & lat_all==lat(ii) & lon_all==lon(ii));
flo=flo_all(ind);

% Scale flo reading to be positive 
ind_neg=find(flo<0);
if ~isempty(ind_neg)
flo=flo+abs(min(flo))+1e-10;
count_neg=count_neg+1;
end 

if length(unique(flo))>1
flo_prof(1:length(ind),count_prof)=flo;
depth_prof(1:length(ind),count_prof)=depth_all(ind);
count_prof=count_prof+1;
else 
bad_prof(count_bad)=ii;
count_bad=count_bad+1;
end 

end

% Remove profiles with all the same values
flo_prof(:,bad_prof)=[];
depth_prof(:,bad_prof)=[];
time(bad_prof)=[];
lat(bad_prof)=[];
lon(bad_prof)=[];
num_obs(bad_prof)=[];
bot_depth(bad_prof)=[];

% Remove casts that started deeper
ind=find(depth_prof(1,:)>30);	% only a few 
flo_prof(:,ind)=[];
depth_prof(:,ind)=[];
time(ind)=[];
lat(ind)=[];
lon(ind)=[];
num_obs(ind)=[];
bot_depth(ind)=[];

% Remove cast where NaN is the first value
ind=find(isnan(flo_prof(1,:))); 
% this occurs because there is some skipped filling since I am using 
% a count variable and not the loop variable
flo_prof(:,ind)=[];
depth_prof(:,ind)=[];
time(ind)=[];
lat(ind)=[];
lon(ind)=[];
num_obs(ind)=[];
bot_depth(ind)=[];

% Remove casts where the max depth is less than 20m
depth_lim=[];
for ii=1:length(time)
temp=depth_prof(:,ii);
depth_lim(ii,1)=min(depth_prof(:,ii));
depth_lim(ii,2)=max(depth_prof(:,ii));
end 

ind=find(depth_lim(:,2)<30);
flo_prof(:,ind)=[];
depth_prof(:,ind)=[];
time(ind)=[];
lat(ind)=[];
lon(ind)=[];
num_obs(ind)=[];
bot_depth(ind)=[];

[row,col]=find(depth_prof>500);
col=unique(col);
if ~isempty(col) 
display('Error');pause
end 

display('Done');pause
% *************************************************** % 
% **************** Done ***************************** % 
% *************************************************** % 

% Plot profiles 
figure(2); set(gca,'ydir','reverse');
cast_std=[];
for ii=1:length(time)
cast_std(ii)=nanstd(flo_prof(:,ii));
f_prof=flo_prof(:,ii);
d_prof=depth_prof(:,ii);
%f_prof=f_prof/f_prof(1);

plot(f_prof,d_prof,'k');
set(gca,'ydir','reverse');axis([0 25 0 200])
title(['Date = ', datestr(time(ii)), ', ', num2str(ii)])
pause(0.25)
end 

DCM_depth=[];
for ii=1:length(time)
flo=flo_prof(:,ii);
depth=depth_prof(:,ii);
ind=find(flo==max(flo));
DCM_depth(ii)=depth(ind(1));
end 

d_test=[0:10:200];
DCM_count=[];d_bin=[];
for zz=1:length(d_test)-1
ind=find(DCM_depth>=d_test(zz) & DCM_depth<d_test(zz+1));
DCM_count(zz)=length(ind);
d_bin(zz)=(d_test(zz)+d_test(zz+1))/2;
end 

figure(3)
barh(d_bin,DCM_count,'facecolor',[0.5 0.5 0.5])
set(gca,'ydir','reverse')
ylabel('Depth');
xlabel('Frequency')
title('Fall')
