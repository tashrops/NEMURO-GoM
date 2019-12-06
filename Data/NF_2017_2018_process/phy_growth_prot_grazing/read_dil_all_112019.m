clear;clc;close all

path1='/home/taylor/my_scripts/projects/cruise_data/data_sets/microzoo_grazing/NF17NF18_dil_112019.xlsx';
path2='/home/taylor/my_scripts/projects/cruise_data/data_sets/microzoo_grazing/Drifter_Tracks_NF17.xlsx';
path3='/home/taylor/my_scripts/projects/cruise_data/data_sets/microzoo_grazing/Drifter_Tracks_NF18.xlsx';
path4='/nexsan/people/taylor/projects/NPZ/data/matlab/eastcoastmask.mat';

% Load Coastline
load(path4,'cc');

% ------------------------------ % 
% Load dilution data
data_all=importdata([path1]);
data=data_all.data;
data_text=data_all.textdata;

p_spec=data(:,12);
z_spec=data(:,13);
date=data(:,15);
depth=data(:,16);

ind=find(isnan(p_spec));
p_spec(ind)=[];
z_spec(ind)=[];
date(ind)=[];
depth(ind)=[];

date_scale=datenum(2017,5,11)-date(1);
dil_date=date+date_scale;
% -------------------------------% 

% ------------ Load data 2017 drifter data ------------------ % 
data_all=importdata([path2]);
data=data_all.data;
data_text=data_all.textdata;

% Define Time
date=data(:,4:5);
date_start=datenum(2017,5,11,9,20,0);
date_scale=date_start-date(1,1);
date1=date+date_scale; clear date;
date2=datestr(date1(:,1),31);
y=str2num(date2(:,1:4));
m=str2num(date2(:,6:7));
d=str2num(date2(:,9:10));

date_17=datenum(y,m,d);
lat_17=data(:,8);
lon_17=data(:,9);
% ----------------------------------------------------------- %

% ------------ Load data 2018 drifter data ------------------ % 
data_all=importdata([path3]);
data=data_all.data;
data_text=data_all.textdata;

% Define Time
date=data(:,4:5);
date_start=datenum(2018,5,5,17,50,0);
date_scale=date_start-date(1,1);
date1=date+date_scale; clear date;
date2=datestr(date1(:,1),31);
y=str2num(date2(:,1:4));
m=str2num(date2(:,6:7));
d=str2num(date2(:,9:10));

date_18=datenum(y,m,d);
lat_18=data(:,8);
lon_18=data(:,9);
% ----------------------------------------------------------- %

lat=[lat_17; lat_18];
lon=[lon_17; lon_18];
date=[date_17;date_18];

% Find coordinates that coorespond to sample experiment date
lat_expt=[]; lon_expt=[];
for ii=1:length(dil_date)
ind=find(date==dil_date(ii));
lat_expt(ii,1)=mean(lat(ind));
lon_expt(ii,1)=mean(lon(ind));
end

figure(1)
plot(cc(:,1),cc(:,2),'k');hold on
plot(lon_expt,lat_expt,'.r');
axis([-98 -78 17 32])

NF_dil_phy_sg_int=p_spec;
NF_dil_zoo_sg_int=z_spec;
NF_dil_lat=lat_expt;
NF_dil_lon=lon_expt;
NF_dil_date=dil_date;
NF_dil_depth=depth;

whos NF_dil_phy_sg_int NF_dil_zoo_sg_int NF_dil_lat NF_dil_lon NF_dil_date NF_dil_depth

% ----save data ----- %
out='/home/taylor/my_scripts/projects/cruise_data/data_sets/observation_positions/';
save([out,'NF_dil_data_112019.mat'],'NF_dil_depth','NF_dil_date','NF_dil_lon',...
'NF_dil_lat','NF_dil_zoo_sg_int','NF_dil_phy_sg_int'])
% ---------------------- % 







