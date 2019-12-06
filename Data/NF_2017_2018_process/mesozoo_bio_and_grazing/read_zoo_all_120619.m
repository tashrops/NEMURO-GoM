% ****************************** % 
% Read NF17 and NF18 mesozooplankton tow data
% Written by: Taylor Shropshire
% Date: 8/13/19
% ****************************** % 
clear;clc;close all

% Paths 
path1='/home/taylor/my_scripts/projects/cruise_data/data_sets/mesozoo_all/Mesozooplankton_FINAL_20190724.xlsx';
path2='/nexsan/people/taylor/projects/NPZ/data/matlab/eastcoastmask.mat';

% Load coastline 
load(path2,'cc');

% Read spreadsheet
data_all=importdata(path1);
data=data_all.data;
data_text=data_all.textdata;

% Flags
pyro_flag=2;
class_flag=1; % have not tested class_flag=2

% ---- Define Variables ---- % 
year=data(:,3);
month=data(:,4);
day=data(:,5);
date=datenum(year,month,day);
lat=data_text(2:end,9);
lon=data_text(2:end,10);
tow_num=data_text(2:end,12);
depth=data(:,17);
s_frac=data_text(2:end,21);
zoo_bio=data(:,31);			% Biomass mg C m-2
zoo_bio_n=data(:,32);			% Biomass mg N m-3 
zoo_graz=data(:,45);			% ug Chl m-2 d-1
zoo_spec_graz=zoo_graz./zoo_bio;	% ug Chl mg C-1 d-1
% ------------------------- % 

% ----- Convert Lat/Lon ----- % 
lat2=[]; lon2=[];
for ii=1:length(lat);
temp=cell2mat(lat(ii,1));
deg=str2num(temp(:,1:2));
dec_min=str2num(temp(:,4:end-2));
dec_deg=dec_min/60;
lat2(ii,1)=deg+dec_deg;

temp=cell2mat(lon(ii,1));
deg=str2num(temp(:,1:3));
dec_min=str2num(temp(:,4:end-2));
dec_deg=dec_min/60;
lon2(ii,1)=(deg+dec_deg)*-1;
end 
lat=lat2; clear lat2;
lon=lon2; clear lon2; 
% --------------------- % 

% ------- Get limits in usable form ------- %
s_frac2=[];
for ii=1:length(s_frac)
temp=cell2mat(s_frac(ii));

if length(temp)>=3
if (strcmp(temp(1:3),'0.2')==1)
s_frac2(ii,1)=str2num(temp(1:3));s_frac2(ii,2)=str2num(temp(5:end));
end 

if strcmp(temp(1:3),'0.5')==1
s_frac2(ii,1)=str2num(temp(1:3));s_frac2(ii,2)=str2num(temp(5));
end
end 

if strcmp(temp(1),'1')==1
s_frac2(ii,1)=str2num(temp(1));s_frac2(ii,2)=str2num(temp(3));
end

if strcmp(temp(1),'2')==1
s_frac2(ii,1)=str2num(temp(1));s_frac2(ii,2)=str2num(temp(3));
end

if strcmp(temp(1),'>')==1
s_frac2(ii,1)=str2num(temp(2));s_frac2(ii,2)=9999;
end

if strcmp(temp(1),'P')==1
s_frac2(ii,1)=NaN; s_frac2(ii,2)=NaN;
end

end 
s_frac=s_frac2;clear s_frac2;
% ------------------------------------ % 

% ----- Get tow number in usable form ----- %
tow_num2=[];
for ii=1:length(s_frac)
temp=cell2mat(tow_num(ii,1));

if length(temp)==3
tow_num2(ii,1)=str2num(temp(end));
end
if length(temp)==4
tow_num2(ii,1)=str2num(temp(end-1:end));
end

end
tow_num=tow_num2;clear tow_num2;
% ----------------------------- % 

% ######################################################## % 
% ##################### Process Tows ##################### % 
% ######################################################## % 

% ------ Filter Data ------ % 

% NaN values for grazing 
ind=find(isnan(zoo_graz));
date(ind)=[]; lat(ind)=[]; lon(ind)=[];
depth(ind)=[]; tow_num(ind)=[]; s_frac(ind,:)=[];
zoo_bio(ind)=[]; zoo_bio_n(ind)=[]; zoo_graz(ind)=[]; zoo_spec_graz(ind)=[];

% Remove Tows with Pyrosomes
%if pyro_flag==1
%ind=find(isnan(s_frac(:,1)));	% pyrosomes
%bad_tows=unique(tow_num(ind));

%for ii=1:length(bad_tows);
%ind=find(tow_num==bad_tows(ii));
%tow_num(ind)=-999;
%end
%ind=find(tow_num==-999);

%date(ind)=[]; lat(ind)=[]; lon(ind)=[];
%depth(ind)=[]; tow_num(ind)=[]; s_frac(ind,:)=[];
%zoo_bio(ind)=[]; zoo_bio_n(ind)=[]; zoo_graz(ind)=[]; zoo_spec_graz(ind)=[];
%end 

% Only remove Pyrosome Sample
if pyro_flag==2
ind=find(isnan(s_frac(:,1)));   % pyrosomes
date(ind)=[]; lat(ind)=[]; lon(ind)=[];
depth(ind)=[]; tow_num(ind)=[]; s_frac(ind,:)=[];
zoo_bio(ind)=[]; zoo_bio_n(ind)=[]; zoo_graz(ind)=[]; zoo_spec_graz(ind)=[];
end 

% Remove size class >5 mm
if class_flag==1
ind=find(s_frac(:,1)==5);
date(ind)=[]; lat(ind)=[]; lon(ind)=[];
depth(ind)=[]; tow_num(ind)=[]; s_frac(ind,:)=[];
zoo_bio(ind)=[]; zoo_bio_n(ind)=[]; zoo_graz(ind)=[]; zoo_spec_graz(ind)=[];
end 

% Only analyze tows with all 5 size classes
%if class_flag==2
%ind=find(s_frac(:,1)==5);
%uni_tows=unique(tow_num(ind));
%flag=tow_num*0;
%for ii=1:length(uni_tows)
%ind=find(tow_num==uni_tows);
%flag(ind)=1;
%end 
%ind=find(flag==0);
%date(ind)=[]; lat(ind)=[]; lon(ind)=[];
%depth(ind)=[]; tow_num(ind)=[]; s_frac(ind,:)=[];
%zoo_bio(ind)=[]; zoo_bio_n(ind)=[]; zoo_graz(ind)=[]; zoo_spec_graz(ind)=[];
%end 

% Remove tow 15 (unique size class 0.5 - 2)
ind=find(tow_num==15);
date(ind)=[]; lat(ind)=[]; lon(ind)=[];
depth(ind)=[]; tow_num(ind)=[]; s_frac(ind,:)=[];
zoo_bio(ind)=[]; zoo_bio_n(ind)=[]; zoo_graz(ind)=[]; zoo_spec_graz(ind)=[];
% ---------------------- % 

uni_tows=unique(tow_num);

tow_count=[];
for ii=1:length(uni_tows);
ind=find(tow_num==uni_tows(ii));
tow_count(ii,1)=length(ind);
end

ind=find(tow_count~=4);
bad_tow=uni_tows(ind);

for ii=1:length(bad_tow);
ind=find(tow_num==bad_tow(ii));
tow_num(ind)=-999;
end
ind=find(tow_num==-999);
date(ind)=[]; lat(ind)=[]; lon(ind)=[];
depth(ind)=[]; tow_num(ind)=[]; s_frac(ind,:)=[];
zoo_bio(ind)=[]; zoo_bio_n(ind)=[]; zoo_graz(ind)=[]; zoo_spec_graz(ind)=[];

uni_tows=unique(tow_num);
tow_count=[]; tow_year=[];
for ii=1:length(uni_tows);
ind=find(tow_num==uni_tows(ii));
tow_count(ii,1)=length(ind);
tow_year(ii,1)=unique(date(ind));
end

display('Done for now');pause 

% ---------- Process Meso Data ----- % 
SL_meso_spec_graz=[]; 
tot_meso_spec_graz=[];
SL_meso_bio=[]; 
zoo_per=[]; 

sfrac_tow=[]; date_tow=[]; lon_tow=[];
lat_tow=[]; depth_tow=[];

for ii=1:length(uni_tows);
ind=find(tow_num==uni_tows(ii));

temp1=zoo_bio(ind);		% used for calculating specific grazing for "LZ" and "PZ"
temp2=zoo_bio_n(ind);		% used for calculating percent biomass - because units of nitrogen are what the model is in
temp3=zoo_graz(ind);		% zooplankton grazing in units of ug Chl m-2 d-1 
temp4=s_frac(ind,:);

% ........... Grazing .................... % 
ind_sc1=find(temp4(:,1)==0.2 & temp4(:,2)==0.5);
ind_sc2=find(temp4(:,1)==0.5 & temp4(:,2)==1);
ind_sc3=find(temp4(:,1)==1 & temp4(:,2)==2);
ind_sc4=find(temp4(:,1)==2 & temp4(:,2)==5);
LZ_graz=[]; LZ_graz=[temp3(ind_sc1), temp3(ind_sc2)];
PZ_graz=[]; PZ_graz=[temp3(ind_sc3), temp3(ind_sc4)];
LZ_bio=[]; LZ_bio=[temp1(ind_sc1), temp1(ind_sc2)];
PZ_bio=[]; PZ_bio=[temp1(ind_sc3), temp1(ind_sc4)];

ind1=find(isnan(LZ_graz) | isnan(LZ_bio));
if isempty(ind1)
SL_meso_spec_graz(ii,1)=sum(LZ_graz)/sum(LZ_bio);
else
SL_meso_spec_graz(ii,1)=NaN;
end 

ind2=find(isnan(PZ_graz) | isnan(PZ_bio));
if isempty(ind2)
SL_meso_spec_graz(ii,2)=sum(PZ_graz)/sum(PZ_bio);
else
SL_meso_spec_graz(ii,2)=NaN;
end
% ............................................. % 

% Total Grazing
if isempty([ind1,ind2])
tot_meso_spec_graz(ii,1)=(sum(LZ_graz)+sum(PZ_graz))/(sum(LZ_bio)+sum(PZ_bio));
else 
tot_meso_spec_graz(ii,1)=NaN;
end 

% Percent Biomass
LZ_bio=[]; LZ_bio=[temp2(ind_sc1), temp2(ind_sc2)];
PZ_bio=[]; PZ_bio=[temp2(ind_sc3), temp2(ind_sc4)];
SL_meso_bio(ii,1)=sum(LZ_bio)/(sum(temp2))*100;
SL_meso_bio(ii,2)=sum(PZ_bio)/(sum(temp2))*100;

zoo_per=cat(1,zoo_per,temp2/sum(temp2)*100);
% ok that the order of size classes is different because I organize by size class below

sfrac_tow(ii,:,:)=temp4;
date_tow(ii,1)=mean(date(ind));
lon_tow(ii,1)=mean(lon(ind));
lat_tow(ii,1)=mean(lat(ind));
depth_tow(ii,1)=mean(depth(ind));
end 
% ----------------------------------- % 

whos SL_meso_spec_graz SL_meso_bio

% ********* Need to save for dilution stuff ******************
%meso_graz=sum(SL_meso_graz,2);
%out='/home/taylor/my_scripts/projects/cruise_data/data_sets/observation_positions/';
%whos meso_graz date_tow lon_tow lat_tow depth_tow
%save([out,'NF_meso_for_dilution.mat'],'meso_graz','date_tow','lon_tow','lat_tow','depth_tow');
% *************************************************** %

% --------- Organize based on size fraction ------------- % 
uni_s_frac=unique(s_frac,'rows');

frac_count=[]; zoo_bio_bin=[]; zoo_spec_graz_bin=[];
for ii=1:size(uni_s_frac,1);
ind=find(s_frac(:,1)==uni_s_frac(ii,1) & s_frac(:,2)==uni_s_frac(ii,2));
frac_count(ii,1)=length(ind);
zoo_bio_bin(:,ii)=zoo_bio_n(ind);
zoo_per_bin(:,ii)=zoo_per(ind);
zoo_spec_graz_bin(:,ii)=zoo_spec_graz(ind);
end 

% Check 
if length(unique(frac_count))>1
display('Error');break
end

% Average percent biomass in each size fraction
figure(1)
boxplot(zoo_per_bin,'labels',{'0.2-0.5';'0.5-1.0';'1.0-2.0';'2.0-5.0'},'color','k','whisker',0.7193,'symbol','+k')
axis([0 5 0 60])

figure(2)
boxplot(SL_meso_bio,'color','k','whisker',0.7193,'symbol','+k')
axis([0 3 0 100])

figure(3)
boxplot(zoo_spec_graz_bin,'labels',{'0.2-0.5';'0.5-1.0';'1.0-2.0';'2.0-5.0'},'color','k','whisker',0.7193,'symbol','+k')
axis([0 5 0 7])
ylabel('ug Chl mg C^-^1 d^-1')

figure(4)
boxplot(SL_meso_spec_graz,'color','k','whisker',0.7193,'symbol','+k')
axis([0 3 0 7])

display('Pause Before Save');pause;pause 
break

% ********* Save Data ********** % 
NF_meso_date=date_tow;
NF_meso_lat=lat_tow;
NF_meso_lon=lon_tow;
NF_meso_depth=depth_tow;
NF_meso_tow=uni_tows;
NF_meso_bio=zoo_bio_bin;
NF_meso_per=zoo_per_bin;
NF_meso_spec_graz_LZPZ=SL_meso_spec_graz;
NF_meso_per_bio_LZPZ=SL_meso_bio;
NF_meso_spec_graz=zoo_spec_graz_bin;
NF_meso_sfrac=uni_s_frac;

whos NF_meso_date NF_meso_lat NF_meso_lon NF_meso_depth NF_meso_tow NF_meso_bio NF_meso_per 
whos NF_meso_spec_graz_LZPZ NF_meso_spec_graz NF_meso_sfrac NF_meso_per_bio_LZPZ

out='/home/taylor/my_scripts/projects/cruise_data/data_sets/observation_positions/';
save([out,'NF_meso_data_092719.mat'],'NF_meso_date','NF_meso_lat','NF_meso_lon','NF_meso_depth','NF_meso_tow','NF_meso_sfrac',...
'NF_meso_bio','NF_meso_per','NF_meso_spec_graz','NF_meso_spec_graz_LZPZ','NF_meso_per_bio_LZPZ');











