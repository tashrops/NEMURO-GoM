% ****************************** % 
% Read Primary Production Data.
% Written by: Taylor Shropshire
% Date: 06/19/19
% ****************************** % 
clear;clc;close all

% Paths 
path='/home/taylor/my_scripts/projects/cruise_data/data_sets/primary_prod/';
file='PP_2017_2018_simplified_data.xlsx';

% Get coastline
load('/nexsan/people/taylor/projects/NPZ/data/matlab/eastcoastmask.mat') % Coastline
clear lat lon

% Read data file 
data_all=importdata([path,file]);
data=data_all.data;

% Define variables 
cycle=data(:,1);
day=data(:,2);
date=data(:,3);
int_depth=data(:,4);
int_pp=data(:,5);

start_date=datenum(2017,5,11);
date_scale=start_date-date(1);
date=date+date_scale;

NF_pp_date=date;
NF_int_pp=int_pp;
NF_pp_depth=int_depth;

path_out='/home/taylor/my_scripts/projects/cruise_data/data_sets/observation_positions/';
save([path_out,'NF_pp_data.mat'],'NF_pp_date','NF_int_pp','NF_pp_depth');



