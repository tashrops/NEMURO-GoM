% *************************************************************************** % 
% North Pacific Ecosystem Model for Understanding Regional Oceanography (NEMURO)
% Description: 1-D Lower Trophic Level Marine Biogeochemical Numerical Model
% Source: (Kishi et al 2007) 
% Written By: Taylor Shropshire 
% Features/Flags: Mixing, Sinking, Daily and Hourly light Forcing, 
% Light limited microzooplankton grazing, linear and quadratic mortality, Chl:C submodel  
% Major Changes from original code: Light limitation functional form, ammomnium inhibition
%
% -------- Example Code  ------- %
% Differential Equation:  dP/dt = -m*P.^2                            (Phytoplankton Quadratic Mortality)
%
% Forward Euler:
%       for t = 1:n
%       P(t+1) = P(t) - dt*(m*P(t).^2);                 (Not positive definite)
%       end 
%
% Backward Euler:
%       for t=1:n
%       P(t+1) = P(t);
%       for iter = 1:10
%       P(t+1) = P(t)/(1+dt*(m*P(t+1)));                 (Positive definite)
%       end 
%       end 
% Note: Forward and Backward Euler are approximately the same at small dt. 
%
% The equations below are written in backward Euler in the form of:
%       for t=1:n
%       P_v1 = P(t); 
%       P_v2 = P_v1;
%       for iter = 1:10
%       P_v2 = P_v1/(1+dt/P_v2*(m*P_v2.^2));              
%       end 
%       P(t+1)=P_v2;
%       end 
% ------------------------- %
% ************************************************************************  %
clear;clc;close all

% Configuration
time_step=3600;				% time step (s)
sim_dur=365*2;				% simulation duration (days)
kz=1.75e-4;                  		% vertical diffusivity (m^2/s)
sink=15;                    		% sinking speed per day - [40] 
sink_flag=1;				% sinking: 1=on, 0=off
mix_flag=1;                 		% mixing: 1=on, 0=off (mixing on slightly disrupts conservation)
season_flag=1;				% seasonal forcing: 1=on, 0=off
dayavg_flag=1;				% daily average forcing: 1=on, 0=off

% Other
dt=time_step/86400;                     % time step used in code (1/day)
kz_day=kz*86400;                        % vertical diffusivity scaled 
tot_iter=sim_dur*86400/time_step;       % total iterations 
tsteps_inhour=3600/time_step;		% total iterations in an hour 
tsteps_inday=86400/time_step;		% total iterations in a day 
tsteps_inyear=tsteps_inday*365;		% total iterations in a year 
start_diag1=tsteps_inyear;          	% start day for hourly diagnostic 
end_diag1=tsteps_inyear*2;          	% end day for hourly diagnostic
start_diag2=365;                        % start day for daily diagnostic
end_diag2=365*2;                        % end day for daily diagnost
dx=4000;                                % does not influence code (m) - conceptual 4km box
dy=4000;                                % does not influence code (m) - conceptual 4km box 


% Initalize
load('prof_info.mat')			% NO, SI, Light, Temperature, Depth (Nutrient Profiles - World Ocean Atlas, Tmp/Rad - HYCOM forcing file)
NO=NO_ini_prof;                         % nitrate
SI=SI_ini_prof;                         % silica
NH=zeros(length(depth_t),1)+0.01;      	% ammonium (mmol N/m3)
PON=NH;                                 % particulate organic nitrogen (mmol N/m3)
DON=NH;                                 % dissolved organic nitrogen (mmol N/m3)
OP=NH;                                  % opal (mmol SI/m3)
sp=NH;                                  % small phytoplankton (mmol N/m3)
lp=NH;                                  % large phytoplankton (mmol N/m3)
sz=NH;                                  % small zooplankton (mmol N/m3)
lz=NH;                                  % large zooplankton (mmol N/m3)
pz=NH;                                  % predatory zooplankton (mmol N/m3)
dz=diff(depth_int);
dz_w=diff(depth_int);
dz_t=diff(depth_t);
vol=dx*dy*dz; 

% Day Averaged Light and Temperature Forcing Flag
if dayavg_flag==1
display('Daily Average Forcing')
display(' ')
light_ts2=[];
tmp_ts2=[];
ind=[0:24:length(light_ts)];
for ii=1:length(ind)-1
light_ts2(ii)=mean(light_ts(1,ind(ii)+1:ind(ii+1)));
tmp_ts2(:,ii)=mean(tmp_ts(:,ind(ii)+1:ind(ii+1)),2);
end 
light_ts=light_ts2;clear light_ts2;
tmp_ts=tmp_ts2;clear tmp_ts2;
else
display('Hourly Forcing')
display(' ')
end

% ********** Parameters ********** %

% Light Atteuation
par_frac=0.43;                  % fraction of photosynthetically active radiation
ext_w=0.03;                     % attenuation due to water - [0.04]
ext_sp = 0.03;                  % atteunation due to small phytoplankton - [0.04]
ext_lp =0.03;                   % atteunation due to large phytoplankton - [0.04] 
rad_constant=900;               % constant irradience (used if seasonality flag = 0)

% Temperature Dependence
tdep1=0.0693; scale_tmp_phyto1=1;        % phytoplankton growth
tdep2=0.0693; scale_tmp_phyto2=1;        % phytoplankton mortality
tdep3=0.0519; scale_tmp_phyto3=1;        % phytolpankton respiration
tdep4=0.0693; scale_tmp_zoo1=1;          % zooplankton grazing (1:0.0693, 0.14:0.1386, 0.016:0.2079)
tdep5=0.0693; scale_tmp_zoo2=1;          % zooplankton mortality (1:0.0693, 0.14:0.1386)
tdep6=0.0693; scale_tmp_decomp=1;        % decomposition (1:0.0693, 0.14:0.1386)

% Small Phytoplankton 
vmax_sp = 0.4;                  % max growth rate at 0 deg C (1/d) - [0.4]
k_NO_sp = 0.5;                  % nitrate half saturation constant - [1.0]
k_NH_sp = 0.1;                  % ammonium half saturation constant - [0.1]
alpha_sp = 0.1;    	    	% PI curve parameter alpha - [0.01]
beta_sp = 4.5e-4;                 % PI curve parmater beta - [4.5e-4], 1.4e-3
ref_resp_sp = 0.03;            % respiration at 0 deg C (1/d) - [0.03]
ref_mort_sp = 0.002;              % mortality at at 0 deg C (1/d) - [0.0585], Linear Mort - ~[0.002]
ext_excr_sp = 0.135;             % extracellular excretion - [0.135]
tc_v_sp = tdep1;                % tmp coefficent for growth - [0.0693]
tc_r_sp = tdep3;                % tmp coefficent for respiration - [0.0519]
tc_m_sp = tdep2;                % tmp coeffiect for mortality - [0.0693]
inh_NH_NO_sp = 0;               % nitrate uptake inhibiton, NOTE: parameter currently not being used - [1.4] 
chl2c_sp_constant=0.01;         % chlorophyll : carbon ratio
sp_mort_type = 1;               % mortality type

% Large Phytoplankton 
vmax_lp = 0.8;                  % max growth rate at 0 deg C (1/d) - [0.8]
k_NO_lp = 3.0;                  % nitrate half saturation constant - [3.0]
k_NH_lp = 0.3;                  % ammonium half saturation constant - [0.3]
k_SI_lp = 6.0;                  % silica half saturation constant - [6.0] 
alpha_lp = 0.1;                 % PI curve parameter alpha - [0.01], 1.4e-3
beta_lp= 4.5e-4;                  % PI curve parmater beta - [4.5e-4]
ref_resp_lp = 0.03;             % respiration at 0 deg C (1/d) - [0.03]
ref_mort_lp = 0.001;              % mortality at at 0 deg C (1/d) - [0.029], Linear Mort - ~[0.001]
ext_excr_lp = 0.135;              % extracellular excretion - [0.0135]
tc_v_lp = tdep1;                % tmp coefficent for growth - [0.0693]
tc_r_lp = tdep3;                % tmp coefficent for respiration - [0.0519]
tc_m_lp = tdep2;                % tmp coeffiect for mortality - [0.0693]
inh_NH_NO_lp = 0;               % nitrate uptake inhibiton, NOTE: parameter currently not being used - [1.4]
chl2c_lp_constant=0.03;         % chlorophyll : carbon ratio
lp_mort_type = 1;               % mortality type 

% Small Zooplankton 
gmax_sz_sp = 0.6;               % max grazing rate on small phytoplankton at 0 deg C (1/d) - [0.04]
gmax2_sz_sp = 0.0;              % minumum grazing under no light.
rad_graz_sz_sp = 0.0;           % half saturation constant for light limited grazing     
tc_g_sz_sp = tdep4;             % tmp coefficent for grazing - [0.0693]
iv_sz_sp = 1.4;                 % ivlev constant - [1.4]
thresh_sz_sp = 0.04 ;           % feeding threshold on small phytoplankton - [0.043]
ref_mort_sz = 0.022;              % mortality at 0 deg C (1/d) - [0.0585], Linear Mort - ~[0.002]
tc_m_sz = tdep5;                % tmp coefficent for mortality - [0.0693]
ae_sz = 0.70;                   % assimilation efficency - [0.7]
gge_sz = 0.30;                  % gross growth efficency - [0.3]
sz_mort_type = 1;               % mortality type 

% Large Zooplankton 
gmax_lz_sp = 0.0;               % max grazing rate on small phytoplankton at 0 deg C (1/d) - [0.1]
gmax_lz_lp = 0.3;               % max grazing rate on large phytoplankton at 0 deg C (1/d) - [0.4]
gmax_lz_sz = 0.3;               % max grazing rate on small zooplankton at 0 deg C (1/d) - [0.4]
tc_g_lz_sp = tdep4;             % tmp coefficent for grazing -[0.0693]
tc_g_lz_lp = tdep4;             % tmp coefficent for grazing -[0.0693]
tc_g_lz_sz = tdep4;             % tmp coefficent for grazing -[0.0693]
iv_lz_sp = 1.4;                 % ivlev constant - [1.4]
iv_lz_lp = 1.4;                 % ivlev constant - [1.4]
iv_lz_sz = 1.4;                 % ivlev constant - [1.4]
thresh_lz_sp = 0.04;            % feeding threshold on small phytoplankton - [0.04]
thresh_lz_lp = 0.04;            % feeding threshold on large phytoplankton - [0.04]
thresh_lz_sz = 0.04;            % feeding threshold on small zooplankton  - [0.04]
ref_mort_lz = 0.022;              % mortality at at 0 deg C (1/d) - [0.0585], Linear Mort - ~[0.002]
tc_m_lz = tdep5;                % tmp coefficent for mortality - [0.0693]
ae_lz = 0.70;                   % assimilation efficency - [0.7]
gge_lz = 0.30;                  % gross growth efficency - [0.3]
lz_mort_type = 1;               % mortality type

% Predatory Zooplankton 
gmax_pz_lp = 0.1;               % max grazing rate on large phytoplankton at 0 deg C (1/d) - [0.2]
gmax_pz_sz = 0.1;               % max grazing rate on small zooplankton at 0 deg C (1/d) - [0.2]
gmax_pz_lz = 0.3;               % max grazing rate on large zooplankton at 0 deg C (1/d) - [0.2]
tc_g_pz_lp = tdep4;             % tmp coefficent for grazing -[0.0693]
tc_g_pz_sz = tdep4;             % tmp coefficent for grazing -[0.0693]
tc_g_pz_lz = tdep4;             % tmp coefficent for grazing -[0.0693]
iv_pz_lp = 1.4;                 % ivlev constant - [1.4]
iv_pz_sz = 1.4;                 % ivlev constant - [1.4]
iv_pz_lz = 1.4;                 % ivlev constant - [1.4]
thresh_pz_lp = 0.04;            % feeding threshold on large phytoplankton - [0.04]
thresh_pz_sz = 0.04;            % feeding threshold on small zooplankton - [0.04]
thresh_pz_lz = 0.04;            % feeding threshold on large zooplankton - [0.04]
inh_szlz_lp = 4.605;            % grazing on large phytoplankton inhibition by small and large zooplankton - [4.605]
inh_lz_sz = 3.01;               % grazing on small zooplankton inhibition by large zooplankton - [3.01]
ref_mort_pz = 0.12;              % mortality at at 0 deg C (1/d) - [0.0585], Linear Mort - ~[0.002]
tc_m_pz = tdep5;                % tmp coefficent for mortality - [0.0693]
ae_pz = 0.70;                    % assimilation efficency - [0.7]
gge_pz = 0.30;                  % gross growth efficency - [0.3] 
pz_mort_type = 2;               % mortality type 

% Nutrients (11)
ref_nitr = 0.003;                % nitrification at 0 deg C (1/d) - [0.03]
ref_dec_PON_NH = 0.01;           % decompositon from PON to NH at 0 deg C (1/d) - [0.1]
ref_dec_PON_DON = 0.05;          % decompositon from PON to DON at 0 deg C (1/d) - [0.1]
ref_dec_DON_NH = 0.02;           % decompositon from DON to NH at 0 deg C (1/d) - [0.2] 
ref_dec_OP_SI = 0.01;            % decompositon from OP to SI at 0 deg C (1/d) - [0.1]
tc_nitr = tdep6;                % tmp coefficent for nitrification - [0.0693] 
tc_dec_PON_NH = tdep6;          % tmp coefficent for decomposition  - [0.0693]
tc_dec_PON_DON = tdep6;         % tmp coefficent for decomposition - [0.0693]
tc_dec_DON_NH = tdep6;          % tmp coefficent for decomposition - [0.0693]
tc_dec_OP_SI = tdep6;           % tmp coefficent for decomposition - [0.0693]
r_SI_N = 1.0;                   % ratio of silica to nitrogen (SI/N) - [2.0]

min_val=1.0e-6;                 % min value 

% Chl2c Submodel (Li et al., 2010)
chl2c_sp_min = 0.0001;          % minimum Chl:C ratio - [0.000]
chl2c_lp_min = 0.005;           % maximum Chl:C ratio - [0.005] 
chl2c_sp_max = 0.015;            % minimum Chl:C ratio - [0.03]
chl2c_lp_max = 0.03;           % maximum Chl:C ratio - [0.061] 
alpha_chl_sp = 0.28;            % chlorophyll specific initial slope of PE curve - [0.28 +- 0.11]
alpha_chl_lp = 0.28;            % chlorophyll specific initial slope of PE curve - [0.28 +- 0.11]

% Other
path3='/nexsan/people/taylor/projects/mit_gcm/simulation_projects/5_grid_point_GoM/GoM5_expt001/details/figures/1D';
print_flag=0;

% Time Series Variables
NO_all=[];NH_all=[];SI_all=[];DON_all=[];PON_all=[];
OP_all=[];sp_all=[];lp_all=[];sz_all=[];lz_all=[];pz_all=[];
NO_all(:,1)=NO; NH_all(:,1)=NH;
PON_all(:,1)=PON; DON_all(:,1)=DON;
sp_all(:,1)=sp; lp_all(:,1)=lp;
sz_all(:,1)=sz; lz_all(:,1)=lz;
pz_all(:,1)=pz; SI_all(:,1)=SI;
OP_all(:,1)=OP;

% Calculate total Nitrogen and Silica in the System
tot_n(1)=sum(NO.*vol)+sum(NH.*vol)+sum(DON.*vol)+sum(PON.*vol)...
+sum(sp.*vol)+sum(sz.*vol)+sum(lz.*vol)+sum(pz.*vol);
tot_s(1)=sum(SI.*vol)+sum(OP.*vol)+sum(lp.*vol.*r_SI_N);

% Print Simulation Information:
display('Finished Initalization and Defining Parameter Values');
display(' ');
display('Begin 1-D NEMURO Simulation');
display(' ');
display(['Number of Days to be Run = ', num2str(sim_dur), ' Days'])
display(' ');
display(['Time Step = ', num2str(time_step), ' Seconds'])
display(' ');
if sp_mort_type==1 & lp_mort_type==1
display('Linear Mortiality')
display(' ')
else if sp_mort_type==2 & lp_mort_type==2
display('Quadtratic Mortality')
display(' ')
end
end
display('############# Simulation Start ################### ');

% ################################### % 
% ############# NEMURO ############## %
% ################################### % 
spec_g_sp_all=[];spec_g_lp_all=[];chl2c_sp=[];chl2c_lp=[];
gpp_all=[];phyto_resp_all=[];phyto_mort_all=[];phyto_graz_all=[];
light_lim_sp_all=[];light_lim_lp_all=[];nit_lim_sp_all=[];amm_lim_sp_all=[];
nit_lim_lp_all=[];amm_lim_lp_all=[];sil_lim_lp_all=[];nut_lim_lp_all=[];
spec_g_sz_all=[];spec_g_lz_all=[];spec_g_pz_all=[];sz_graz_lim_all=[];
lz_graz_lim_all=[];pz_graz_lim_all=[];phy_excr_all=[];sz_graz_all=[];
lz_graz_all=[];pz_graz_all=[];tmp_lim_v_sp_all=[];tmp_lim_v_lp_all=[];
graz_on_sp_all=[];graz_on_lp_all=[];ns_lim_comp_all=[];

gpp_all=zeros(length(depth_t),tot_iter)*NaN;
spec_g_sp_all=gpp_all;spec_g_lp_all=gpp_all;chl2c_sp=gpp_all;
chl2c_lp=gpp_all;phyto_resp_all=gpp_all;phyto_mort_all=gpp_all;
phyto_graz_all=gpp_all;light_lim_sp_all=gpp_all;light_lim_lp_all=gpp_all;
nit_lim_sp_all=gpp_all;amm_lim_sp_all=gpp_all;nit_lim_lp_all=gpp_all;
amm_lim_lp_all=gpp_all;sil_lim_lp_all=gpp_all;nut_lim_lp_all=gpp_all;
spec_g_sz_all=gpp_all;spec_g_lz_all=gpp_all;spec_g_pz_all=gpp_all;
sz_graz_lim_all=gpp_all;lz_graz_lim_all=gpp_all;pz_graz_lim_all=gpp_all;
phy_excr_all=gpp_all;sz_graz_all=gpp_all;lz_graz_all=gpp_all;
pz_graz_all=gpp_all;tmp_lim_v_sp_all=gpp_all;tmp_lim_v_lp_all=gpp_all;
graz_on_sp_all=gpp_all;graz_on_lp_all=gpp_all;

h_count=1;
d_count=1;
tic

% >>>> Main Loop <<<<< %
for tt=1:tot_iter
for zz=1:length(depth_t)

NO_v1=NO_all(zz,tt); NO_v2=NO_v1;
NH_v1=NH_all(zz,tt); NH_v2=NH_v1;
SI_v1=SI_all(zz,tt); SI_v2=SI_v1;
DON_v1=DON_all(zz,tt); DON_v2=DON_v1;
PON_v1=PON_all(zz,tt); PON_v2=PON_v1;
OP_v1=OP_all(zz,tt); OP_v2=OP_v1; 
sp_v1=sp_all(zz,tt); sp_v2=sp_v1;
lp_v1=lp_all(zz,tt); lp_v2=lp_v1;
sz_v1=sz_all(zz,tt); sz_v2=sz_v1;
lz_v1=lz_all(zz,tt); lz_v2=lz_v1;
pz_v1=pz_all(zz,tt); pz_v2=pz_v1;

% ---- Determine Temperature and Light ---- % 
% Check for seasonality 
if season_flag==1

% Check for daily average
if dayavg_flag==1
count_type=d_count;
else 
count_type=h_count;
end
% Compute Light and tmp availability 
if zz==1
rad=light_ts(count_type)*par_frac*exp(-1*(ext_sp*sp_v1+ext_lp*lp_v1+ext_w)*depth_t(zz));
else 
rad=rad*exp(-1*(ext_sp*sp_v1+ext_lp*lp_v1+ext_w)*(depth_t(zz)-depth_t(zz-1)));
end 
tmp=tmp_ts(zz,count_type);
% If Constant Light and tmp
else
if zz==1
rad=rad_constant*par_frac*exp(-1*(ext_sp*sp_v1+ext_lp*lp_v1+ext_w)*depth_t(zz));
else
rad=rad*exp(-1*(ext_sp*sp_v1+ext_lp*lp_v1+ext_w)*(depth_t(zz)-depth_t(zz-1)));
end
if tt==1
tmp_ts2=mean(tmp_ts,2);
end 

tmp=(tmp_ts2(zz));
end  
% ------------------------------------ %

%%%%%%%%%% Backward Euler Loop %%%%%%%%%%%%%%
for iter=1:5 

% ### Functional Groups: ### %

% *** Nutrient Uptake Terms *** %

% --- Small Phytoplankton --- %
% Calculate Limitations
t_lim_v_sp = exp(tc_v_sp*tmp)*scale_tmp_phyto1;
l_lim_sp = (1.0-exp(-1.0*alpha_sp*rad/vmax_sp))*exp(-1.0*beta_sp*rad/vmax_sp);
NO_lim_sp = NO_v2/(NO_v2+k_NO_sp)*(1/(1+NH_v2/k_NH_sp));
NH_lim_sp = NH_v2/(NH_v2+k_NH_sp);

% sp NO Uptake 
delta_NO = NO_v1-(NO_v1/(1.0+dt/NO_v2*(vmax_sp*t_lim_v_sp...
*l_lim_sp*NO_lim_sp*sp_v2)));

% sp NH Uptake 
delta_NH = NH_v1-(NH_v1/(1.0+dt/NH_v2*(vmax_sp*t_lim_v_sp...      
*l_lim_sp*NH_lim_sp*sp_v2)));

% sp Excretion 
delta_sp = ext_excr_sp*(delta_NO+delta_NH);

% Update
NO_v2 = NO_v1-delta_NO;
NH_v2 = NH_v1-delta_NH;
DON_v2 = DON_v1+delta_sp;
sp_v2 = sp_v1+delta_NO+delta_NH-delta_sp;
% Note: v2 is updated here because with multiple loss terms its possible that 
% the loss terms could remove all of a state variable (i.e. grazing removes 70% and mortality removes 40%)

% Other
np_frac_sp = delta_NO/(max([min_val,(delta_NO+delta_NH)]));

% Diagnostic
gpp_sp=delta_NO+delta_NH;
spec_g_sp=vmax_sp*t_lim_v_sp*l_lim_sp*(NO_lim_sp+NH_lim_sp);
light_lim_sp=l_lim_sp;
nit_lim_sp=NO_lim_sp;
amm_lim_sp=NH_lim_sp;
nut_lim_sp=NO_lim_sp+NH_lim_sp;
excr_sp=delta_sp;
% ---------------------- %

% --- Large Phytoplankton ---%
% Calculate Limitations
t_lim_v_lp = exp(tc_v_lp*tmp)*scale_tmp_phyto1;
l_lim_lp = (1.0-exp(-1.0*alpha_lp*rad/vmax_lp))*exp(-1.0*beta_lp*rad/vmax_lp); 
NO_lim_lp = NO_v2/(NO_v2+k_NO_lp)*(1/(1+NH_v2/k_NH_lp));
NH_lim_lp = NH_v2/(NH_v2+k_NH_lp);
SI_lim_lp = SI_v2/(SI_v2+k_SI_lp);
NS_lim_comp = min([1.0,((SI_lim_lp)/(max([min_val,(NO_lim_lp+NH_lim_lp)])))]);

% lp NO Uptake 
delta_NO = NO_v2-(NO_v2/(1.0+dt/NO_v2*(vmax_lp*t_lim_v_lp*l_lim_lp...
*NO_lim_lp*NS_lim_comp*lp_v2)));

% lp NH Uptake 
delta_NH = NH_v2-(NH_v2/(1.0+dt/NH_v2*(vmax_lp*t_lim_v_lp*l_lim_lp...
*NH_lim_lp*NS_lim_comp*lp_v2)));

% lp Excretion
delta_lp = ext_excr_lp*(delta_NO+delta_NH);

% Update 
NO_v2 = NO_v2-delta_NO;
NH_v2 = NH_v2-delta_NH;
SI_v2 = SI_v1-r_SI_N*(delta_NO+delta_NH)+r_SI_N*delta_lp;
DON_v2 = DON_v2+delta_lp;
lp_v2 = lp_v1+delta_NO+delta_NH-delta_lp;

% Other 
np_frac_lp = delta_NO/(max([min_val,(delta_NO+delta_NH)]));

% Diagnostic
gpp_lp=delta_NO+delta_NH;
spec_g_lp=vmax_lp*t_lim_v_lp*l_lim_lp*(NO_lim_lp+NH_lim_lp)*NS_lim_comp;
light_lim_lp=l_lim_lp;
nit_lim_lp=NO_lim_lp;
amm_lim_lp=NH_lim_lp;
sil_lim_lp=SI_lim_lp;
nut_lim_lp=(NO_lim_lp+NH_lim_lp)*NS_lim_comp;
excr_lp=delta_lp;
% ---------------------- %

% *** END Nutrient Uptake Terms *** %

% ---- Phytoplankton Respiration & Mortality Terms --- %

% sp respiration
t_lim_r_sp = exp(tc_r_sp*tmp)*scale_tmp_phyto3;
delta_sp = sp_v2-(sp_v2/(1.0+dt/sp_v2*(ref_resp_sp*t_lim_r_sp*sp_v2)));

% lp respiration
t_lim_r_lp = exp(tc_r_lp*tmp)*scale_tmp_phyto3;
delta_lp = lp_v2-(lp_v2/(1.0+dt/lp_v2*(ref_resp_lp*t_lim_r_lp*lp_v2)));

% Update
sp_v2 = sp_v2-delta_sp;
lp_v2 = lp_v2-delta_lp;
NO_v2 = NO_v2+delta_sp*np_frac_sp+delta_lp*np_frac_lp;
NH_v2 = NH_v2+delta_sp*(1.0-np_frac_sp)+delta_lp*(1.0-np_frac_lp);
SI_v2 = SI_v2+delta_lp*r_SI_N;

% Diagnostic
phyto_resp=delta_sp+delta_lp;

% sp Mortality
t_lim_m_sp = exp(tc_m_sp*tmp)*scale_tmp_phyto2;
delta_sp = sp_v2-(sp_v2/(1+dt/sp_v2*(ref_mort_sp*t_lim_m_sp*sp_v2.^sp_mort_type)));

% lp Mortality 
t_lim_m_lp = exp(tc_m_lp*tmp)*scale_tmp_phyto2;
delta_lp = lp_v2-(lp_v2/(1+dt/lp_v2*(ref_mort_lp*t_lim_m_lp*lp_v2.^lp_mort_type)));

% Update
sp_v2 = sp_v2-delta_sp;
lp_v2 = lp_v2-delta_lp;
PON_v2 = PON_v1+delta_sp+delta_lp;
OP_v2 = OP_v1+delta_lp*r_SI_N;

% Diagnostics
phyto_mort=delta_sp+delta_lp;

% ------------------------- %

% ******* Zooplankton Grazing Terms ********  %

% --- Small Zooplankton Grazing --- %

% sz grazing on sp
t_lim_g_sz_sp = exp(tc_g_sz_sp*tmp)*scale_tmp_zoo1;
g_lim_sz_sp = max([0,(1-exp(iv_sz_sp*(thresh_sz_sp-sp_v2)))]);
l_lim_graz_sz_sp = rad/(rad+rad_graz_sz_sp);
delta_sp = sp_v2-(sp_v2/(1+dt/sp_v2*((gmax2_sz_sp+gmax_sz_sp*l_lim_graz_sz_sp)*t_lim_g_sz_sp*g_lim_sz_sp*sz_v2)));

% Update
sp_v2 = sp_v2-delta_sp;
sz_v2 = sz_v1+gge_sz*delta_sp;
NH_v2 = NH_v2+(ae_sz-gge_sz)*delta_sp;
PON_v2 = PON_v2+(1-ae_sz)*delta_sp;

% Diagnostics
phyto_graz_sz=delta_sp;
sz_graz=delta_sp;
graz_on_sp=delta_sp;
spec_g_sz=gmax_sz_sp*t_lim_g_sz_sp*g_lim_sz_sp;
sz_graz_lim=g_lim_sz_sp;
% ------------------------% 

% --- Large Zooplankton Grazing --- %

% lz grazing on sp
t_lim_g_lz_sp = exp(tc_g_lz_sp*tmp)*scale_tmp_zoo1;
g_lim_lz_sp = max([0,(1-exp(iv_lz_sp*(thresh_lz_sp-sp_v2)))]);
delta_sp = sp_v2-(sp_v2/(1+dt/sp_v2*(gmax_lz_sp*t_lim_g_lz_sp*g_lim_lz_sp*lz_v2)));

% lz grazing on lp
t_lim_g_lz_lp =exp(tc_g_lz_lp*tmp)*scale_tmp_zoo1;
g_lim_lz_lp = max([0,(1-exp(iv_lz_lp*(thresh_lz_lp-lp_v2)))]);
delta_lp = lp_v2-(lp_v2/(1+dt/lp_v2*(gmax_lz_lp*t_lim_g_lz_lp*g_lim_lz_lp*lz_v2)));

% lz grazing on sz
t_lim_g_lz_sz = exp(tc_g_lz_sz*tmp)*scale_tmp_zoo1;
g_lim_lz_sz = max([0,(1-exp(iv_lz_sz*(thresh_lz_sz-sz_v2)))]);
delta_sz = sz_v2-(sz_v2/(1+dt/sz_v2*(gmax_lz_sz*t_lim_g_lz_sz*g_lim_lz_sz*lz_v2)));

% Update
sp_v2 = sp_v2-delta_sp;
lp_v2 = lp_v2-delta_lp;
sz_v2 = sz_v2-delta_sz;
lz_v2 = lz_v1+gge_lz*(delta_sp+delta_lp+delta_sz);
NH_v2 = NH_v2+(ae_lz-gge_lz)*(delta_sp+delta_lp+delta_sz);
PON_v2 = PON_v2+(1-ae_lz)*(delta_sp+delta_lp+delta_sz);
OP_v2 = OP_v2+delta_lp*r_SI_N;

% Diagnostics
phyto_graz_lz=delta_sp+delta_lp;
lz_graz=phyto_graz_lz+delta_sz;
graz_on_sp=graz_on_sp+delta_sp;
graz_on_lp=delta_lp;
spec_g_lz=(gmax_lz_sp*t_lim_g_lz_sp*g_lim_lz_sp)+(gmax_lz_lp*t_lim_g_lz_lp*g_lim_lz_lp)+(gmax_lz_sz*t_lim_g_lz_sz*g_lim_lz_sz);
lz_graz_lim=g_lim_lz_sp+g_lim_lz_lp+g_lim_lz_sz;
% -------------------------% 

% --- Predatory Zooplankton Grazing --- % 

% pz grazing on lp
t_lim_g_pz_lp = exp(tc_g_pz_lp*tmp)*scale_tmp_zoo1;
g_lim_pz_lp = max([0,((1-exp(iv_pz_lp*(thresh_pz_lp-lp_v2)))*exp(-inh_szlz_lp*(sz_v2+lz_v2)))]);
delta_lp = lp_v2-(lp_v2/(1+dt/lp_v2*(gmax_pz_lp*t_lim_g_pz_lp*g_lim_pz_lp*pz_v2)));

% pz Grazing on sz
t_lim_g_pz_sz = exp(tc_g_pz_sz*tmp)*scale_tmp_zoo1;
g_lim_pz_sz = max([0,((1-exp(iv_pz_sz*(thresh_pz_sz-sz_v2)))*exp(-inh_lz_sz*lz_v2))]);
delta_sz = sz_v2-(sz_v2/(1+dt/sz_v2*(gmax_pz_sz*t_lim_g_pz_sz*g_lim_pz_sz*pz_v2)));

% pz Grazing on lz
t_lim_g_pz_lz = exp(tc_g_pz_lz*tmp)*scale_tmp_zoo1;
g_lim_pz_lz = max([0,(1-exp(iv_pz_lz*(thresh_pz_lz-lz_v2)))]);
delta_lz = lz_v2-(lz_v2/(1+dt/lz_v2*(gmax_pz_lz*t_lim_g_pz_lz*g_lim_pz_lz*pz_v2)));

% Update
lp_v2 = lp_v2-delta_lp;
sz_v2 = sz_v2-delta_sz;
lz_v2 = lz_v2-delta_lz;
pz_v2 = pz_v1+gge_pz*(delta_lp+delta_sz+delta_lz);
NH_v2 = NH_v2+(ae_pz-gge_pz)*(delta_lp+delta_sz+delta_lz);
PON_v2 = PON_v2+(1-ae_pz)*(delta_lp+delta_sz+delta_lz);
OP_v2 = OP_v2+delta_lp*r_SI_N;

% Diagnostics
phyto_graz_pz=delta_lp;
pz_graz=delta_lp+delta_sz+delta_lz;
graz_on_lp=graz_on_lp+delta_lp;
spec_g_pz=(gmax_pz_lp*t_lim_g_pz_lp*g_lim_pz_lp)+(gmax_pz_sz*t_lim_g_pz_sz*g_lim_pz_sz)+(gmax_pz_lz*t_lim_g_pz_lz*g_lim_pz_lz);
pz_graz_lim=g_lim_pz_lp+g_lim_pz_sz+g_lim_pz_lz;
% -----------------------------%

% *** END of Grazing Terms *** %

% *** Zooplankton Mortality Terms *** % 

% sz Mortality 
t_lim_m_sz = exp(tc_m_sz*tmp)*scale_tmp_zoo2;
delta_sz = sz_v2-(sz_v2/(1+dt/sz_v2*(ref_mort_sz*t_lim_m_sz*sz_v2.^sz_mort_type)));

% lz Mortality 
t_lim_m_lz = exp(tc_m_lz*tmp)*scale_tmp_zoo2;
delta_lz = lz_v2-(lz_v2/(1+dt/lz_v2*(ref_mort_lz*t_lim_m_lz*lz_v2.^lz_mort_type)));

% lz Mortality 
t_lim_m_pz = exp(tc_m_pz*tmp)*scale_tmp_zoo2;
delta_pz = pz_v2-(pz_v2/(1+dt/pz_v2*(ref_mort_pz*t_lim_m_pz*pz_v2.^pz_mort_type)));

% Update 
sz_v2 = sz_v2-delta_sz;
lz_v2 = lz_v2-delta_lz;
pz_v2 = pz_v2-delta_pz;
PON_v2 = PON_v2+delta_sz+delta_lz+delta_pz;
% ----------------- %

% *** END of Mortality Terms *** %

% ### END Functional Groups ### % 

% ### Chemical Transformations ### %

% Nitrification
t_lim_nit =exp(tc_nitr*tmp)*scale_tmp_decomp;
delta_NH = NH_v2-(NH_v2/(1+dt/NH_v2*(ref_nitr*t_lim_nit*NH_v2)));
NH_v2 = NH_v2-delta_NH;
NO_v2 = NO_v2+delta_NH;

% Decomp of PON to NH 
t_lim_dec_PON_NH = exp(tc_dec_PON_NH*tmp)*scale_tmp_decomp;
delta_PON = PON_v2-(PON_v2/(1+dt/PON_v2*(ref_dec_PON_NH*t_lim_dec_PON_NH*PON_v2)));
PON_v2 = PON_v2-delta_PON;
NH_v2 = NH_v2+delta_PON;

% Decomp of DON to NH 
t_lim_dec_DON_NH = exp(tc_dec_DON_NH*tmp)*scale_tmp_decomp;
delta_DON = DON_v2-(DON_v2/(1+dt/DON_v2*(ref_dec_DON_NH*t_lim_dec_DON_NH*DON_v2)));
DON_v2 = DON_v2-delta_DON;
NH_v2 = NH_v2+delta_DON;

% Decomp of PON to DON
t_lim_dec_PON_DON = exp(tc_dec_PON_DON*tmp)*scale_tmp_decomp;
delta_PON = PON_v2-(PON_v2/(1+dt/PON_v2*(ref_dec_PON_DON*t_lim_dec_PON_DON*PON_v2)));
PON_v2 = PON_v2-delta_PON;
DON_v2 = DON_v2+delta_PON;

% Decomp of OP to SI (First OP Loss Term)
t_lim_dec_OP_SI = exp(tc_dec_OP_SI*tmp)*scale_tmp_decomp;
delta_OP = OP_v2-(OP_v2/(1+dt/OP_v2*(ref_dec_OP_SI*t_lim_dec_OP_SI*OP_v2)));
OP_v2 = OP_v2-delta_OP;
SI_v2 = SI_v2+delta_OP;

% ### END Chemical Transformations ### %

end % Iter for loop 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Update 
NO_all(zz,tt+1)=NO_v2; NH_all(zz,tt+1)=NH_v2;
PON_all(zz,tt+1)=PON_v2; DON_all(zz,tt+1)=DON_v2;
sp_all(zz,tt+1)=sp_v2; lp_all(zz,tt+1)=lp_v2;
sz_all(zz,tt+1)=sz_v2; lz_all(zz,tt+1)=lz_v2;
pz_all(zz,tt+1)=pz_v2; SI_all(zz,tt+1)=SI_v2;
OP_all(zz,tt+1)=OP_v2;
% Note: The final update is done outside the backward Euler loop

% Check for NaN or Zeros
x=[NO_v2 NH_v2 PON_v2 DON_v2 sp_v2 lp_v2 sz_v2 lz_v2 pz_v2...
SI_v2 OP_v2];
ind=find(isnan(x) | x<0);
if isempty(ind)==0
disp('NaN Error')
display(['NO_v2(tt+1) = ', num2str(NO_all(zz,tt+1))])
display(['NH_v2(tt+1) = ', num2str(NH_all(zz,tt+1))])
display(['SI_v2(tt+1) = ', num2str(SI_all(zz,tt+1))])
display(['DON_v2(tt+1) = ', num2str(DON_all(zz,tt+1))])
display(['PON_v2(tt+1) = ', num2str(PON_all(zz,tt+1))])
display(['OP_v2(tt+1) = ', num2str(OP_all(zz,tt+1))])
display(['sp_v2(tt+1) = ', num2str(sp_all(zz,tt+1))])
display(['lp_v2(tt+1) = ', num2str(lp_all(zz,tt+1))])
display(['sz_v2(tt+1) = ', num2str(sz_all(zz,tt+1))])
display(['lz_v2(tt+1) = ', num2str(lz_all(zz,tt+1))])
display(['pz_v2(tt+1) = ', num2str(pz_all(zz,tt+1))])
display('')
display(['NO_v2(tt) = ', num2str(NO_all(zz,tt))])
display(['NH_v2(tt) = ', num2str(NH_all(zz,tt))])
display(['SI_v2(tt) = ', num2str(SI_all(zz,tt))])
display(['DON_v2(tt) = ', num2str(DON_all(zz,tt))])
display(['PON_v2(tt) = ', num2str(PON_all(zz,tt))])
display(['OP_v2(tt) = ', num2str(OP_all(zz,tt))])
display(['sp_v2(tt) = ', num2str(sp_all(zz,tt))])
display(['lp_v2(tt) = ', num2str(lp_all(zz,tt))])
display(['sz_v2(tt) = ', num2str(sz_all(zz,tt))])
display(['lz_v2(tt) = ', num2str(lz_all(zz,tt))])
display(['pz_v2(tt) = ', num2str(pz_all(zz,tt))])
pause 
break
end

% **** Save Diagnostics **** &
gpp_all(zz,tt)=gpp_sp+gpp_lp; 
phyto_resp_all(zz,tt)=phyto_resp;
phyto_mort_all(zz,tt)=phyto_mort;
phyto_graz_all(zz,tt)=phyto_graz_sz+phyto_graz_lz+phyto_graz_pz;
graz_on_sp_all(zz,tt)=graz_on_sp;
graz_on_lp_all(zz,tt)=graz_on_lp;
spec_g_sp_all(zz,tt)=spec_g_sp;
spec_g_lp_all(zz,tt)=spec_g_lp;
light_lim_sp_all(zz,tt)=light_lim_sp;
light_lim_lp_all(zz,tt)=light_lim_lp;
nut_lim_sp_all(zz,tt)=nut_lim_sp;
tmp_lim_v_sp_all(zz,tt)=t_lim_v_sp;
nut_lim_lp_all(zz,tt)=nut_lim_lp;
ns_lim_comp_all(zz,tt)=NS_lim_comp;
tmp_lim_v_lp_all(zz,tt)=t_lim_v_lp;
spec_g_sz_all(zz,tt)=spec_g_sz;
spec_g_lz_all(zz,tt)=spec_g_lz;
spec_g_pz_all(zz,tt)=spec_g_pz;
sz_graz_lim_all(zz,tt)=sz_graz_lim;
lz_graz_lim_all(zz,tt)=lz_graz_lim;
pz_graz_lim_all(zz,tt)=pz_graz_lim;
phyto_excr_all(zz,tt)=excr_sp+excr_lp;
sz_graz_all(zz,tt)=sz_graz;
lz_graz_all(zz,tt)=lz_graz;
pz_graz_all(zz,tt)=pz_graz;

% ****** Chl:C Sub Model ******* %
% Small Phytoplankton 
t_lim_v_sp = exp(tc_v_sp*tmp);
NO_lim_sp = NO_v2/(NO_v2+k_NO_sp)*(1/(1+NH_v2/k_NH_sp));
NH_lim_sp = NH_v2/(NH_v2+k_NH_sp);
v_chl_sp=vmax_sp*t_lim_v_sp*(NO_lim_sp+NH_lim_sp)*(alpha_sp/(alpha_sp+beta_sp))*(beta_sp/(alpha_sp+beta_sp))^(beta_sp/alpha_sp);
chl2c_sp(zz,tt)=chl2c_sp_max/(1.0+0.5*chl2c_sp_max*alpha_chl_sp/v_chl_sp*rad);

if chl2c_sp(zz,tt)<chl2c_sp_min
chl2c_sp(zz,tt)=chl2c_sp_min;
end 

% Large Phytoplankton
t_lim_v_lp = exp(tc_v_lp*tmp);
NO_lim_lp = NO_v2/(NO_v2+k_NO_lp)*(1/(1+NH_v2/k_NH_lp));
NH_lim_lp = NH_v2/(NH_v2+k_NH_lp);
SI_lim_lp = SI_v2/(SI_v2+k_SI_lp);
NS_lim_comp = min([1.0,((SI_lim_lp)/(max([min_val,(NO_lim_lp+NH_lim_lp)])))]);
v_chl_lp=vmax_lp*t_lim_v_lp*(NO_lim_lp+NH_lim_lp)*NS_lim_comp*(alpha_lp/(alpha_lp+beta_lp))*(beta_lp/(alpha_lp+beta_lp))^(beta_lp/alpha_lp);
chl2c_lp(zz,tt)=chl2c_lp_max/(1.0+0.5*chl2c_lp_max*alpha_chl_lp/v_chl_lp*rad);

if chl2c_lp(zz,tt)<chl2c_lp_min
chl2c_lp(zz,tt)=chl2c_lp_min;
end 

end 				% zz (End NEMURO spatial loop) 

% ****** Sinking ****** %
if sink_flag==1
PON_v1=PON_all(:,tt+1); PON_v2=[];
OP_v1=OP_all(:,tt+1); OP_v2=[];

for zz=1:length(depth_t)
if zz==1
PON_v2(zz,1)=PON_v1(zz)-PON_v1(zz)*sink*dt/dz(zz);
OP_v2(zz,1)=OP_v1(zz)-OP_v1(zz)*sink*dt/dz(zz);
end 
if zz==length(depth_t)
PON_v2(zz,1)=(PON_v1(zz)*vol(zz)+(PON_v1(zz-1)*vol(zz-1)*sink*dt/dz(zz-1)))/vol(zz);
OP_v2(zz,1)=(OP_v1(zz)*vol(zz)+(OP_v1(zz-1)*vol(zz-1)*sink*dt/dz(zz-1)))/vol(zz);
end
if zz>1 & zz<length(depth_t)
PON_v2(zz,1)=(PON_v1(zz)*vol(zz)-(PON_v1(zz)*vol(zz)*sink*dt/dz(zz))+(PON_v1(zz-1)*vol(zz-1)*sink*dt/dz(zz-1)))/vol(zz);
OP_v2(zz,1)=(OP_v1(zz)*vol(zz)-(OP_v1(zz)*vol(zz)*sink*dt/dz(zz))+(OP_v1(zz-1)*vol(zz-1)*sink*dt/dz(zz-1)))/vol(zz);
end 
end 

% Update 
PON_all(:,tt+1)=PON_v2;
OP_all(:,tt+1)=OP_v2;
end

% ****** Mixing ******* %
if mix_flag==1
NO_v1=NO_all(:,tt+1); NO_v2=[];
NH_v1=NH_all(:,tt+1); NH_v2=[];
SI_v1=SI_all(:,tt+1); SI_v2=[];
DON_v1=DON_all(:,tt+1); DON_v2=[];
PON_v1=PON_all(:,tt+1); PON_v2=[];
OP_v1=OP_all(:,tt+1); OP_v2=[];
sp_v1=sp_all(:,tt+1); sp_v2=[];
lp_v1=lp_all(:,tt+1); lp_v2=[];
sz_v1=sz_all(:,tt+1); sz_v2=[];
lz_v1=lz_all(:,tt+1); lz_v2=[];
pz_v1=pz_all(:,tt+1); pz_v2=[];

for zz=1:length(depth_t)
if zz==1
NO_sbc=NO_v1(zz); NH_sbc=NH_v1(zz);		% surface boundary condtion
SI_sbc=SI_v1(zz); DON_sbc=DON_v1(zz);		% surface boundary condtion
PON_sbc=PON_v1(zz); OP_sbc=OP_v1(zz);		% surface boundary condtion
sp_sbc=sp_v1(zz); lp_sbc=lp_v1(zz);		% surface boundary condtion
sz_sbc=sz_v1(zz); lz_sbc=lz_v1(zz);		% surface boundary condtion
pz_sbc=pz_v1(zz); dz_sbc=dz_t(zz);		% surface boundary condtion
NO_v2(zz)=NO_v1(zz)+dt*kz_day*((NO_v1(zz+1)-NO_v1(zz))/dz_t(zz) - (NO_v1(zz)-NO_sbc)/dz_sbc)/dz_w(zz);
NH_v2(zz)=NH_v1(zz)+dt*kz_day*((NH_v1(zz+1)-NH_v1(zz))/dz_t(zz) - (NH_v1(zz)-NH_sbc)/dz_sbc)/dz_w(zz);
SI_v2(zz)=SI_v1(zz)+dt*kz_day*((SI_v1(zz+1)-SI_v1(zz))/dz_t(zz) - (SI_v1(zz)-SI_sbc)/dz_sbc)/dz_w(zz);
DON_v2(zz)=DON_v1(zz)+dt*kz_day*((DON_v1(zz+1)-DON_v1(zz))/dz_t(zz) - (DON_v1(zz)-DON_sbc)/dz_sbc)/dz_w(zz);
PON_v2(zz)=PON_v1(zz)+dt*kz_day*((PON_v1(zz+1)-PON_v1(zz))/dz_t(zz) - (PON_v1(zz)-PON_sbc)/dz_sbc)/dz_w(zz);
OP_v2(zz)=OP_v1(zz)+dt*kz_day*((OP_v1(zz+1)-OP_v1(zz))/dz_t(zz) - (OP_v1(zz)-OP_sbc)/dz_sbc)/dz_w(zz);
sp_v2(zz)=sp_v1(zz)+dt*kz_day*((sp_v1(zz+1)-sp_v1(zz))/dz_t(zz) - (sp_v1(zz)-sp_sbc)/dz_sbc)/dz_w(zz);
lp_v2(zz)=lp_v1(zz)+dt*kz_day*((lp_v1(zz+1)-lp_v1(zz))/dz_t(zz) - (lp_v1(zz)-lp_sbc)/dz_sbc)/dz_w(zz);
sz_v2(zz)=sz_v1(zz)+dt*kz_day*((sz_v1(zz+1)-sz_v1(zz))/dz_t(zz) - (sz_v1(zz)-sz_sbc)/dz_sbc)/dz_w(zz);
lz_v2(zz)=lz_v1(zz)+dt*kz_day*((lz_v1(zz+1)-lz_v1(zz))/dz_t(zz) - (lz_v1(zz)-lz_sbc)/dz_sbc)/dz_w(zz);
pz_v2(zz)=pz_v1(zz)+dt*kz_day*((pz_v1(zz+1)-pz_v1(zz))/dz_t(zz) - (pz_v1(zz)-pz_sbc)/dz_sbc)/dz_w(zz);
end 
if zz==length(depth_t)
NO_bbc=NO_v1(zz);NH_bbc=NH_v1(zz);		% bottom boundary condtion
SI_bbc=SI_v1(zz);DON_bbc=DON_v1(zz);		% bottom boundary condtion
PON_bbc=PON_v1(zz);OP_bbc=OP_v1(zz);		% bottom boundary condtion
sp_bbc=sp_v1(zz);lp_bbc=lp_v1(zz);		% bottom boundary condtion
sz_bbc=sz_v1(zz);lz_bbc=lz_v1(zz);		% bottom boundary condtion
pz_bbc=pz_v1(zz); dz_bbc=dz_t(zz-1);		% bottom boundary condtion
NO_v2(zz)=NO_v1(zz)+dt*kz_day*((NO_bbc-NO_v1(zz))/dz_bbc - (NO_v1(zz)-NO_v1(zz-1))/dz_t(zz-1))/dz_w(zz);
NH_v2(zz)=NH_v1(zz)+dt*kz_day*((NH_bbc-NH_v1(zz))/dz_bbc - (NH_v1(zz)-NH_v1(zz-1))/dz_t(zz-1))/dz_w(zz);
SI_v2(zz)=SI_v1(zz)+dt*kz_day*((SI_bbc-SI_v1(zz))/dz_bbc - (SI_v1(zz)-SI_v1(zz-1))/dz_t(zz-1))/dz_w(zz);
DON_v2(zz)=DON_v1(zz)+dt*kz_day*((DON_bbc-DON_v1(zz))/dz_bbc - (DON_v1(zz)-DON_v1(zz-1))/dz_t(zz-1))/dz_w(zz);
PON_v2(zz)=PON_v1(zz)+dt*kz_day*((PON_bbc-PON_v1(zz))/dz_bbc - (PON_v1(zz)-PON_v1(zz-1))/dz_t(zz-1))/dz_w(zz);
OP_v2(zz)=OP_v1(zz)+dt*kz_day*((OP_bbc-OP_v1(zz))/dz_bbc - (OP_v1(zz)-OP_v1(zz-1))/dz_t(zz-1))/dz_w(zz);
sp_v2(zz)=sp_v1(zz)+dt*kz_day*((sp_bbc-sp_v1(zz))/dz_bbc - (sp_v1(zz)-sp_v1(zz-1))/dz_t(zz-1))/dz_w(zz);
lp_v2(zz)=lp_v1(zz)+dt*kz_day*((lp_bbc-lp_v1(zz))/dz_bbc - (lp_v1(zz)-lp_v1(zz-1))/dz_t(zz-1))/dz_w(zz);
sz_v2(zz)=sz_v1(zz)+dt*kz_day*((sz_bbc-sz_v1(zz))/dz_bbc - (sz_v1(zz)-sz_v1(zz-1))/dz_t(zz-1))/dz_w(zz);
lz_v2(zz)=lz_v1(zz)+dt*kz_day*((lz_bbc-lz_v1(zz))/dz_bbc - (lz_v1(zz)-lz_v1(zz-1))/dz_t(zz-1))/dz_w(zz);
pz_v2(zz)=pz_v1(zz)+dt*kz_day*((pz_bbc-pz_v1(zz))/dz_bbc - (pz_v1(zz)-pz_v1(zz-1))/dz_t(zz-1))/dz_w(zz);
end 
if  zz>1 & zz<length(depth_t)
NO_v2(zz)=NO_v1(zz)+dt*kz_day*((NO_v1(zz+1)-NO_v1(zz))/dz_t(zz) - (NO_v1(zz)-NO_v1(zz-1))/dz_t(zz-1))/dz_w(zz);
NH_v2(zz)=NH_v1(zz)+dt*kz_day*((NH_v1(zz+1)-NH_v1(zz))/dz_t(zz) - (NH_v1(zz)-NH_v1(zz-1))/dz_t(zz-1))/dz_w(zz);
SI_v2(zz)=SI_v1(zz)+dt*kz_day*((SI_v1(zz+1)-SI_v1(zz))/dz_t(zz) - (SI_v1(zz)-SI_v1(zz-1))/dz_t(zz-1))/dz_w(zz);
DON_v2(zz)=DON_v1(zz)+dt*kz_day*((DON_v1(zz+1)-DON_v1(zz))/dz_t(zz) - (DON_v1(zz)-DON_v1(zz-1))/dz_t(zz-1))/dz_w(zz);
PON_v2(zz)=PON_v1(zz)+dt*kz_day*((PON_v1(zz+1)-PON_v1(zz))/dz_t(zz) - (PON_v1(zz)-PON_v1(zz-1))/dz_t(zz-1))/dz_w(zz);
OP_v2(zz)=OP_v1(zz)+dt*kz_day*((OP_v1(zz+1)-OP_v1(zz))/dz_t(zz) - (OP_v1(zz)-OP_v1(zz-1))/dz_t(zz-1))/dz_w(zz);
sp_v2(zz)=sp_v1(zz)+dt*kz_day*((sp_v1(zz+1)-sp_v1(zz))/dz_t(zz) - (sp_v1(zz)-sp_v1(zz-1))/dz_t(zz-1))/dz_w(zz);
lp_v2(zz)=lp_v1(zz)+dt*kz_day*((lp_v1(zz+1)-lp_v1(zz))/dz_t(zz) - (lp_v1(zz)-lp_v1(zz-1))/dz_t(zz-1))/dz_w(zz);
sz_v2(zz)=sz_v1(zz)+dt*kz_day*((sz_v1(zz+1)-sz_v1(zz))/dz_t(zz) - (sz_v1(zz)-sz_v1(zz-1))/dz_t(zz-1))/dz_w(zz);
lz_v2(zz)=lz_v1(zz)+dt*kz_day*((lz_v1(zz+1)-lz_v1(zz))/dz_t(zz) - (lz_v1(zz)-lz_v1(zz-1))/dz_t(zz-1))/dz_w(zz);
pz_v2(zz)=pz_v1(zz)+dt*kz_day*((pz_v1(zz+1)-pz_v1(zz))/dz_t(zz) - (pz_v1(zz)-pz_v1(zz-1))/dz_t(zz-1))/dz_w(zz);
end 
end 
% Note original diffusion term:
% kz_day*(NO_v1(zz+1)+NO_v1(zz-1)-2*NO_v1(zz))/dz(zz)^2 
% only works if dz is constant in the domain that is why the mixing code above is written differently.

% Update
NO_all(:,tt+1)=NO_v2;
NH_all(:,tt+1)=NH_v2;
SI_all(:,tt+1)=SI_v2;
DON_all(:,tt+1)=DON_v2;
PON_all(:,tt+1)=PON_v2;
OP_all(:,tt+1)=OP_v2;
sp_all(:,tt+1)=sp_v2;
lp_all(:,tt+1)=lp_v2;
sz_all(:,tt+1)=sz_v2;
lz_all(:,tt+1)=lz_v2;
pz_all(:,tt+1)=pz_v2;

end
% ****************** %

% Check Conservation
temp1=NO_all(:,tt+1).*vol; temp2=NH_all(:,tt+1).*vol;
temp3=SI_all(:,tt+1).*vol; temp4=DON_all(:,tt+1).*vol;
temp5=PON_all(:,tt+1).*vol; temp6=OP_all(:,tt+1).*vol;
temp7=sp_all(:,tt+1).*vol; temp8=lp_all(:,tt+1).*vol;
temp88=lp_all(:,tt+1).*vol.*r_SI_N; temp9=sz_all(:,tt+1).*vol;
temp10=lz_all(:,tt+1).*vol; temp11=pz_all(:,tt+1).*vol;
tot_n(tt+1)=sum(temp1)+sum(temp2)+sum(temp4)+sum(temp5)+sum(temp7)+sum(temp8)+sum(temp9)+sum(temp10)+sum(temp11);
tot_s(tt+1)=sum(temp3)+sum(temp6)+sum(temp88);

% Count Updates
if mod(tt,tsteps_inhour)==0
h_count=h_count+1;
end 
if mod(tt,tsteps_inday)==0
d_count=d_count+1;
display(['Day Number = ', num2str(tt/tsteps_inday)]);
end 
if d_count==366
d_count=1;
h_count=1;
end

end 				 % tt (end temporal loop)
toc
%  >>>> End Main Loop <<<<< % 
% ##################################### %
% ############ END NEMURO ############# %
% ##################################### %

display('############# Simulation End ################### ');
display(' ')

% --------- Process Diagnostics ---------- %
% Chl
n2c_conv=106/16*1/1000*12.0107*1000;
chl2c_sp2(:,2:size(chl2c_sp,2)+1)=chl2c_sp;chl2c_sp2(:,1)=chl2c_sp(:,1);chl2c_sp=chl2c_sp2;clear chl2c_sp2; % make array sizes agree
chl2c_lp2(:,2:size(chl2c_lp,2)+1)=chl2c_lp;chl2c_lp2(:,1)=chl2c_lp(:,1);chl2c_lp=chl2c_lp2;clear chl2c_lp2; % make array sizes agree	
chl_sp=sp_all.*n2c_conv.*chl2c_sp;
chl_lp=lp_all.*n2c_conv.*chl2c_lp;
chl_sp2=sp_all.*n2c_conv.*chl2c_sp_constant;
chl_lp2=lp_all.*n2c_conv.*chl2c_lp_constant;
chl_all=chl_sp+chl_lp;
chl2_all=chl_sp2+chl_lp2;

% Create Daily Diagnostics 
ind=[0:tsteps_inday:size(gpp_all,2)];
gpp_day=[];phyto_resp_day=[];phyto_mort_day=[];phyto_graz_day=[];
spec_g_sp_day=[];spec_g_lp_day=[];spec_g_sz_day=[];spec_g_lz_day=[];
spec_g_pz_day=[];light_lim_sp_day=[];light_lim_lp_day=[];
nut_lim_sp_day=[];nut_lim_lp_day=[];sz_graz_lim_day=[];
lz_graz_lim_day=[];pz_graz_lim_day=[];sz_graz_day=[];lz_graz_day=[];
pz_graz_day=[];chl_day=[];NO_day=[];NH_day=[];SI_day=[];DON_day=[];
PON_day=[];OP_day=[];sp_day=[];lp_day=[];sz_day=[];lz_day=[];pz_day=[];
chl2_day=[];tmp_lim_v_sp_day=[];tmp_lim_v_lp_day=[];
graz_on_sp_day=[]; graz_on_lp_day=[]; ns_lim_comp_day=[];

for ii=1:length(ind)-1
for zz=1:size(gpp_all,1)
gpp_day(zz,ii)=sum(gpp_all(zz,ind(ii)+1:ind(ii+1)));
phyto_resp_day(zz,ii)=sum(phyto_resp_all(zz,ind(ii)+1:ind(ii+1)));
phyto_mort_day(zz,ii)=sum(phyto_mort_all(zz,ind(ii)+1:ind(ii+1)));
phyto_graz_day(zz,ii)=sum(phyto_graz_all(zz,ind(ii)+1:ind(ii+1)));
phyto_excr_day(zz,ii)=sum(phyto_excr_all(zz,ind(ii)+1:ind(ii+1)));
spec_g_sp_day(zz,ii)=mean(spec_g_sp_all(zz,ind(ii)+1:ind(ii+1)));
spec_g_lp_day(zz,ii)=mean(spec_g_lp_all(zz,ind(ii)+1:ind(ii+1)));
spec_g_sz_day(zz,ii)=mean(spec_g_sz_all(zz,ind(ii)+1:ind(ii+1)));
spec_g_lz_day(zz,ii)=mean(spec_g_lz_all(zz,ind(ii)+1:ind(ii+1)));
spec_g_pz_day(zz,ii)=mean(spec_g_pz_all(zz,ind(ii)+1:ind(ii+1)));
light_lim_sp_day(zz,ii)=mean(light_lim_sp_all(zz,ind(ii)+1:ind(ii+1)));
light_lim_lp_day(zz,ii)=mean(light_lim_lp_all(zz,ind(ii)+1:ind(ii+1)));
nut_lim_sp_day(zz,ii)=mean(nut_lim_sp_all(zz,ind(ii)+1:ind(ii+1)));
nut_lim_lp_day(zz,ii)=mean(nut_lim_lp_all(zz,ind(ii)+1:ind(ii+1)));
sz_graz_lim_day(zz,ii)=mean(sz_graz_lim_all(zz,ind(ii)+1:ind(ii+1)));
lz_graz_lim_day(zz,ii)=mean(lz_graz_lim_all(zz,ind(ii)+1:ind(ii+1)));
pz_graz_lim_day(zz,ii)=mean(pz_graz_lim_all(zz,ind(ii)+1:ind(ii+1)));
sz_graz_day(zz,ii)=mean(sz_graz_all(zz,ind(ii)+1:ind(ii+1)));
lz_graz_day(zz,ii)=mean(lz_graz_all(zz,ind(ii)+1:ind(ii+1)));
pz_graz_day(zz,ii)=mean(pz_graz_all(zz,ind(ii)+1:ind(ii+1)));
NO_day(zz,ii)=mean(NO_all(zz,ind(ii)+1:ind(ii+1)));
NH_day(zz,ii)=mean(NH_all(zz,ind(ii)+1:ind(ii+1)));
SI_day(zz,ii)=mean(SI_all(zz,ind(ii)+1:ind(ii+1)));
DON_day(zz,ii)=mean(DON_all(zz,ind(ii)+1:ind(ii+1)));
PON_day(zz,ii)=mean(PON_all(zz,ind(ii)+1:ind(ii+1)));
OP_day(zz,ii)=mean(OP_all(zz,ind(ii)+1:ind(ii+1)));
sp_day(zz,ii)=mean(sp_all(zz,ind(ii)+1:ind(ii+1)));
lp_day(zz,ii)=mean(lp_all(zz,ind(ii)+1:ind(ii+1)));
sz_day(zz,ii)=mean(sz_all(zz,ind(ii)+1:ind(ii+1)));
lz_day(zz,ii)=mean(lz_all(zz,ind(ii)+1:ind(ii+1)));
pz_day(zz,ii)=mean(pz_all(zz,ind(ii)+1:ind(ii+1)));
chl_day(zz,ii)=mean(chl_all(zz,ind(ii)+1:ind(ii+1)));
chl2_day(zz,ii)=mean(chl2_all(zz,ind(ii)+1:ind(ii+1)));
tmp_lim_v_sp_day(zz,ii)=mean(tmp_lim_v_sp_all(zz,ind(ii)+1:ind(ii+1)));
tmp_lim_v_lp_day(zz,ii)=mean(tmp_lim_v_lp_all(zz,ind(ii)+1:ind(ii+1)));
graz_on_sp_day(zz,ii)=mean(graz_on_sp_all(zz,ind(ii)+1:ind(ii+1)));
graz_on_lp_day(zz,ii)=mean(graz_on_lp_all(zz,ind(ii)+1:ind(ii+1)));
ns_lim_comp_day(zz,ii)=mean(ns_lim_comp_all(zz,ind(ii)+1:ind(ii+1)));
end 
end 

display('Surface State Variable Average Values:');
display(['NO = ',num2str(mean(NO_all(1,start_diag1:end_diag1)))]);
display(['NH = ',num2str(mean(NH_all(1,start_diag1:end_diag1)))]);
display(['SI = ',num2str(mean(SI_all(1,start_diag1:end_diag1)))]);
display(['DON = ',num2str(mean(DON_all(1,start_diag1:end_diag1)))]);
display(['PON = ',num2str(mean(PON_all(1,start_diag1:end_diag1)))]);
display(['OP = ',num2str(mean(OP_all(1,start_diag1:end_diag1)))]);
display(['sp = ',num2str(mean(sp_all(1,start_diag1:end_diag1)))]);
display(['lp = ',num2str(mean(lp_all(1,start_diag1:end_diag1)))]);
display(['sz = ',num2str(mean(sz_all(1,start_diag1:end_diag1)))]);
display(['lz = ',num2str(mean(lz_all(1,start_diag1:end_diag1)))]);
display(['pz = ',num2str(mean(pz_all(1,start_diag1:end_diag1)))]);
temp1= mean(lp_all(1,start_diag1:end_diag1));
temp2= mean(sp_all(1,start_diag1:end_diag1));
temp=temp1/(temp1+temp2);
display('-----');
display(['% Large Phytoplankton =  ', num2str(temp)]);

% Conservation
figure(20)
plot([1:tt+1],tot_n/tot_n(1)*100,'k','linewidth',3)
axis([1 tt+1 90 110])
hold on
plot([1:tt+1],tot_s/tot_s(1)*100,'r','linewidth',3)
title('Conservation Check')
legend('Total N', 'Total Si')
xlabel('Time Step')
ylabel('Percent Inital')

% Average Profiles
prof_gpp=mean(gpp_day(:,start_diag2:end_diag2),2);
prof_resp=mean(phyto_resp_day(:,start_diag2:end_diag2),2);
prof_mort=mean(phyto_mort_day(:,start_diag2:end_diag2),2);
prof_graz=mean(phyto_graz_day(:,start_diag2:end_diag2),2);
prof_excr=mean(phyto_excr_day(:,start_diag2:end_diag2),2);
spec_g_sp_prof=mean(spec_g_sp_day(:,start_diag2:end_diag2),2);
spec_g_lp_prof=mean(spec_g_lp_day(:,start_diag2:end_diag2),2);
spec_g_sz_prof=mean(spec_g_sz_day(:,start_diag2:end_diag2),2);
spec_g_lz_prof=mean(spec_g_lz_day(:,start_diag2:end_diag2),2);
spec_g_pz_prof=mean(spec_g_pz_day(:,start_diag2:end_diag2),2);
nut_lim_sp_prof=mean(nut_lim_sp_day(:,start_diag2:end_diag2),2);
nut_lim_lp_prof=mean(nut_lim_lp_day(:,start_diag2:end_diag2),2);
light_lim_sp_prof=mean(light_lim_sp_day(:,start_diag2:end_diag2),2);
light_lim_lp_prof=mean(light_lim_lp_day(:,start_diag2:end_diag2),2);
sz_graz_lim_prof=mean(sz_graz_lim_day(:,start_diag2:end_diag2),2);
lz_graz_lim_prof=mean(lz_graz_lim_day(:,start_diag2:end_diag2),2);
pz_graz_lim_prof=mean(pz_graz_lim_day(:,start_diag2:end_diag2),2);
sz_graz_prof=mean(sz_graz_day(:,start_diag2:end_diag2),2);
lz_graz_prof=mean(lz_graz_day(:,start_diag2:end_diag2),2);
pz_graz_prof=mean(pz_graz_day(:,start_diag2:end_diag2),2);
graz_on_sp_prof=mean(graz_on_sp_day(:,start_diag2:end_diag2),2);
graz_on_lp_prof=mean(graz_on_lp_day(:,start_diag2:end_diag2),2);
ns_lim_comp_prof=mean(ns_lim_comp_day(:,start_diag2:end_diag2),2);
NO_prof=mean(NO_all(:,start_diag1:end_diag1),2);
NH_prof=mean(NH_all(:,start_diag1:end_diag1),2);
SI_prof=mean(SI_all(:,start_diag1:end_diag1),2);
DON_prof=mean(DON_all(:,start_diag1:end_diag1),2);
PON_prof=mean(PON_all(:,start_diag1:end_diag1),2);
OP_prof=mean(OP_all(:,start_diag1:end_diag1),2);
sp_prof=mean(sp_all(:,start_diag1:end_diag1),2);
lp_prof=mean(lp_all(:,start_diag1:end_diag1),2);
sz_prof=mean(sz_all(:,start_diag1:end_diag1),2);
lz_prof=mean(lz_all(:,start_diag1:end_diag1),2);
pz_prof=mean(pz_all(:,start_diag1:end_diag1),2);
chl_prof=mean(chl_all(:,start_diag1:end_diag1),2);
chl2_prof=mean(chl2_all(:,start_diag1:end_diag1),2);
NO_prof_norm=NO_prof/NO_prof(1,1);
NH_prof_norm=NH_prof/NH_prof(1,1);
SI_prof_norm=SI_prof/SI_prof(1,1);
DON_prof_norm=DON_prof/DON_prof(1,1);
PON_prof_norm=PON_prof/PON_prof(1,1);
OP_prof_norm=OP_prof/OP_prof(1,1);
sp_prof_norm=sp_prof/sp_prof(1,1);
lp_prof_norm=lp_prof/lp_prof(1,1);
sz_prof_norm=sz_prof/sz_prof(1,1);
lz_prof_norm=lz_prof/lz_prof(1,1);
pz_prof_norm=pz_prof/pz_prof(1,1);

% Zoo 
zoo_all=pz_all+lz_all;
zoo_prof=mean(zoo_all(:,start_diag1:end_diag1),2);
ind=find(depth_t<=200);
zoo_200=mean(zoo_prof(ind));

% Convert to m^2
prof_gpp_C=prof_gpp.*n2c_conv.*dz;
prof_resp_C=prof_resp.*n2c_conv.*dz;
prof_mort_C=prof_mort.*n2c_conv.*dz;
prof_graz_C=prof_graz.*n2c_conv.*dz;
prof_excr_C=prof_excr.*n2c_conv.*dz;

% Rates
figure(1)
subplot('position',[0.1 0.15 0.4 0.8])
hold on
set(gca,'ydir','reverse')
plot(prof_gpp_C,depth_t,'g','linewidth',3)
plot(prof_excr_C,depth_t,'m','linewidth',3)
plot(prof_resp_C,depth_t,'k','linewidth',3)
plot(prof_mort_C,depth_t,'r','linewidth',3)
plot(prof_graz_C,depth_t,'linewidth',3)
xlabel('Rates')
ylabel('Depth')
axis([0 250 0 200])
legend('gpp','excr','resp','mort','graz')
xlabel('mg C m^-^2 d^-^1')
ylabel('Depth')
box on

growth_graz_sp_prof=spec_g_sp_prof-graz_on_sp_prof./sp_prof;
growth_graz_lp_prof=spec_g_lp_prof-graz_on_lp_prof./lp_prof;

subplot('position',[0.55 0.15 0.4 0.8])
plot(spec_g_sp_prof,depth_t,'g','linewidth',3)
hold on
plot(spec_g_lp_prof,depth_t,'linewidth',3)
plot(spec_g_sz_prof,depth_t,'k','linewidth',3)
plot(spec_g_lz_prof,depth_t,'r','linewidth',3)
plot(spec_g_pz_prof,depth_t,'m','linewidth',3)
set(gca,'ydir','reverse')
xlabel('Specific Growth d^-^1')
ylabel('Depth')
axis([0 1.5 0 200])
leg=legend('SP','LP','SZ','LZ','PZ');
box on
if print_flag==1
print('-dpng','-r75',[path3,'/1_rates.png']);
cd (path3)
!convert -trim 1_rates.png 1_rates.png
end

% Normalized Rates 
prof_resp_norm=prof_resp./prof_gpp*100;
prof_mort_norm=prof_mort./prof_gpp*100;
prof_graz_norm=prof_graz./prof_gpp*100;
prof_excr_norm=prof_excr./prof_gpp*100;

figure(2)
subplot('position',[0.1 0.15 0.4 0.8])
plot(prof_resp_norm,depth_t,'k','linewidth',3)
hold on
plot(prof_mort_norm,depth_t,'r','linewidth',3)
plot(prof_graz_norm,depth_t,'linewidth',3)
plot(prof_excr_norm,depth_t,'m','linewidth',3)
set(gca,'ydir','reverse')
xlabel('normalized rates')
ylabel('Depth')
axis([0 120 0 200])
xlabel('% GPP')
ylabel('Depth')
leg=legend('resp','mort','graz','excr');
box on
subplot('position',[0.55 0.15 0.4 0.8])
prof_loss=prof_resp_norm+prof_mort_norm+prof_graz_norm+prof_excr_norm;
plot(prof_loss,depth_t,'k','linewidth',3)
hold on
plot(depth_t*0+100,depth_t,'--k')
xlabel('Sum of Loss Terms')
set(gca,'ydir','reverse')
axis([0 200 0 200])
box on
if print_flag==1
print('-dpng','-r75',[path3,'/2_rates_norm.png']);
cd (path3)
!convert -trim 2_rates_norm.png 2_rates_norm.png
end

% Chl Profile
figure(3)
subplot('position',[0.1 0.1 0.4 0.8])
plot(NO_prof,depth_t,'k','linewidth',3)
set(gca,'ydir','reverse')
axis([0 20 0 200])
xlabel('Nitrate mmol N m^-^3')
ylabel('Depth')
text(5,40,['Surface Nitrate = ',num2str(NO_prof(1,1))])
box on
subplot('position',[0.55 0.1 0.4 0.8])
plot(chl_prof,depth_t,'g','linewidth',3)
hold on
plot(chl2_prof,depth_t,'--g')
xlabel('Chl mg Chl m^-^3')
set(gca,'ydir','reverse')
axis([0 2.5 0 200])
text(0.75,35,['Surface Chl = ',num2str(chl_prof(1,1))])
box on
if print_flag==1
print('-dpng','-r75',[path3,'/3_nitrate_chl.png']);
cd (path3)
!convert -trim 3_nitrate_chl.png 3_nitrate_chl.png
end

display(['Surface Nitrate = ',num2str(NO_prof(1,1))]);
display(['Surface Chl = ',num2str(chl_prof(1,1))])

% Nutrient and Light Limitation
figure(4)
subplot('position',[0.1 0.1 0.4 0.8])
plot(nut_lim_sp_prof,depth_t,'g','linewidth',3)
hold on
plot(nut_lim_lp_prof,depth_t,'linewidth',3)
plot(ns_lim_comp_prof,depth_t,'--b')
set(gca,'ydir','reverse')
axis([0 1 0 200])
ylabel('Depth')
xlabel('Nutrient Limitation')
leg=legend('SP','LP');
box on
subplot('position',[0.55 0.1 0.4 0.8])
plot(light_lim_sp_prof,depth_t,'g','linewidth',3)
hold on
plot(light_lim_lp_prof,depth_t,'linewidth',3)
set(gca,'ydir','reverse')
xlabel('Light Limitation')
axis([0 1 0 200])

rad_test=[];
for zz=1:length(depth_t)
if zz==1
rad_test(zz)=max(light_ts)*par_frac*exp(-1*(ext_sp*sp_prof(zz)+ext_lp*lp_prof(zz)+ext_w)*depth_t(zz));
else
rad_test(zz)=rad_test(zz-1)*exp(-1*(ext_sp*sp_prof(zz)+ext_lp*lp_prof(zz)+ext_w)*(depth_t(zz)-depth_t(zz-1)));
end
end
l_lim_sp = (1.0-exp(-1.0.*alpha_sp.*rad_test./vmax_sp)).*exp(-1.0.*beta_sp.*rad_test./vmax_sp);
l_lim_lp = (1.0-exp(-1.0.*alpha_lp.*rad_test./vmax_lp)).*exp(-1.0.*beta_lp.*rad_test./vmax_lp);
light_prof_norm=rad_test/rad_test(1,1);
ind=find(light_prof_norm>=0.01);
text(0.4,180,['Euphotic Zone = ', num2str(depth_t(max(ind)))])
plot(light_prof_norm,depth_t,'--k')
plot(l_lim_sp,depth_t,'--g')
plot(l_lim_lp,depth_t,'b--')
box on 

if print_flag==1
print('-dpng','-r75',[path3,'/4_phyto_lim.png']);
cd (path3)
!convert -trim 4_phyto_lim.png 4_phyto_lim.png
end

% Zooplankton Limitation 
graz_tot=sz_graz_prof+lz_graz_prof+pz_graz_prof;
sz_graz_prof_norm=sz_graz_prof./graz_tot.*100;
lz_graz_prof_norm=lz_graz_prof./graz_tot.*100;
pz_graz_prof_norm=pz_graz_prof./graz_tot.*100;
sz_graz_prof_norm(find(isnan(sz_graz_prof_norm)))=0;
lz_graz_prof_norm(find(isnan(lz_graz_prof_norm)))=0;
pz_graz_prof_norm(find(isnan(pz_graz_prof_norm)))=0;

figure(5)
subplot('position',[0.1 0.1 0.4 0.8])
hold on
plot(sz_graz_lim_prof,depth_t,'k','linewidth',3)
plot(lz_graz_lim_prof,depth_t,'r','linewidth',3)
plot(pz_graz_lim_prof,depth_t,'m','linewidth',3)
set(gca,'ydir','reverse')
axis([0 1 0 200])
leg=legend('SZ','LZ','PZ');
xlabel('Grazing Limitation')
ylabel('Depth')
box on
subplot('position',[0.55 0.1 0.4 0.8])
plot(sz_graz_prof_norm,depth_t,'k','linewidth',3)
hold on
plot(lz_graz_prof_norm,depth_t,'r','linewidth',3)
plot(pz_graz_prof_norm,depth_t,'m','linewidth',3)
set(gca,'ydir','reverse')
axis([0 100 0 200])
xlabel('% Total Grazing')
box on
if print_flag==1
print('-dpng','-r75',[path3,'/5_zoo_lim.png']);
cd (path3)
!convert -trim 5_zoo_lim.png 5_zoo_lim.png
end

% Temperature Limitation
figure(6)
subplot('position',[0.1 0.1 0.4 0.8])
x=[0:33];
prof1=exp(tdep1*x)*scale_tmp_phyto1;
prof2=exp(tdep2*x)*scale_tmp_phyto2;
prof3=exp(tdep3*x)*scale_tmp_phyto3;
prof4=exp(tdep4*x)*scale_tmp_zoo1;
prof5=exp(tdep5*x)*scale_tmp_zoo2;
prof6=exp(tdep6*x)*scale_tmp_decomp;
plot(x,prof2,'m','linewidth',3)
hold on
plot(x,prof5,'c','linewidth',3)
plot(x,prof6,'g','linewidth',3)
plot(x,prof1,'k','linewidth',3)
plot(x,prof3,'linewidth',3)
plot(x,prof4,'r','linewidth',3)
axis([0 30 0 10])
ylabel('Limitation')
xlabel('Temperature')
leg=legend('Pmort','Zmort','Decomp','Growth','Resp','Graz');
set(leg,'location','northwest')
box on
subplot('position',[0.575 0.1 0.4 0.8])
x=mean(tmp_ts,2);
prof1=exp(tdep1*x)*scale_tmp_phyto1;
prof2=exp(tdep2*x)*scale_tmp_phyto2;
prof3=exp(tdep3*x)*scale_tmp_phyto3;
prof4=exp(tdep4*x)*scale_tmp_zoo1;
prof5=exp(tdep5*x)*scale_tmp_zoo2;
prof6=exp(tdep6*x)*scale_tmp_decomp;
plot(prof2,depth_t,'m','linewidth',3)
hold on
plot(prof5,depth_t,'c','linewidth',3)
plot(prof6,depth_t,'g','linewidth',3)
plot(prof1,depth_t,'k','linewidth',3)
plot(prof3,depth_t,'linewidth',3)
plot(prof4,depth_t,'r','linewidth',3)
set(gca,'ydir','reverse')
axis([1 10 0 200])
title(['Graz / Growth @ 95m = ', num2str(prof4(10,1)/prof1(10,1)*100),' %'])
xlabel('Temperature Limitation')
ylabel('Depth')
box on
if print_flag==1
print('-dpng','-r75',[path3,'/6_tmp_lim.png']);
cd (path3)
!convert -trim 6_tmp_lim.png 6_tmp_lim.png
end 

%  Nutrient Profiles 
figure(7)
subplot('position',[0.1 0.1 0.4 0.8])
plot(NO_prof,depth_t,'k','linewidth',3)
hold on
plot(NH_prof,depth_t,'r','linewidth',3)
plot(SI_prof,depth_t,'linewidth',3)
plot(DON_prof,depth_t,'y','linewidth',3)
plot(PON_prof,depth_t,'g','linewidth',3)
plot(OP_prof,depth_t,'m','linewidth',3)
set(gca,'ydir','reverse')
axis([0 15 0 200])
xlabel('Nutrients')
ylabel('Depth')
box on
subplot('position',[0.55 0.1 0.4 0.8])
hold on
plot(NO_prof_norm,depth_t,'k','linewidth',3)
plot(NH_prof_norm,depth_t,'r','linewidth',3)
plot(SI_prof_norm,depth_t,'linewidth',3)
plot(DON_prof_norm,depth_t,'y','linewidth',3)
plot(PON_prof_norm,depth_t,'g','linewidth',3)
plot(OP_prof_norm,depth_t,'m','linewidth',3)
set(gca,'ydir','reverse')
axis([0 200 0 200])
legend('NO','NH','SI','DON','PON','OP')
xlabel('Normalized to Surface')
box on
if print_flag==1
print('-dpng','-r75',[path3,'/7_nutr_prof.png']);
cd (path3)
!convert -trim 7_nutr_prof.png 7_nutr_prof.png
end

% Biomass Profiles
figure(8)
subplot('position',[0.1 0.1 0.4 0.8])
plot(sp_prof,depth_t,'g','linewidth',3)
hold on
plot(lp_prof,depth_t,'linewidth',3)
plot(sz_prof,depth_t,'k','linewidth',3)
plot(lz_prof,depth_t,'r','linewidth',3)
plot(pz_prof,depth_t,'m','linewidth',3)
set(gca,'ydir','reverse')
axis([0 0.5 0 200])
xlabel('Biomass')
text(0.15,100,['200m Zoo Avg = ', num2str(mean(zoo_200))]);
box on 
ylabel('Depth')
subplot('position',[0.55 0.1 0.4 0.8])
hold on
plot(sp_prof_norm,depth_t,'g','linewidth',3)
plot(lp_prof_norm,depth_t,'linewidth',3)
plot(sz_prof_norm,depth_t,'k','linewidth',3)
plot(lz_prof_norm,depth_t,'r','linewidth',3)
plot(pz_prof_norm,depth_t,'m','linewidth',3)
set(gca,'ydir','reverse')
axis([0 10 0 200])
legend('SP','LP','SZ','LZ','PZ')
xlabel('Normalized to Surface')
box on
if print_flag==1
print('-dpng','-r75',[path3,'/8_bio_prof.png']);
cd (path3)
!convert -trim 8_bio_prof.png 8_bio_prof.png
end

display(['Zoo 200m Avg = ', num2str(mean(zoo_200))])
display('Finished Main Figures') 
% ------------------------------------------------------------ %

more_flag=0;
if more_flag==1

% Time Series
NO_ts=NO_day(1,start_diag2:end_diag2); NH_ts=NH_day(1,start_diag2:end_diag2);
SI_ts=SI_day(1,start_diag2:end_diag2); DON_ts=DON_day(1,start_diag2:end_diag2);
PON_ts=PON_day(1,start_diag2:end_diag2); OP_ts=OP_day(1,start_diag2:end_diag2);
sp_ts=sp_day(1,start_diag2:end_diag2); lp_ts=lp_day(1,start_diag2:end_diag2);
sz_ts=sz_day(1,start_diag2:end_diag2); lz_ts=lz_day(1,start_diag2:end_diag2);
pz_ts=pz_day(1,start_diag2:end_diag2); chl_ts=chl_day(1,start_diag2:end_diag2);
chl2_ts=chl2_day(1,start_diag2:end_diag2);
NO_ts2=NO_ts/mean(NO_ts); NH_ts2=NH_ts/mean(NH_ts);
SI_ts2=SI_ts/mean(SI_ts); DON_ts2=DON_ts/mean(DON_ts);
PON_ts2=PON_ts/mean(PON_ts); OP_ts2=OP_ts/mean(OP_ts);
sp_ts2=sp_ts/mean(sp_ts); lp_ts2=lp_ts/mean(lp_ts);
sz_ts2=sz_ts/mean(sz_ts); lz_ts2=lz_ts/mean(lz_ts);
pz_ts2=pz_ts/mean(pz_ts); chl_ts2=chl_ts/mean(chl_ts);
chl2_ts2=chl2_ts/mean(chl2_ts);
nut_lim_sp_ts=nut_lim_sp_day(1,start_diag2:end_diag2);
nut_lim_lp_ts=nut_lim_lp_day(1,start_diag2:end_diag2);
light_lim_sp_ts=light_lim_sp_day(1,start_diag2:end_diag2);
light_lim_lp_ts=light_lim_lp_day(1,start_diag2:end_diag2);
tmp_lim_v_sp_ts=tmp_lim_v_sp_day(1,start_diag2:end_diag2);
tmp_lim_v_lp_ts=tmp_lim_v_lp_day(1,start_diag2:end_diag2);
sz_graz_lim_ts=sz_graz_lim_day(1,start_diag2:end_diag2);
lz_graz_lim_ts=lz_graz_lim_day(1,start_diag2:end_diag2);
pz_graz_lim_ts=pz_graz_lim_day(1,start_diag2:end_diag2);
spec_g_sp_ts=spec_g_sp_day(1,start_diag2:end_diag2);
spec_g_lp_ts=spec_g_lp_day(1,start_diag2:end_diag2);
spec_g_sz_ts=spec_g_sz_day(1,start_diag2:end_diag2);
spec_g_lz_ts=spec_g_lz_day(1,start_diag2:end_diag2);
spec_g_pz_ts=spec_g_pz_day(1,start_diag2:end_diag2);
spec_g_sp_ts_norm=spec_g_sp_ts/mean(spec_g_sp_ts);
spec_g_lp_ts_norm=spec_g_lp_ts/mean(spec_g_lp_ts);
spec_g_sz_ts_norm=spec_g_sz_ts/mean(spec_g_sz_ts);
spec_g_lz_ts_norm=spec_g_lz_ts/mean(spec_g_lz_ts);
spec_g_pz_ts_norm=spec_g_pz_ts/mean(spec_g_pz_ts);
tot_lim_sp=nut_lim_sp_ts.*light_lim_sp_ts.*tmp_lim_v_sp_ts;
tot_lim_lp=nut_lim_lp_ts.*light_lim_lp_ts.*tmp_lim_v_lp_ts;
graz_on_sp_ts=graz_on_sp_day(1,start_diag2:end_diag2);
graz_on_lp_ts=graz_on_lp_day(1,start_diag2:end_diag2);
time=[1:length(NO_ts)];

% Nutrient Time Series 
figure(9)
subplot('position',[0.1 0.55 0.8 0.4])
hold on
plot(time,NO_ts,'k','linewidth',3)
plot(time,NH_ts,'r','linewidth',3)
plot(time,SI_ts,'linewidth',3)
plot(time,DON_ts,'y','linewidth',3)
plot(time,PON_ts,'g','linewidth',3)
plot(time,OP_ts,'m','linewidth',3)
ylabel('mmol N/m3')
axis([1 365 0 0.5])
set(gca,'xticklabel',[]);
box on
title('Surface (0-10m)  Nutrient Time Series')
leg=legend('NO','NH','SI','DON','PON','OP');
set(leg,'orientation','horizontal')
set(leg,'location','north')
subplot('position',[0.1 0.1 0.8 0.4])
hold on
plot(time,NO_ts2,'k','linewidth',3)
plot(time,NH_ts2,'r','linewidth',3)
plot(time,SI_ts2,'linewidth',3)
plot(time,DON_ts2,'y','linewidth',3)
plot(time,PON_ts2,'g','linewidth',3)
plot(time,OP_ts2,'m','linewidth',3)
plot(time,time*0+1,'--k')
xlabel('Days')
ylabel('Scaled 1/mean value')
axis([1 365 0 4])
box on
if print_flag==1
print('-dpng','-r75',[path3,'/9_nutr_ts.png']);
cd (path3)
!convert -trim 9_nutr_ts.png 9_nutr_ts.png
end

% Box Biomass Time Series
figure(10)
subplot('position',[0.1 0.55 0.8 0.4])
hold on
plot(time,sp_ts,'g','linewidth',3)
plot(time,lp_ts,'linewidth',3)
plot(time,sz_ts,'k','linewidth',3)
plot(time,lz_ts,'r','linewidth',3)
plot(time,pz_ts,'m','linewidth',3)
ylabel('mmol N/m3')
leg=legend('sp','lp','sz','lz','pz');
set(leg,'orientation','horizontal')
axis([1 365 0 0.25])
set(gca,'xticklabel',[]);
box on
title('Surface Biomass Time Series')
subplot('position',[0.1 0.1 0.8 0.4])
hold on
plot(time,sp_ts2,'g','linewidth',3)
plot(time,lp_ts2,'linewidth',3)
plot(time,sz_ts2,'k','linewidth',3)
plot(time,lz_ts2,'r','linewidth',3)
plot(time,pz_ts2,'m','linewidth',3)
plot(time,time*0+1,'--k')
ylabel('Scaled 1/mean value')
xlabel('Days')
axis([1 365 0.5 2.5])
box on
if print_flag==1
print('-dpng','-r75',[path3,'/10_bio_ts.png']);
cd (path3)
!convert -trim 10_bio_ts.png 10_bio_ts.png
end

% Chla Time Series
ind=find(depth_t<=100);
chl_100m_ts=chl_day(ind,start_diag2:end_diag2);
chl_100m_ts=mean(chl_100m_ts,1);
chl_100m_ts_norm=chl_100m_ts./mean(chl_100m_ts);
chl2_100m_ts=chl2_day(ind,start_diag2:end_diag2);
chl2_100m_ts=mean(chl2_100m_ts,1);
chl2_100m_ts_norm=chl2_100m_ts./mean(chl2_100m_ts);
chl_ts_norm=chl_ts./mean(chl_ts);
chl2_ts_norm=chl2_ts./mean(chl2_ts);

temp=[1:5:366];
chl_ts_norm2=[];
time2=[];
for ii=1:length(temp)-1
chl_ts_norm2(1,ii)=mean(chl_ts_norm(temp(ii):temp(ii+1)-1));
time2(1,ii)=mean([temp(ii),temp(ii+1)-1]);
end

figure(11)
subplot('position',[0.1 0.55 0.8 0.4])
hold on
plot(time,chl2_ts,'k','linewidth',3)
plot(time,chl_ts,'r','linewidth',3)
plot(time,chl2_100m_ts,'k')
plot(time,chl_100m_ts,'r')
leg=legend('constant Chl:C','Chl:C model');
ylabel('mg Chl/m^3')
axis([1 365 0 1.5])
box on
title('Chl Time Series')
subplot('position',[0.1 0.1 0.8 0.4])
plot(time,chl2_ts_norm,'k','linewidth',3)
hold on
plot(time2,chl_ts_norm2,'r','linewidth',3)
plot(time,chl_100m_ts_norm,'r')
plot(time,chl2_100m_ts_norm,'k')
plot(time,time*0+1,'--k')
axis([1 365 0 3])
box on
if print_flag==1
print('-dpng','-r75',[path3,'/11_chl_ts.png']);
cd (path3)
!convert -trim 11_chl_ts.png 11_chl_ts.png
end

% Zooplankton Time Series
ind=find(depth_t<=200);
zoo_200m_ts=lz_day(ind,:)+pz_day(ind,:);
zoo_200m_ts=mean(zoo_200m_ts,1);
zoo_200m_ts=zoo_200m_ts(1,start_diag2:end_diag2);
lz_200m_ts=mean(lz_day(ind,:),1);
lz_200m_ts=lz_200m_ts(1,start_diag2:end_diag2);
pz_200m_ts=mean(pz_day(ind,:),1);
pz_200m_ts=pz_200m_ts(1,start_diag2:end_diag2);
zoo_200m_ts_norm=zoo_200m_ts./mean(zoo_200m_ts);
lz_200m_ts_norm=lz_200m_ts./mean(lz_200m_ts);
pz_200m_ts_norm=pz_200m_ts./mean(pz_200m_ts);

figure(12)
subplot('position',[0.1 0.55 0.8 0.4])
hold on
plot(time,zoo_200m_ts,'k','linewidth',3)
plot(time,lz_200m_ts,'r','linewidth',3)
plot(time,pz_200m_ts,'m','linewidth',3)
xlabel('Day')
ylabel('mmol N m^-^3')
axis([1 365 0 0.12])
box on
legend('200m Zoo Avg', 'LZ', 'PZ')
title('Zooplankton Time Series')
subplot('position',[0.1 0.1 0.8 0.4])
hold on
plot(time,zoo_200m_ts_norm,'k','linewidth',3)
plot(time,lz_200m_ts_norm,'r','linewidth',3)
plot(time,pz_200m_ts_norm,'m','linewidth',3)
plot(time,time*0+1,'--k')
axis([1 365 0.5 1.5])
box on
if print_flag==1
print('-dpng','-r75',[path3,'/12_zoo_ts.png']);
cd (path3)
!convert -trim 12_zoo_ts.png 12_zoo_ts.png
end

% Phytoplankton Limitaiton TS
nut_lim_sp_ts_norm=nut_lim_sp_ts/mean(nut_lim_sp_ts);
nut_lim_lp_ts_norm=nut_lim_lp_ts/mean(nut_lim_lp_ts);
light_lim_sp_ts_norm=light_lim_sp_ts/mean(light_lim_sp_ts);
light_lim_lp_ts_norm=light_lim_lp_ts/mean(light_lim_lp_ts);

temp=[1:5:366];
ts1=[];
ts2=[];
ts3=[];
for ii=1:length(temp)-1
ts1(1,ii)=mean(light_lim_sp_ts_norm(temp(ii):temp(ii+1)-1));
ts2(1,ii)=mean(light_lim_lp_ts_norm(temp(ii):temp(ii+1)-1));
ts3(1,ii)=mean([temp(ii),temp(ii+1)-1]);
end 

figure(13)
subplot('position',[0.1 0.55 0.8 0.4])
plot(time,nut_lim_sp_ts,'g','linewidth',3)
hold on
plot(time,nut_lim_lp_ts,'linewidth',3)
plot(time,light_lim_sp_ts,'color',[0 0.5 0],'linewidth',3)
plot(time,light_lim_lp_ts,'color',[0 0 0.5],'linewidth',3)
axis([1 365 0 1])
leg=legend('SP','LP');
set(leg,'orientation','horizontal')
title('Phytoplankton Limitation')
box on
subplot('position',[0.1 0.1 0.8 0.4])
plot(time,nut_lim_sp_ts_norm,'g','linewidth',3)
hold on
plot(time,nut_lim_lp_ts_norm,'linewidth',3)
plot(ts3,ts1,'color',[0 0.5 0],'linewidth',3)
plot(ts3,ts2,'color',[0 0 0.5],'linewidth',3)
plot(time,time*0+1,'--k')
axis([1 365 0.5 1.5])
box on
if print_flag==1
print('-dpng','-r75',[path3,'/13_phyto_lim_ts.png']);
cd (path3)
!convert -trim 13_phyto_lim_ts.png 13_phyto_lim_ts.png
end

% Zooplankton Limitation TS 
sz_graz_lim_ts_norm=sz_graz_lim_ts./mean(sz_graz_lim_ts);
lz_graz_lim_ts_norm=lz_graz_lim_ts./mean(lz_graz_lim_ts);
pz_graz_lim_ts_norm=pz_graz_lim_ts./mean(pz_graz_lim_ts);

figure(14)
subplot('position',[0.1 0.55 0.8 0.4])
plot(time,sz_graz_lim_ts,'k','linewidth',3)
hold on
plot(time,lz_graz_lim_ts,'r','linewidth',3)
plot(time,pz_graz_lim_ts,'m','linewidth',3)
axis([1 365 0 0.25]);
leg=legend('SZ','LZ','PZ');
set(leg,'orientation','horizontal')
box on
title('Zooplankton Limitation')
subplot('position',[0.1 0.1 0.8 0.4])
plot(time,sz_graz_lim_ts_norm,'k','linewidth',3)
hold on
plot(time,lz_graz_lim_ts_norm,'r','linewidth',3)
plot(time,pz_graz_lim_ts_norm,'m','linewidth',3)
plot(time,time*0+1,'--k')
axis([1 365 0 5])
box on
if print_flag==1
print('-dpng','-r75',[path3,'/14_zoo_lim_ts.png']);
cd (path3)
!convert -trim 14_zoo_lim_ts.png 14_zoo_lim_ts.png
end

figure(15)
subplot('position',[0.1 0.55 0.8 0.4])
plot(time,spec_g_sp_ts,'g','linewidth',3)
hold on
plot(time,spec_g_lp_ts,'linewidth',3)
plot(time,spec_g_sz_ts,'k','linewidth',3)
plot(time,spec_g_lz_ts,'r','linewidth',3)
plot(time,spec_g_pz_ts,'m','linewidth',3)
axis([1 365 0 0.25])
ylabel('Specific Growth Rate d^-1')
leg=legend('SP','LP','SZ','LZ','PZ');
set(leg,'location','north','orientation','horizontal')
box on
title('Specific Growth Rates')
subplot('position',[0.1 0.1 0.8 0.4])
plot(time,spec_g_sp_ts_norm,'g','linewidth',3)
hold on
plot(time,spec_g_lp_ts_norm,'linewidth',3)
plot(time,spec_g_sz_ts_norm,'k','linewidth',3)
plot(time,spec_g_lz_ts_norm,'r','linewidth',3)
plot(time,spec_g_pz_ts_norm,'m','linewidth',3)
plot(time,time*0+1,'--k')
axis([1 365 0 5])
xlabel('Day')
ylabel('Specific Growth Rate 1/mean')
box on
if print_flag==1
print('-dpng','-r75',[path3,'/15_spec_g_ts.png']);
cd (path3)
!convert -trim 15_spec_g_ts.png 15_spec_g_ts.png
end

end % more_flag

% ---------------------------------------------------------------------------

%whos NO_prof NH_prof SI_prof DON_prof PON_prof OP_prof sp_prof lp_prof sz_prof lz_prof pz_prof
%save bgcm_profiles NO_prof NH_prof SI_prof DON_prof PON_prof OP_prof sp_prof lp_prof sz_prof lz_prof pz_prof  




















