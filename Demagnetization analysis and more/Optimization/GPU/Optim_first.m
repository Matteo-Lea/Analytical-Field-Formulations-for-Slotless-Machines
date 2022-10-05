%% ANALYTICAL FIELD SOLUTION FROM PM (SCRIPT)

% Matteo Leandro
% matteo.leandro@ntnu.no
% VERSION: 12Feb2020

% This code finds the field solution in the airgap region, magnets region
% and iron region of a slotless machine (inrunner or outrunner) and plots
% the field map in all the regions
%% Machine parameters
clearvars
clc
close all
tic
% machine/magnets configuration parameters
el = 100; % population for magnets configuration (set it a multiple of 4)

params_mag.Halbach = randi([0,1],el,1); % 0 for even segments and 1 for odd segmens arrays
params_mag.Halbach_1 = params_mag.Halbach+2*randi([1,2],el,1); % number of segments per pole
params_mag.B_r = 0.9+0.5*rand(el,1); % remanent flux density [T] (mid-magnet)
params_mag.B_rs = 0.9+0.5*rand(el,1); % remanent flux density [T] (side-magnets)
params_mag.alpha_p =  rand(el,1)*(0.9-0.1)+0.1; % mid-magnet to pole ratio [-]
params_mag.alpha_p1 =  rand(el,1).*(1-params_mag.alpha_p)+params_mag.alpha_p;% side-magnet + mid-magnet to pole ratio [-]
params_mag.theta_m_end = rand(el,1).*90;  % orientation angle end-side-magnet [deg]
params_mag.theta_m_side = rand(el,1).*(90-params_mag.theta_m_end)+params_mag.theta_m_end; % orientation angle side-magnet [deg]
params_mag.p = randi([1,25],el,1); % pole pairs
params_mag.In_Out = sign(rand(el,1).*(1+1)-1); % 1 for inrunner -1 for outrunner (this shouldn't vary in the optimizatio n process)

m_PM = 10; % number of harmonics or series components for the magnetization functions

[M_r_n,M_theta_n] = Magnetization(el,m_PM,params_mag);

% Winding/machine parameters
params_mach.N_tc = 1; % number of conductors per coil
params_mach.q = 1; %coils per pole and per phase
params_mach.b = 1; % parallel ways per phase
params_mach.I = 28; % phase current peak value [A]
params_mach.omega = 1000*2*pi/60*params_mag.p; % mechanical angular frequency [rad/s]
params_mach.rho_mag = 7500; % PM density [kg/m3] (NdFeB=7500;Ferrite=5000;BondedNdFeB=5100;SmCo=8400;Alnico=7300)
J_Ph_max = 10*1e6; % maximum peak current density [A/m2]
J_Ph_min = 3*1e6; % minimum peak current density [A/m2]

% parameters to be optimized
params_geo.In_Out = params_mag.In_Out; 
params_geo.l_a = rand(el,1).*(1-0.01)+0.01; % active lenght [m]
params_geo.R_s = rand(el,1).*(1-0.05)+0.05; % stator radius (facing winding) [m]
params_geo.t_bi = rand(el,1).*(0.1-0.001)+0.001; % stator back-iron thickness [m]
params_geo.t_wind = rand(el,1).*(params_geo.In_Out.*params_geo.R_s-params_geo.In_Out.*sqrt(params_geo.R_s.^2-params_mach.N_tc*params_mach.q/params_mach.b*params_mach.I*6*params_mag.p./(pi*J_Ph_min).*params_geo.In_Out)-(params_geo.In_Out.*params_geo.R_s-params_geo.In_Out.*sqrt(params_geo.R_s.^2-params_mach.N_tc*params_mach.q/params_mach.b*params_mach.I*6*params_mag.p./(pi*J_Ph_max).*params_geo.In_Out)))+params_geo.In_Out.*params_geo.R_s-params_geo.In_Out.*sqrt(params_geo.R_s.^2-params_mach.N_tc*params_mach.q/params_mach.b*params_mach.I*6*params_mag.p./(pi*J_Ph_max).*params_geo.In_Out); % stator winding thickness [m]
params_geo.g = rand(el,1).*(0.007-0.001)+0.001; % air-gap thickness [m]
params_geo.t_m = rand(el,1).*(0.1-0.001)+0.001; % magnets thickness [m]
params_geo.t_sleeve = rand(el,1).*(0.1-0.001)+0.001; % magnets sleeve thickness (when is there) [m]
% derived parameters
params_geo.R_sleeve = params_geo.R_s - params_geo.In_Out.*params_geo.t_m- params_geo.In_Out.*params_geo.g- params_geo.In_Out.*params_geo.t_wind- params_geo.In_Out.*params_geo.t_sleeve; % magnets array radius (facing magnets support)[m]
% The chosen initialization allows R_r to assume negative values. Theese
% values are therefore corrected
params_geo.R_s(params_geo.R_sleeve<0) = params_geo.R_s(params_geo.R_sleeve<0) + 2*abs(params_geo.R_sleeve(params_geo.R_sleeve<0));
params_geo.R_se = params_geo.R_s + params_geo.In_Out.*params_geo.t_bi; % outer stator radius [m]
params_geo.R_wi = params_geo.R_s; % winding radius (facing back-iron) [m]
params_geo.R_w = params_geo.R_wi - params_geo.In_Out.*params_geo.t_wind; % winding radius (facing air-gap)[m]
params_geo.R_m = params_geo.R_w - params_geo.In_Out.*params_geo.g; % magnets array radius (facing air-gap)[m]
params_geo.R_r = params_geo.R_m - params_geo.In_Out.*params_geo.t_m; % magnets array radius (facing magnets support)[m]

% This definition of R_i shouldn't vary iteratively in the optimiation
% process to keep the geometries with and without backing iron
params_geo.R_i = params_geo.R_r; % infinite permeability boundary behind the magnets [m] (magnetic support is default)
params_geo.R_i(params_geo.In_Out==1) = params_geo.R_i(params_geo.In_Out==1).*randi([0,1],nnz(params_geo.In_Out==1),1); % non-magnetic support set for some in-runner topologies (randomly)
params_geo.R_i(params_geo.In_Out==-1) = params_geo.R_i(params_geo.In_Out==-1).*1./randi([0,1],nnz(params_geo.In_Out==-1),1); % non-magnetic support set for some out-runner topologies (randomly)
params_geo.R_sleeve(:,1) = 0;
params_geo.R_sleeve(params_geo.R_i==params_geo.R_r) = params_geo.R_r(params_geo.R_i==params_geo.R_r)- params_geo.In_Out(params_geo.R_i==params_geo.R_r).*params_geo.t_sleeve(params_geo.R_i==params_geo.R_r); % external magnetic support radius [m]

params_geo.mu_r = ones(el,1); % PM recoil permeability [-]
params_geo.p = params_mag.p; % pole pairs
params_geo.In_Out = params_mag.In_Out; 



[E_n,T_avg,B_bi,B_sleeve] = Field_Solution(el,m_PM,params_geo,M_r_n,M_theta_n,params_mach);
toc