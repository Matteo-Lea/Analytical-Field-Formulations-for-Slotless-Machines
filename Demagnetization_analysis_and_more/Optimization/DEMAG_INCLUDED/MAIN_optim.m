%% PSO parameters limits and objectives definition

% Matteo Leandro
% matteo.leandro@ntnu.no
% VERSION: 

%
%% Machine parameters
clearvars
clc
close all
tic
% machine/magnets configuration parameters
el = 100; % population for magnets configuration (set it a multiple of 4)
    
% Definition of non-integers variables in the optimization
alpha_p_lim = [0.3 ,0.7];
l_a_lim = [1.2 ,5]; % possible scaling factor (L/D)
R_s_lim = [0.025, 0.15];
t_wind_lim = [0.0015, 0.06];
g_lim = [0.0007, 0.0007];
t_m_lim = [0.004, 0.04];
J_lim = [7*1e6, 7*1e6]*sqrt(2);

Var_min = [ alpha_p_lim(1), l_a_lim(1), R_s_lim(1), ...
            t_wind_lim(1), g_lim(1), t_m_lim(1), ...
                                    J_lim(1) ].*ones(el,1);
Var_max = [ alpha_p_lim(2), l_a_lim(2), R_s_lim(2), ...
            t_wind_lim(2), g_lim(2), t_m_lim(2), ...
                                    J_lim(2) ].*ones(el,1);
                                                
% Definition of the integers variables in the optimization
In_Out_lim = [1, 1];
In_Out_adm = 1; % number of admissible values for the variable In_Out
p_lim = [15, 15]; % singularity on p=2 was not solved!!!!
p_adm = p_lim(2)-p_lim(1)+1; % number of admissible values for the variable p

Int_min = [In_Out_lim(1), p_lim(1)].*ones(el,1);
Int_max = [In_Out_lim(2), p_lim(2)].*ones(el,1);

adm.Int = [In_Out_adm, p_adm];


% Definition of constant parameters
Halbach = 0;
Halbach_1 = 2;
B_r = 1.3112;
B_rs = B_r;
theta_m_end = 0;
theta_m_side = 0;
Const = [Halbach, Halbach_1, B_r, B_rs, theta_m_end, theta_m_side]...
                                                            .*ones(el,1); 
                                                        
% Winding/machine parameters
    params_mach.N_tc = 1; % number of conductors per coil
    params_mach.q = 1; %coils per pole and per phase
    params_mach.b = 1; % parallel ways per phase
    params_mach.K_fill = 0.6; % number of conductors per coil
    params_mach.omega = 3000*2*pi/60; % mechanical angular frequency [rad/s]
    params_mach.res_copp = 0.025; % copper resisitivity @100 C [Ohm mm2/m] 
    params_mach.rho_mag = 8300; % PM density [kg/m3] (NdFeB=7500;Ferrite=5000;BondedNdFeB=5100;SmCo=8400;Alnico=7300)
    params_mach.rho_iron = 7700; % iron density [kg/m3] 
    params_mach.rho_copper = 8960; % copper density [kg/m3] 
    params_mach.rho_epoxy = 1300; % epoxy density [kg/m3] 
    params_mach.rho_alu = 2700; % aluminum density [kg/m3]
    params_mach.B_sat = 1.8; % saturation point iron [T]
    params_mach.mu_r = 1.067; % PM relative permeability

% Optiization Parameters
params.Np = el;         % Population size
params.Nr = 200;        % Repository size
params.maxgen = 100;    % Maximum number of generations
params.W = 0.5;         % Inertia weight
params.C1 = 2;          % Individual confidence factor
params.C2 = 2;          % Swarm confidence factor
params.C3 = 1.5;  % Social cognitive coefficient
params.C4 = 1.2;  % self-cognitive coefficient
% params.C3 = 1+1.1/(p_adm);  % Social cognitive coefficient
% params.C4 = 1+1/(p_adm);  % self-cognitive coefficient
params.ngrid = 20;      % Number of grids in each dimension
params.maxvel = 1;      % Maxmium vel in percentage
params.u_mut = 0.5;     % Uniform mutation percentage
params.m_PM = 400; % number of harmonics or series components for the magnetization functions
params.r_dis = 55; % radial discretization of magnets region

Torque = 5; % target torque [Nm]
Back_emf = 50; % maximum back emf [V]
Dem = 1; % admissible demagnetized volume [%]
Constraints = [Torque, Back_emf, Dem]; % constraints

% REP = Optimizer(params,params_mach,Var_max,Var_min,Int_min,Int_max,Const,adm,Constraints);
REP = MAIN_PSO(params,params_mach,Var_max,Var_min,Int_min,Int_max,Const,adm,Constraints);

toc
% Display info
display('Repository fitness values are stored in REP.pos_fit');
display('Repository particles positions are store in REP.pos');
