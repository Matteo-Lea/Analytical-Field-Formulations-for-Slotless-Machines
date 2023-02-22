%% ANALYTICAL CORE LOSS CALCULATION (SCRIPT)

% Matteo Leandro
% matteo.leandro@ntnu.no
% VERSION: 12Feb2020

% This code finds the field solution in the stator core region of an 
% inrunner slotless PM machine, in terms of flux density distribution,
% and performs the iron losses analysis based on the input loss data
% given by the user

%% Machine parameters
% clearvars
% clc
% close all

rho_Fe = 7650; % steel specific weight (mass density)[kg/m^3]
mu_0 = 4*pi*1e-7; % air permeability

Inrunner

CORR = 0.82; % correction factor for axial leakage based on some measurement on a different motor
STACK = 0.97; % laminations stacking factor (from 0.93 to 0.97

% CORR = 1; % correction factor for axial leakage  
% STACK = 1; % laminations stacking factor (from 0.93 to 0.97

%% useful indices for series harmonics definition
m_PM = 10; % number of harmonics or series components tionfor the magnetization functions
x = 0:1:m_PM;
% x = linspace(0,40,m_PM+1);
n = p*(2*x+1); % existing harmonic series components (n=1,3,5,...)

%% circumferential discretization
m_th = 1000; % points along the circumference/arc
theta = linspace(0,2*pi/p,m_th);
Theta = repmat(theta,m_PM+1,1);

%% Magnetization distribution

% PARALLEL MAGNETIZED PM
% Coefficients for defining the magnetization distribution
if Halbach_1 == 2 && alpha_p ==0
% Mid-magnets coefficients
M_1_0 = 0;
M_2_0 = 0;
% Side-magnets coefficients
M_1_2 = cos(theta_m_side+delta_m)./((n-1)*pi/(2*p));
M_2_2 = cos(theta_m_side+delta_m)./((n+1)*pi/(2*p));
M_3_2 = cos((n-1)*pi/(2*p)+theta_m_side+delta_m)./((n-1)*pi/(2*p));
M_4_2 = cos((n+1)*pi/(2*p)-theta_m_side-delta_m)./((n+1)*pi/(2*p));
% End-side-magnets coefficients
M_1_3 = 0;
M_2_3 = 0;
M_3_3 = 0;
M_4_3 = 0;

else
% Mid-magnets coefficients
M_1_0 = sin((n-1)*alpha_p*pi/(2*p))./((n-1)*alpha_p*pi/(2*p));
M_2_0 = sin((n+1)*alpha_p*pi/(2*p))./((n+1)*alpha_p*pi/(2*p));
% Side-magnets coefficients
M_1_2 = cos((n-1)*alpha_p*pi/(2*p)+theta_m_side+delta_m)./((n-1)*alpha_p*pi/(2*p));
M_2_2 = cos((n+1)*alpha_p*pi/(2*p)-theta_m_side-delta_m)./((n+1)*alpha_p*pi/(2*p));
M_3_2 = cos((n-1)*alpha_p1*pi/(2*p)+theta_m_side+delta_m)./((n-1)*alpha_p1*pi/(2*p));
M_4_2 = cos((n+1)*alpha_p1*pi/(2*p)-theta_m_side-delta_m)./((n+1)*alpha_p1*pi/(2*p));
% End-Side-magnets coefficients
M_1_3 = cos((n-1)*alpha_p1*pi/(2*p)+theta_m_end+delta_m1)./((n-1)*alpha_p1*pi/(2*p));
M_2_3 = cos((n+1)*alpha_p1*pi/(2*p)-theta_m_end-delta_m1)./((n+1)*alpha_p1*pi/(2*p));
M_3_3 = cos((n-1)*pi/(2*p)+theta_m_end+delta_m1)./((n-1)*pi/(2*p));
M_4_3 = cos((n+1)*pi/(2*p)-theta_m_end-delta_m1)./((n+1)*pi/(2*p));

end
    

if p==1
   if Halbach_1 == 2 && alpha_p ==0
       M_1_2(1) = 0;
    %    M_1_2(1) = -sin(theta_m_side+delta_m);
       M_3_2(1) = -sin(theta_m_side+delta_m);
   else
       M_1_0(1) = 1; 
       M_1_2(1) = -sin(theta_m_side+delta_m);
       M_3_2(1) = -sin(theta_m_side+delta_m);
       M_1_3(1) = -sin(theta_m_end+delta_m1);
       M_3_3(1) = -sin(theta_m_end+delta_m1);
   end
end

% Magnetization harmonics amplitude
if R_w>R_m % inrunner case
        M_r_n_par_mid = B_r/mu_0*alpha_p*(M_1_0+M_2_0);
        M_theta_n_par_mid = -B_r/mu_0*alpha_p*(M_1_0-M_2_0);
        M_r_n_par_end_side = B_rs/mu_0*(alpha_p1*(M_1_3-M_2_3)+(M_4_3-M_3_3));
        M_theta_n_par_end_side = -B_rs/mu_0*(alpha_p1*(M_2_3+M_1_3)-(M_3_3+M_4_3));
        M_r_n_par_side = B_rs/mu_0*(alpha_p*(M_1_2-M_2_2)+alpha_p1*(M_4_2-M_3_2));
        M_theta_n_par_side = -B_rs/mu_0*(alpha_p*(M_2_2+M_1_2)-alpha_p1*(M_3_2+M_4_2));
elseif R_w<R_m % outrunner case
    if p == 1 && R_w<R_m % outrunner case
        M_r_n_par_mid = B_r/mu_0*alpha_p*(M_1_0+M_2_0);
        M_theta_n_par_mid = -B_r/mu_0*alpha_p*(M_1_0-M_2_0);
        M_r_n_par_end_side = B_rs/mu_0*(alpha_p1*(M_1_3-M_2_3)+(M_4_3-M_3_3));
        M_theta_n_par_end_side = -B_rs/mu_0*(alpha_p1*(M_2_3+M_1_3)-(M_3_3+M_4_3));
        M_r_n_par_side = -B_rs/mu_0*(alpha_p*(M_1_2-M_2_2)+alpha_p1*(M_4_2-M_3_2));
        M_theta_n_par_side = B_rs/mu_0*(alpha_p*(M_2_2+M_1_2)-alpha_p1*(M_3_2+M_4_2));
    elseif R_w<R_m % outrunner case
        M_r_n_par_mid = B_r/mu_0*alpha_p*(M_1_0+M_2_0);
        M_theta_n_par_mid = B_r/mu_0*alpha_p*(M_1_0-M_2_0);
        M_r_n_par_end_side = B_rs/mu_0*(alpha_p1*(M_1_3-M_2_3)+(M_4_3-M_3_3));
        M_theta_n_par_end_side = B_rs/mu_0*(alpha_p1*(M_2_3+M_1_3)-(M_3_3+M_4_3));
        M_r_n_par_side = B_rs/mu_0*(alpha_p*(M_1_2-M_2_2)+alpha_p1*(M_4_2-M_3_2));
        M_theta_n_par_side = B_rs/mu_0*(alpha_p*(M_2_2+M_1_2)-alpha_p1*(M_3_2+M_4_2));
    end
end
    
M_r_n = M_r_n_par_mid+M_r_n_par_end_side+M_r_n_par_side;
M_theta_n = (M_theta_n_par_mid+M_theta_n_par_end_side+M_theta_n_par_side);


%% HARMONICS COEFFICIENTS FROM THE MAGNETIC VECTOR POTENTIAL
% (WINDINGS/AIR-GAP REGIONS)
% coefficient from the particular solution of Poisson's equation accounting
% for the magnetization distribution
A_zm_n = -mu_0*(n.*M_r_n+M_theta_n)./(1-(n).^2);

if p==1
   
    A_zm_n(1) = -mu_0/2*(M_r_n(1)+M_theta_n(1));
    
end

% K_B_n = (((A_zm_n-n.*A_zm_n+mu_0*M_theta_n)-2*(R_r/R_m).^(n+1).*(A_zm_n+mu_0*M_theta_n)+(A_zm_n+n.*A_zm_n+mu_0*M_theta_n).*(R_r/R_m).^(2*n))./((mu_r+1)*(1-(R_r/R_s).^(2*n))-(mu_r-1)*((R_m/R_s).^(2*n)-(R_r/R_m).^(2*n))));

DEN_I = ((mu_r^2-1)*(-(R_i/R_r).^(2*n)-(R_m/R_s).^(2*n)+(R_r/R_s).^(2*n)+(R_i/R_m).^(2*n))+(mu_r+1)^2*(1-(R_i/R_s).^(2*n))+(mu_r-1)^2*((R_m/R_s).^(2*n).*(R_i/R_r).^(2*n)-(R_r/R_m).^(2*n)));


K_B_n = (((-A_zm_n+n.*A_zm_n-mu_0*M_theta_n).*((mu_r+1)-(R_i/R_r).^(2*n)*(mu_r-1))+(A_zm_n+n.*A_zm_n+mu_0*M_theta_n).*((R_r/R_m).^(2*n)*(mu_r-1)-(R_i/R_m).^(2*n)*(mu_r+1))+2*(R_r/R_m).^(n+1).*((A_zm_n+mu_0*M_theta_n-mu_r*n.*A_zm_n)+(R_i/R_r).^(2*n).*(A_zm_n+mu_0*M_theta_n+mu_r*n.*A_zm_n)))./DEN_I);


%% Setting up data for the iron loss analysis
n_rpm = 500; % rotational speed [rpm]
f = p*n_rpm/60; % fundamental frequency [Hz]
f_h = f*n./p; % harmonc frequency components

% LINEAR STATOR CORE DISCRETIZATION METHOD
dis = 40; % points along the radial direction of the stator core (R_s and R_s included);
r_dis = linspace(R_s,R_se,dis); % discretized stator core line
delta_r = r_dis(2)-r_dis(1); % width of each sub-segment
r_dis = linspace(R_s+delta_r/2,R_se-delta_r/2,dis-1); % mid-point radius for each sub-segment
dis_weight = l_a*rho_Fe*(pi*((r_dis+delta_r/2).^2-(r_dis-delta_r/2).^2)); % vector holding the weight of the annuli

B_rad_n = CORR*abs(2*K_B_n./((R_s/R_se).^(2*n)-1).*((r_dis'/R_se).^(n-1).*(R_m/R_se).^(n+1)-(R_m./r_dis)'.^(n+1)));
B_tan_n = CORR*abs(2*K_B_n./((R_s/R_se).^(2*n)-1).*((r_dis'/R_se).^(n-1).*(R_m/R_se).^(n+1)+(R_m./r_dis)'.^(n+1)));

% Iron losses contributions with the  linear stator core discretization
% method
B_max_r = max(B_rad_n*cos(n'.*Theta),[],2);
B_max_theta = max(-B_tan_n*sin(n'.*Theta),[],2);
B_max = sqrt(B_max_r.^2+B_max_theta.^2);

%% Iron loss analysis based on the three component equation CCM:
% Ke*B^2*f^2 + Kh*B^alpha*f + Ka*B^1.5*f^1.5

% the different coefficients and exponents are defined in the following for
% the sake of generality
COEFF = [1.534203933539936e-05 0.007612428379577 6.717665659241482e-04]; % namely: eddy currents(Ke), hysteresis(Kh), excess loss (Ka) coefficients 
B_exp = [2 4.606837286647030 1.5]; % flux density exponents from the above written equation [2 2.64688 1.5]
f_exp = [2 1 1.5]; % frequency exponents from the above written equation


B_rad_n = CORR*abs(2*K_B_n./((R_s/R_se).^(2*n)-1).*((r_dis'/R_se).^(n-1).*(R_m/R_se).^(n+1)-(R_m./r_dis)'.^(n+1)));
B_tan_n = CORR*abs(2*K_B_n./((R_s/R_se).^(2*n)-1).*((r_dis'/R_se).^(n-1).*(R_m/R_se).^(n+1)+(R_m./r_dis)'.^(n+1)));

% Iron losses contributions with the  linear stator core discretization
% method
% Eddy currents contribution
EC_loss_dis_rad = dis_weight*COEFF(1)*(B_rad_n).^B_exp(1)*(f_h.^f_exp(1))';
EC_loss_dis_tan = dis_weight*COEFF(1)*(B_tan_n).^B_exp(1)*(f_h.^f_exp(1))';


% Hysteresis contribution
B_max_r = max(B_rad_n*cos(n'.*Theta),[],2);
B_max_theta = max(-B_tan_n*sin(n'.*Theta),[],2);
HY_loss_dis_rad = dis_weight*COEFF(2)*(B_max_r).^B_exp(2)*f;
HY_loss_dis_tan = dis_weight*COEFF(2)*(B_max_theta).^B_exp(2)*f;

% Excess loss contribution
EX_loss_dis_rad = dis_weight*COEFF(3)*(B_rad_n).^B_exp(3)*(f_h.^f_exp(3))';
EX_loss_dis_tan = dis_weight*COEFF(3)*(B_tan_n).^B_exp(3)*(f_h.^f_exp(3))';

% time domain
EC_loss_dis_rad_diff = mean(dis_weight*COEFF(1)*abs(diff((B_rad_n*cos(n'.*Theta)),1,2)/(theta(2)-theta(1))).^B_exp(1)*(2*pi*f/p)^B_exp(1)/(2*pi^2));
EC_loss_dis_tan_diff = mean(dis_weight*COEFF(1)*abs(diff((-B_tan_n*sin(n'.*Theta)),1,2)/(theta(2)-theta(1))).^B_exp(1)*(2*pi*f/p)^B_exp(1)/(2*pi^2));

EX_loss_dis_rad_diff = mean(dis_weight*COEFF(3)*abs(diff((B_rad_n*cos(n'.*Theta)),1,2)/(theta(2)-theta(1))).^B_exp(3)*(2*pi*f/p)^B_exp(3)/8.7631);
EX_loss_dis_tan_diff = mean(dis_weight*COEFF(3)*abs(diff((-B_tan_n*sin(n'.*Theta)),1,2)/(theta(2)-theta(1))).^B_exp(3)*(2*pi*f/p)^B_exp(3)/8.7631);


CCM_loss_time = STACK*(EC_loss_dis_rad_diff+EC_loss_dis_tan_diff+HY_loss_dis_rad+HY_loss_dis_tan+EX_loss_dis_rad_diff+EX_loss_dis_tan_diff);
CCM_loss_freq = STACK*(EC_loss_dis_rad+EC_loss_dis_tan+HY_loss_dis_rad+HY_loss_dis_tan+EX_loss_dis_rad+EX_loss_dis_tan);


