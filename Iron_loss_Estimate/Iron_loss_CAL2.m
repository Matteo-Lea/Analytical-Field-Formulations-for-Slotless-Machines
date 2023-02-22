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

    %% Iron loss analysis based on CAL2 model (frequency domain):

    % Iron loss analysis based on the two components equation with variable 
    % coefficients and constant flux density exponent for the hysteresis 
    % component:
    % c1(f,B)*B^2*f^2 + c2(f,B)*B^2*f 

    % The following part of the code, evaluates iron losses based on the above
    % written equation. The model is based on experimental loss data, which are
    % splitted in two frequency ranges (to get a better fit of the data
    % themselves). For this reason, the maximum frequency of the first range
    % (fr1) has to be set (the model can be extended to more ranges by
    % followingthe same coding fashion).
    % In each frequency range the the dependency on the frequency itself is
    % neglected.

    fr1 = 200; % break-point frequency
    f = p*n_rpm/60; % fundamental frequency [Hz]
    index_1 = round((fr1/f)/2);

    CAL2_loss1 = 0;

    if index_1>0
        n1 = n(1:index_1);
        f_h1 = f*n1./p; % harmonic frequency components

        Kh_pol1 = [0.0268410355056720,-0.00763928948562653,-0.0496283787257046,0.0952874581566632,-0.0662449579608217,0.0171928957187306];
        Ke_pol1 = [1.89140553722200e-05,0.000277575435533729,-0.000838867014213555,0.00113091388730795,-0.000722107419157703,0.000177850114015975];

        % Iron losses contributions with the  linear stator core discretization
        % method
        % Eddy currents contribution
        B_rad_n1 = CORR*abs(2*K_B_n(1:index_1)./((R_s/R_se).^(2*n1)-1).*((r_dis'/R_se).^(n1-1).*(R_m/R_se).^(n1+1)-(R_m./r_dis)'.^(n1+1)));
        B_tan_n1 = CORR*abs(2*K_B_n(1:index_1)./((R_s/R_se).^(2*n1)-1).*((r_dis'/R_se).^(n1-1).*(R_m/R_se).^(n1+1)+(R_m./r_dis)'.^(n1+1)));

        EC_loss_dis_rad1 = dis_weight*((Ke_pol1(1) +Ke_pol1(2)*B_rad_n1 +Ke_pol1(3)*B_rad_n1.^2 +Ke_pol1(4)*B_rad_n1.^3 +Ke_pol1(5)*B_rad_n1.^4 +Ke_pol1(6)*B_rad_n1.^5).*abs(B_rad_n1).^2)*(f_h1.^2)';
        EC_loss_dis_tan1 = dis_weight*((Ke_pol1(1) +Ke_pol1(2)*B_tan_n1 +Ke_pol1(3)*B_tan_n1.^2 +Ke_pol1(4)*B_tan_n1.^3 +Ke_pol1(5)*B_tan_n1.^4 +Ke_pol1(6)*B_tan_n1.^5).*abs(B_tan_n1).^2)*(f_h1.^2)';

        % Hysteresis contribution (single hysteresis loop)
        B_max_r = CORR*max((2*K_B_n./((R_s/R_se).^(2*n)-1).*((r_dis'/R_se).^(n-1).*(R_m/R_se).^(n+1)-(R_m./r_dis)'.^(n+1)))*cos(n'.*Theta),[],2);
        B_max_theta = CORR*max(-(2*K_B_n./((R_s/R_se).^(2*n)-1).*((r_dis'/R_se).^(n-1).*(R_m/R_se).^(n+1)+(R_m./r_dis)'.^(n+1)))*sin(n'.*Theta),[],2);
        HY_loss_dis_rad_main = dis_weight*((Kh_pol1(1) +Kh_pol1(2)*B_max_r +Kh_pol1(3)*B_max_r.^2 +Kh_pol1(4)*B_max_r.^3 +Kh_pol1(5)*B_max_r.^4 +Kh_pol1(6)*B_max_r.^5).*(B_max_r).^2)*f;
        HY_loss_dis_tan_main = dis_weight*((Kh_pol1(1) +Kh_pol1(2)*B_max_theta +Kh_pol1(3)*B_max_theta.^2 +Kh_pol1(4)*B_max_theta.^3 +Kh_pol1(5)*B_max_theta.^4 +Kh_pol1(6)*B_max_theta.^5).*(B_max_theta).^2)*f;

        % Hysteresis contribution (individual loop per harmonic)
        HY_loss_dis_rad1 = dis_weight*((Kh_pol1(1) +Kh_pol1(2)*B_rad_n1 +Kh_pol1(3)*B_rad_n1.^2 +Kh_pol1(4)*B_rad_n1.^3 +Kh_pol1(5)*B_rad_n1.^4 +Kh_pol1(6)*B_rad_n1.^5).*(B_rad_n1).^2)*(f_h1)';
        HY_loss_dis_tan1 = dis_weight*((Kh_pol1(1) +Kh_pol1(2)*B_tan_n1 +Kh_pol1(3)*B_tan_n1.^2 +Kh_pol1(4)*B_tan_n1.^3 +Kh_pol1(5)*B_tan_n1.^4 +Kh_pol1(6)*B_tan_n1.^5).*(B_tan_n1).^2)*(f_h1)';

        CAL2_loss1 = EC_loss_dis_rad1+EC_loss_dis_tan1+HY_loss_dis_rad_main+HY_loss_dis_tan_main;

    end


    % second frequency range
    n2 = n(index_1+1:end);
    f_h2 = f*n2; % harmonc frequency components

    Kh_pol2 = [0.028101397733866 0.027933569345486 -0.129962460975250 0.182176020474115 -0.116246753245310 0.029377444778697]; % fifth order polynomial
    Ke_pol2 = [3.570394037072378e-05 -2.985395210619579e-05 1.934492016089640e-06 5.224345043919332e-05 -5.474156578122447e-05 1.788517108580745e-05];
    % Edddy currents loss with direct analytical integration




    % Iron losses contributions with the  linear stator core discretization
    % method
    % Eddy currents contribution
    B_rad_n2 = CORR*abs(2*K_B_n(index_1+1:end)./((R_s/R_se).^(2*n2)-1).*((r_dis'/R_se).^(n2-1).*(R_m/R_se).^(n2+1)-(R_m./r_dis)'.^(n2+1)));
    B_tan_n2 = CORR*abs(2*K_B_n(index_1+1:end)./((R_s/R_se).^(2*n2)-1).*((r_dis'/R_se).^(n2-1).*(R_m/R_se).^(n2+1)+(R_m./r_dis)'.^(n2+1)));

    EC_loss_dis_rad2 = dis_weight*((Ke_pol2(1) +Ke_pol2(2)*B_rad_n2 +Ke_pol2(3)*B_rad_n2.^2 +Ke_pol2(4)*B_rad_n2.^3+Ke_pol2(5)*B_rad_n2.^4+Ke_pol2(6)*B_rad_n2.^5).*abs(B_rad_n2).^2)*(f_h2.^2)';
    EC_loss_dis_tan2 = dis_weight*((Ke_pol2(1) +Ke_pol2(2)*B_tan_n2 +Ke_pol2(3)*B_tan_n2.^2 +Ke_pol2(4)*B_tan_n2.^3+Ke_pol2(5)*B_tan_n2.^4+Ke_pol2(6)*B_tan_n2.^5).*abs(B_tan_n2).^2)*(f_h2.^2)';

    % Hysteresis contribution (individual loop per harmonic)
    HY_loss_dis_rad2 = dis_weight*((Kh_pol2(1) +Kh_pol2(2)*B_rad_n2 +Kh_pol2(3)*B_rad_n2.^2 +Kh_pol2(4)*B_rad_n2.^3+Kh_pol2(5)*B_rad_n2.^4+Kh_pol2(6)*B_rad_n2.^5).*abs(B_rad_n2).^2)*(f_h2)';
    HY_loss_dis_tan2 = dis_weight*((Kh_pol2(1) +Kh_pol2(2)*B_tan_n2 +Kh_pol2(3)*B_tan_n2.^2 +Kh_pol2(4)*B_tan_n2.^3+Kh_pol2(5)*B_tan_n2.^4+Kh_pol2(6)*B_tan_n2.^5).*abs(B_tan_n2).^2)*(f_h2)';

    CAL2_loss2 = EC_loss_dis_rad2+EC_loss_dis_tan2+HY_loss_dis_rad2+HY_loss_dis_tan2;


    CAL2_loss = STACK*(CAL2_loss1 + CAL2_loss2);


   



