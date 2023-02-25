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

    %% Iron loss analysis based on VARCOext model (frequency domain):

    % Iron loss analysis based on the three component equation with variable coefficients VARCO:
    % c1(B)*B^2*f^2 + c2(f)*B^d(B,f)*f + c3(B)*B^1.5*f^1.5

    % The following part of the code, evaluates iron losses based on the above
    % written equation. The model is based on experimental loss data, which are
    % splitted in two frequency ranges (to get a better fit of the data
    % themselves). For this reason, the maximum frequency of the first range
    % (fr1) has to be set (the model can be extended to more ranges by
    % followingthe same coding fashion).
    % the different coefficients and exponents are defined in the following for
    % the sake of generality
    fr1 = 200; % max frequency of first range
    f = p*n_rpm/60; % fundamental frequency [Hz]
    index_1 = round((fr1/f)/2);

    VARCO_loss1 = 0;
    VARCO_loss2 = 0;

    % first frequency range
    if index_1>0
        n1 = n(1:index_1);
        f_h1 = f*n1./p; % harmonc frequency components

        Ka_pol1 = [-0.000684053327551060,0.00705300085577693,-0.0264530877405776,0.0511894417752783,-0.0529348376233235,0.0280915283079930,-0.00599300797039526];
        Ke_pol1 = [0.000143553970503323,-0.000974107646688384,0.00371002775066561,-0.00711572773726470,0.00720184370162939,-0.00369935784664406,0.000761173577766608]; 
        alpha_01f = 1.77374174459906;
        alpha_11f = -0.449288221156261;
        alpha_21f = 0.507332014433408;
        alpha_31f = 0.103107349390035;
        K_h1 = 0.0149;

        % Iron losses contributions with the  linear stator core discretization
        % method
        % Eddy currents contribution
        B_rad_n1 = CORR*abs(2*K_B_n(1:index_1)./((R_s/R_se).^(2*n1)-1).*((r_dis'/R_se).^(n1-1).*(R_m/R_se).^(n1+1)-(R_m./r_dis)'.^(n1+1)));
        B_tan_n1 = CORR*abs(2*K_B_n(1:index_1)./((R_s/R_se).^(2*n1)-1).*((r_dis'/R_se).^(n1-1).*(R_m/R_se).^(n1+1)+(R_m./r_dis)'.^(n1+1)));

        EC_loss_dis_rad1 = dis_weight*((Ke_pol1(1) +Ke_pol1(2)*B_rad_n1 +Ke_pol1(3)*B_rad_n1.^2 +Ke_pol1(4)*B_rad_n1.^3 +Ke_pol1(5)*B_rad_n1.^4 +Ke_pol1(6)*B_rad_n1.^5 +Ke_pol1(7)*B_rad_n1.^6).*abs(B_rad_n1).^2)*(f_h1.^2)';
        EC_loss_dis_tan1 = dis_weight*((Ke_pol1(1) +Ke_pol1(2)*B_tan_n1 +Ke_pol1(3)*B_tan_n1.^2 +Ke_pol1(4)*B_tan_n1.^3 +Ke_pol1(5)*B_tan_n1.^4 +Ke_pol1(6)*B_tan_n1.^5 +Ke_pol1(7)*B_tan_n1.^6).*abs(B_tan_n1).^2)*(f_h1.^2)';

        % Hysteresis contribution (single hysteresis loop)
        B_max_r = CORR*max((2*K_B_n./((R_s/R_se).^(2*n)-1).*((r_dis'/R_se).^(n-1).*(R_m/R_se).^(n+1)-(R_m./r_dis)'.^(n+1)))*cos(n'.*p.*Theta),[],2);
        B_max_theta = CORR*max(-(2*K_B_n./((R_s/R_se).^(2*n)-1).*((r_dis'/R_se).^(n-1).*(R_m/R_se).^(n+1)+(R_m./r_dis)'.^(n+1)))*sin(n'.*p.*Theta),[],2);
        HY_loss_dis_rad_main = l_a*rho_Fe*K_h1*(pi*((r_dis+delta_r/2).^2-(r_dis-delta_r/2).^2))*(B_max_r).^((alpha_01f)+(alpha_11f).*B_max_r+(alpha_21f).*B_max_r.^2+(alpha_31f).*B_max_r.^3)*f;
        HY_loss_dis_tan_main = l_a*rho_Fe*K_h1*(pi*((r_dis+delta_r/2).^2-(r_dis-delta_r/2).^2))*(B_max_theta).^((alpha_01f)+(alpha_11f).*B_max_theta+(alpha_21f).*B_max_theta.^2+(alpha_31f).*B_max_theta.^3)*f;

        % Hysteresis contribution (individual loop per harmonic)
        HY_loss_dis_rad1 = (l_a*rho_Fe*(pi*((r_dis+delta_r/2).^2-(r_dis-delta_r/2).^2)))*(K_h1.*abs(B_rad_n1).^((alpha_01f)+(alpha_11f).*B_rad_n1+(alpha_21f).*B_rad_n1.^2+(alpha_31f).*B_rad_n1.^3))*(f_h1)';
        HY_loss_dis_tan1 = (l_a*rho_Fe*(pi*((r_dis+delta_r/2).^2-(r_dis-delta_r/2).^2)))*(K_h1.*abs(B_tan_n1).^((alpha_01f)+(alpha_11f).*B_tan_n1+(alpha_21f).*B_tan_n1.^2+(alpha_31f).*B_tan_n1.^3))*(f_h1)';

        % Excess loss contribution
        EX_loss_dis_rad1 = (l_a*rho_Fe*(pi*((r_dis+delta_r/2).^2-(r_dis-delta_r/2).^2)))*((Ka_pol1(1) +Ka_pol1(2)*B_rad_n1 +Ka_pol1(3)*B_rad_n1.^2 +Ka_pol1(4)*B_rad_n1.^3 +Ka_pol1(5)*B_rad_n1.^4 +Ka_pol1(6)*B_rad_n1.^5 +Ka_pol1(7)*B_rad_n1.^6).*abs(B_rad_n1).^1.5)*(f_h1.^1.5)';
        EX_loss_dis_tan1 = (l_a*rho_Fe*(pi*((r_dis+delta_r/2).^2-(r_dis-delta_r/2).^2)))*((Ka_pol1(1) +Ka_pol1(2)*B_tan_n1 +Ka_pol1(3)*B_tan_n1.^2 +Ka_pol1(4)*B_tan_n1.^3 +Ka_pol1(5)*B_tan_n1.^4 +Ka_pol1(6)*B_tan_n1.^5 +Ka_pol1(7)*B_tan_n1.^6).*abs(B_tan_n1).^1.5)*(f_h1.^1.5)';

    %     VARCO_loss1 = EC_loss_dis_rad1+EC_loss_dis_tan1+HY_loss_dis_rad_main+HY_loss_dis_tan_main+EX_loss_dis_rad1+EX_loss_dis_tan1;
        VARCO_loss1 = EC_loss_dis_rad1+EC_loss_dis_tan1+HY_loss_dis_rad1+HY_loss_dis_tan1+EX_loss_dis_rad1+EX_loss_dis_tan1;

    end

    % second frequency range frequency range
        n2 = n(index_1+1:end);
        f_h2 = f*n2./p; % harmonc frequency components

        Ka_pol2 = [0.000361021540203931,-0.00183522385958509,0.0145976593845193,-0.0328413327185457,0.0350554346644131,-0.0185946800279470,0.00394282361816964]; 
        Ke_pol2 = [1.08095426908372e-05,0.000128260820873277,-0.000763756851638181,0.00170339793507077,-0.00182561295888373,0.000952923251909155,-0.000194056477894611]; 
        alpha_02f =2.00403993679833; % quadratic polynomial in f to avoid singluarity for high frequency and low flux density
        alpha_12f =-2.06776899886409;
        alpha_22f =4.18447459806868;
        alpha_32f =-1.60370911854258;
        K_h2 = 0.0120140326304892;

        % Iron losses contributions with the  linear stator core discretization
        % method
        % Eddy currents contribution
        B_rad_n2 = CORR*abs(2*K_B_n(index_1+1:end)./((R_s/R_se).^(2*n2)-1).*((r_dis'/R_se).^(n2-1).*(R_m/R_se).^(n2+1)-(R_m./r_dis)'.^(n2+1)));
        B_tan_n2 = CORR*abs(2*K_B_n(index_1+1:end)./((R_s/R_se).^(2*n2)-1).*((r_dis'/R_se).^(n2-1).*(R_m/R_se).^(n2+1)+(R_m./r_dis)'.^(n2+1)));


        EC_loss_dis_rad2 = dis_weight*((Ke_pol2(1) +Ke_pol2(2)*B_rad_n2 +Ke_pol2(3)*B_rad_n2.^2 +Ke_pol2(4)*B_rad_n2.^3 +Ke_pol2(5)*B_rad_n2.^4 +Ke_pol2(6)*B_rad_n2.^5 +Ke_pol2(7)*B_rad_n2.^6).*abs(B_rad_n2).^2)*(f_h2.^2)';
        EC_loss_dis_tan2 = dis_weight*((Ke_pol2(1) +Ke_pol2(2)*B_tan_n2 +Ke_pol2(3)*B_tan_n2.^2 +Ke_pol2(4)*B_tan_n2.^3 +Ke_pol2(5)*B_tan_n2.^4 +Ke_pol2(6)*B_tan_n2.^5 +Ke_pol2(7)*B_tan_n2.^6).*abs(B_tan_n2).^2)*(f_h2.^2)';

        % Hysteresis contribution (single hysteresis loop)
        B_max_r = CORR*max((2*K_B_n./((R_s/R_se).^(2*n)-1).*((r_dis'/R_se).^(n-1).*(R_m/R_se).^(n+1)-(R_m./r_dis)'.^(n+1)))*cos(n'.*p.*Theta),[],2);
        B_max_theta = CORR*max(-(2*K_B_n./((R_s/R_se).^(2*n)-1).*((r_dis'/R_se).^(n-1).*(R_m/R_se).^(n+1)+(R_m./r_dis)'.^(n+1)))*sin(n'.*p.*Theta),[],2);
        HY_loss_dis_rad_main = dis_weight*(K_h2)*(B_max_r).^((alpha_02f)+(alpha_12f).*B_max_r +(alpha_22f).*B_max_r.^2 +(alpha_32f).*B_max_r.^3)*f;
        HY_loss_dis_tan_main = dis_weight*(K_h2)*(B_max_theta).^((alpha_02f)+(alpha_12f).*B_max_theta+(alpha_22f).*B_max_theta.^2 +(alpha_32f).*B_max_theta.^3)*f;

        % Hysteresis contribution (individual loop per harmonic)
        HY_loss_dis_rad2 = dis_weight*(((K_h2).*abs(B_rad_n2).^((alpha_02f)+(alpha_12f).*B_rad_n2+(alpha_22f).*B_rad_n2.^2 +(alpha_32f).*B_rad_n2.^3))*(f_h2)');
        HY_loss_dis_tan2 = dis_weight*(((K_h2).*abs(B_tan_n2).^((alpha_02f)+(alpha_12f).*B_tan_n2+(alpha_22f).*B_tan_n2.^2 +(alpha_32f).*B_tan_n2.^2))*(f_h2)');

        % Excess loss contribution
        EX_loss_dis_rad2 = dis_weight*((Ka_pol2(1) +Ka_pol2(2)*B_rad_n2 +Ka_pol2(3)*B_rad_n2.^2 +Ka_pol2(4)*B_rad_n2.^3 +Ka_pol2(5)*B_rad_n2.^4 +Ka_pol2(6)*B_rad_n2.^5 +Ka_pol2(7)*B_rad_n2.^6).*abs(B_rad_n2).^1.5)*(f_h2.^1.5)';
        EX_loss_dis_tan2 = dis_weight*((Ka_pol2(1) +Ka_pol2(2)*B_tan_n2 +Ka_pol2(3)*B_tan_n2.^2 +Ka_pol2(4)*B_tan_n2.^3 +Ka_pol2(5)*B_tan_n2.^4 +Ka_pol2(6)*B_tan_n2.^5 +Ka_pol2(7)*B_tan_n2.^6).*abs(B_tan_n2).^1.5)*(f_h2.^1.5)';

    %     VARCO_loss2 = EC_loss_dis_rad2+EC_loss_dis_tan2+EX_loss_dis_rad2+EX_loss_dis_tan2+HY_loss_dis_rad_main+HY_loss_dis_tan_main;
    VARCO_loss2 = EC_loss_dis_rad2+EC_loss_dis_tan2+EX_loss_dis_rad2+EX_loss_dis_tan2+HY_loss_dis_rad2+HY_loss_dis_tan2;


    VARCO_loss = STACK*(VARCO_loss1 + VARCO_loss2);





