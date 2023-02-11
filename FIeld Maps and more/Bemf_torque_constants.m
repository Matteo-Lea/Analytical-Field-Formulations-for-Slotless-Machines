%% ANALYTICAL FIELD SOLUTION FROM PM (SCRIPT)

% Matteo Leandro
% matteo.leandro@ntnu.no
% VERSION: 06Nov2020

% This code finds the field solution in the whole domain of slotess
% machines. Inrunner and outrunner cases are treated separately as the OFF
% is implemented. The analytical computation of back-emf and torque is
% performed as well and the results are plotted. Finally the field map is
% plotted. 
% For improved accuracy, the field solution in every region is multiplyied
% by the Lanczos sigma cubed (variable sigma in the code)  to remove the 
% Gibbs phenomenon where it may occur.
% After the geometry is loaded, the problem is conveniently discretized:
% radial discretization, number of harmonics to be considered,  and
% circumferential discretization.

%% Machine parameters

mu_0 = 4*pi*1e-7; % air permeability

%% inrunner geometry
Inrunner

%% outrunner geometry
% Outrunner

%% RADIAL DISCRETIZATION 
if R_s>R_m %  inrunner
r = linspace(R_m,R_s,50)'; % radius at which the flux density is to be evaluated
r_m = linspace(R_r,R_m,50)'; % radius at which the flux density is to be evaluated
r_s = linspace(R_s,R_se,50)';
elseif R_s<R_m %  outrunner
r = linspace(R_s,R_m,50)'; % radius at which the flux density is to be evaluated
r_m = linspace(R_m,R_r,50)'; % radius at which the flux density is to be evaluated
r_s = linspace(R_se,R_s,50)';
end

if R_i == R_r %iron backing
    if R_s >R_m % inrunner
        R_ie = 0.9*R_r;
        r_ext = linspace(R_ie,R_r,50)';
    else % outrunner
        R_ie = 1.1*R_r; 
        r_ext = linspace(R_r,R_ie,50)';
    end
else % no iron backing
    if R_s >R_m % inrunner
        R_ie = R_r-0.5*pi*R_r/p;
        r_ext = linspace(R_ie,R_r,50)';
    else % outrunner
        R_ie = R_r+0.5*pi*R_r/p;
        r_ext = linspace(R_r,R_ie,50)';
    end
end
%% useful indices for series harmonics definition
m_PM = 300; % number of harmonics or series components tionfor the magnetization functions
x = 0:1:m_PM;
% x = linspace(0,40,m_PM+1);
n = p*(2*x+1); % existing harmonic series components (n=1,3,5,...)

sigma = (sin(pi*n./n(end))./(pi*n./n(end))).^3; % Lanczos sigma for Gibbs phenomenon reduction
%% circumferential discretization
sec = 1;                                                                    % number of poles to be modeled
m_th = 100*sec;                                                            % points along the modeled sector
% mechanical angle
theta = linspace(0,sec*pi/(p),m_th)-pi/(2*p)*sec;

% m_th = 2000; % points along the circumference/arc
% theta = linspace(0,2*pi/(2*p),m_th)-pi/(2*p); % circumferential discretization (leave the -pi/(2*p))
theta_i = p*theta(1)+pi/2*sec; % current displacement angle [rad]
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
        M_r_n_par_end_side = -B_rs/mu_0*(alpha_p1*(M_1_3-M_2_3)+(M_4_3-M_3_3));
        M_theta_n_par_end_side = B_rs/mu_0*(alpha_p1*(M_2_3+M_1_3)-(M_3_3+M_4_3));
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


% magnetization reconstruction from harmonics components (not needed)
% M_theta_par_mid = M_theta_n_par_mid*sin(n'.*Theta);
% M_r_par_mid = M_r_n_par_mid*cos(n'.*Theta);
% 
% 
% M_r_par_end_side = M_r_n_par_end_side*cos(n'.*Theta);
% M_theta_par_end_side = M_theta_n_par_end_side*sin(n'.*Theta);
% 
% M_r_par_side = M_r_n_par_side*cos(n'.*Theta);
% M_theta_par_side = M_theta_n_par_side*sin(n'.*Theta);
% 
% M_r_par = (M_r_n_par_mid+M_r_n_par_end_side+M_r_n_par_side)*cos(n'.*Theta);
% M_theta_par = (M_theta_n_par_mid+M_theta_n_par_end_side+M_theta_n_par_side)*sin(n'.*Theta);

    
M_r_n = M_r_n_par_mid+M_r_n_par_end_side+M_r_n_par_side;
M_theta_n = (M_theta_n_par_mid+M_theta_n_par_end_side+M_theta_n_par_side);

%% Field solution

% HARMONICS COEFFICIENTS FROM THE MAGNETIC VECTOR POTENTIAL
% (WINDINGS/AIR-GAP REGIONS)
% coefficient from the particular solution of Poisson's equation accounting
% for the magnetization distribution
A_zm_n = mu_0*(n.*M_r_n+M_theta_n)./(n.^2-1);
A_zm_1 = -mu_0/2*(M_r_n(1)+M_theta_n(1)); % singularity correction when n*p=1


DEN_I = ((mu_r^2-1)*(-(R_i/R_r).^(2*n)-(R_m/R_s).^(2*n)+(R_r/R_s).^(2*n)+(R_i/R_m).^(2*n))+(mu_r+1)^2*(1-(R_i/R_s).^(2*n))+(mu_r-1)^2*((R_m/R_s).^(2*n).*(R_i/R_r).^(2*n)-(R_r/R_m).^(2*n)));

DEN_O = ((mu_r^2-1)*(-(R_s/R_r).^(2*n)-(R_m/R_i).^(2*n)+(R_r/R_i).^(2*n)+(R_s/R_m).^(2*n))+(mu_r+1)^2*((R_s/R_i).^(2*n)-1)+(mu_r-1)^2*((R_m/R_r).^(2*n)-(R_r/R_i).^(2*n).*(R_s/R_m).^(2*n)));

K_Bn = (((-A_zm_n+n.*A_zm_n-mu_0*M_theta_n).*((mu_r+1)-(R_i/R_r).^(2*n)*(mu_r-1))+(A_zm_n+n.*A_zm_n+mu_0*M_theta_n).*((R_r/R_m).^(2*n)*(mu_r-1)-(R_i/R_m).^(2*n)*(mu_r+1))+2*(R_r/R_m).^(n+1).*((A_zm_n+mu_0*M_theta_n-mu_r*n.*A_zm_n)+(R_i/R_r).^(2*n).*(A_zm_n+mu_0*M_theta_n+mu_r*n.*A_zm_n)))./DEN_I);

K_Bn_out = (((-A_zm_n+n.*A_zm_n-mu_0*M_theta_n).*((R_m/R_i).^(2*n)*(mu_r+1)-(R_m/R_r).^(2*n)*(mu_r-1))+(A_zm_n+n.*A_zm_n+mu_0*M_theta_n).*((R_r/R_i).^(2*n)*(mu_r-1)-(mu_r+1))+2*(R_r/R_m)*(R_r/R_i).^(n).*(R_m/R_i).^(n).*(A_zm_n+mu_0*M_theta_n-mu_r*n.*A_zm_n)+2*(R_m/R_r).^(n-1).*(A_zm_n+mu_0*M_theta_n+mu_r*n.*A_zm_n))./DEN_O);

if p == 1
   DEN_I(1) = -((R_s^2+R_m^2)*(R_m^2-R_r^2)*(R_r^2+R_i^2)-2*(R_m^2+R_r^2)*(R_m^2*R_i^2-R_s^2*R_r^2)*mu_r+(R_s^2-R_m^2)*(R_m^2-R_r^2)*(R_r^2-R_i^2)*mu_r^2);
   DEN_O(1) = DEN_I(1);
   K_Bn(1) = R_s^2*((R_m^2-R_r^2)*(A_zm_1+mu_0*M_theta_n(1))*(R_r^2*(mu_r+1)-R_i^2*(mu_r-1))+2*R_r^2*A_zm_1*(R_r^2*(mu_r-1)-R_i^2*(mu_r+1))*log(R_r/R_m))/DEN_I(1);
   K_Bn_out(1) = K_Bn(1)*R_m^2/R_s^2;   
   if R_i == Inf
       DEN_O(1) = -((R_s^2+R_m^2)*(R_m^2-R_r^2)-2*(R_m^2+R_r^2)*R_m^2*mu_r-(R_s^2-R_m^2)*(R_m^2-R_r^2)*mu_r^2);
       K_Bn_out(1) =  -R_m^2*((R_m^2-R_r^2)*(A_zm_1+mu_0*M_theta_n(1))*(mu_r-1)+2*R_r^2*A_zm_1*(mu_r+1)*log(R_r/R_m))/DEN_O(1);
   end   
end

C = (R_m/R_s).^(n).*(R_ws^2*(R_ws/R_s).^(n)-R_w^2*(R_w/R_s).^(n))./((n+2)*R_m)+(R_ws^2*(R_m/R_ws).^(n)-R_w^2*(R_m/R_w).^(n))./((2-n)*R_m);

C_out = (R_w^2*(R_w/R_m).^(n)-R_ws^2*(R_ws/R_m).^(n))./((n+2)*R_m)+(R_w^2/R_s*(R_s/R_w).^(n)-R_ws^2/R_s*(R_s/R_ws).^(n))./(2-n).*(R_s/R_m).^(n+1);

% the back-emf integral exhibits a singularity for np = 2
if p == 2
    C(1) = R_m.*log(R_ws./R_w) + R_m./(4*R_s.^4).*(R_ws.^4-R_w.^4);
    C_out(1) = 1./R_m.^3.*(R_s.^4.*log(R_w./R_ws)+R_w.^4/4-R_ws.^4/4);
end


%% WINDING-AIRGAP REGION
if R_s>R_m % INRUNNER winding/air-gap region field compuatation
% FLUX DENSITY IN THE WINDINGS/AIR-GAP REGIONS

% Back-emf and Torque constants computation

E_n = -(n_cs*l_a/S_Ph*p*omega*4*sin(n/p*pi/6).*C.*(K_Bn*R_m^2./(n)));
emf_a = E_n*cos(n'.*theta./sec);
k_v = max(abs(emf_a)) / omega; % [V/rad/s]
k_t = 3/2*k_v; % [Nm/A_pk]

else % OUTRUNNER winding/air-gap region field compuatation
% FLUX DENSITY IN THE WINDINGS/AIR-GAP REGIONS

% Back-emf and Torque computation

E_n = -(n_cs*l_a/S_Ph*p*omega*4*sin(n/p*pi/6).*C_out.*(K_Bn_out*R_m^2./(n)));
emf_a = E_n*cos(n'.*theta./sec);
k_v = max(abs(emf_a)) / omega; % [V/rad/s]
k_t = 3/2*k_v; % [Nm/A_pk]

end

