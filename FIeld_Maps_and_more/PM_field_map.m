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
clearvars -except runs
clc
close all
mu_0 = 4*pi*1e-7; % air permeability


%% inrunner geometry
Inrunner

%% outrunner geometry
% Outrunner

%% RADIAL DISCRETIZATION 

mapping = "no";
torque = "no";
back_emf = "no";
if R_s>R_m %  inrunner
r = linspace(R_m,R_s,20)'; % radius at which the flux density is to be evaluated
% r = (R_m:0.1*1e-3:R_s)'; % radius at which the flux density is to be evaluated
r_m = linspace(R_r,R_m,20)'; % radius at which the flux density is to be evaluated
r_s = linspace(R_s,R_se,10)';
elseif R_s<R_m %  outrunner
r = linspace(R_s,R_m,20)'; % radius at which the flux density is to be evaluated
r_m = linspace(R_m,R_r,20)'; % radius at which the flux density is to be evaluated
r_s = linspace(R_se,R_s,10)';
end

if R_i == R_r %iron backing
    if R_s >R_m % inrunner
        R_ie = 0.95*R_r;
        r_ext = linspace(R_ie,R_r,20)';
    else % outrunner
        R_ie = 1.1*R_r; 
        r_ext = linspace(R_r,R_ie,20)';
    end
else % no iron backing
    if R_s >R_m % inrunner
        R_ie = R_r-0.25*pi*R_r/p;
        r_ext = linspace(R_ie,R_r,50)';
    else % outrunner
        R_ie = R_r+0.5*pi*R_r/p;
        r_ext = linspace(R_r,R_ie,20)';
    end
end
%% useful indices for series harmonics definition
m_PM = 1000; % number of harmonics or series components tionfor the magnetization functions
x = 0:1:m_PM;
% x = linspace(0,40,m_PM+1);
n = p*(2*x+1); % existing harmonic series components (n=1,3,5,...)

sigma = (sin(pi*n./n(end))./(pi*n./n(end))).^3; % Lanczos sigma for Gibbs phenomenon reduction
%% ANGULAR DISCRETIZATION
sec = 1;                                                                    % number of poles to be modeled
m_th = 500*sec;                                                            % points along the modeled sector
% mechanical angle
theta = linspace(0,sec*pi/(p),m_th)-pi/(2*p)*sec; % circumferential discretization (leave the -pi/(2*p))

% m_th = 2000; % points along the circumference/arc
% theta = linspace(0,2*pi/(2*p),m_th)-pi/(2*p); 
theta_i = p*theta(1)+pi/2*sec; % current displacement angle [rad]

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

G_Sn = 2*K_Bn./((R_s/R_se).^(2*n)-1);

G_Sn_out = 2*K_Bn_out./(1-(R_se/R_s).^(2*n));

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

% radial component
Amp_r_m = (K_Bn.*((r./R_s).^(n-1).*(R_m/R_s).^(n+1)+(R_m./r).^(n+1)));
B_g_r_m = sigma.*Amp_r_m*cos(n'.*theta);

Amp_Az = ((K_Bn./(n)*R_m).*((r./R_s).^(n).*(R_m/R_s).^(n)+(R_m./r).^(n)));
A_z = sigma.*Amp_Az*sin(n'.*theta);

% circumferential component
Amp_theta_m = (K_Bn.*(-(r./R_s).^(n-1).*(R_m/R_s).^(n+1)+(R_m./r).^(n+1)));
B_g_theta_m = sigma.*Amp_theta_m*sin(n'.*theta);

% average maximum flux density in the stator core region
B_bi_A = abs(A_z(end,1)/(R_se-R_s));

Flux_in = dot((K_Bn./(n).*R_m).*2.*(R_m./R_s).^(n),sin((n.*pi./(2*p))),2);

% Back-emf and Torque computation

E_n = -(n_cs*l_a/S_Ph*omega*4*sin(n/p*pi/6).*C.*(K_Bn*R_m^2./(n)));
emf_a = E_n*cos(n'.*theta./sec);
emf_b = E_n*cos(n'.*theta./sec-2*pi/3*(n/p)');
emf_c = E_n*cos(n'.*theta./sec-4*pi/3*(n/p)');

% resulting radial component along the circumference
% emf_a(isnan(emf_a))=0; % higher order harmonics are not represented in Matlab double precision floating point representation (NaN numbers are therefore zeroed out)
% emf_a = sum(emf_a,1);
% emf_b(isnan(emf_b))=0; % higher order harmonics are not represented in Matlab double precision floating point representation (NaN numbers are therefore zeroed out)
% emf_b = sum(emf_b,1);
% emf_c(isnan(emf_c))=0; % higher order harmonics are not represented in Matlab double precision floating point representation (NaN numbers are therefore zeroed out)
% emf_c = sum(emf_c,1);

T = -p/omega*(emf_a*I.*cos(p*theta./sec+theta_i)+emf_b*I.*cos(p*theta./sec-2*pi/3+theta_i)+emf_c*I.*cos(p*theta./sec-4*pi/3+theta_i));

T_avg = -p/omega*3/2*E_n(1)*I;

if contains(back_emf, 'yes')
    figure
    hold on
    plot(theta./sec,emf_a)
    plot(theta./sec,emf_b)
    plot(theta./sec,emf_c)
    title('E_{ph}')
    grid on
end
if contains(torque, 'yes')
    figure
    hold on
    plot(theta./sec,T)
    plot(theta./sec,linspace(T_avg,T_avg,m_th))
    title('T_{em}')
    grid on
end


else % OUTRUNNER winding/air-gap region field compuatation
% FLUX DENSITY IN THE WINDINGS/AIR-GAP REGIONS

% radial component
Amp_r_m = K_Bn_out.*((r./R_m).^(n-1)+(R_s/R_m).^(n-1).*(R_s./r).^(n+1));
B_g_r_m = sigma.*Amp_r_m*cos(n'.*theta);

% circumferential component
Amp_theta_m = K_Bn_out.*(-(r./R_m).^(n-1)+(R_s/R_m).^(n-1).*(R_s./r).^(n+1));
B_g_theta_m = sigma.*Amp_theta_m*sin(n'.*theta);

Amp_Az = ((K_Bn_out./(n)*R_m).*((r./R_m).^(n)+(R_s/R_m).^(n).*(R_s./r).^(n)));
A_z = sigma.*Amp_Az*sin(n'.*theta);

% Back-emf and Torque computation

E_n = -(n_cs*l_a/S_Ph*omega*4*sin(n/p*pi/6).*C_out.*(K_Bn_out*R_m^2./(n)));
emf_a = E_n*cos(n'.*theta./sec);
emf_b = E_n*cos(n'.*theta./sec-2*pi/3*(n/p)');
emf_c = E_n*cos(n'.*theta./sec-4*pi/3*(n/p)');

T = -p/omega*(emf_a*I.*cos(p*theta./sec+theta_i)+emf_b*I.*cos(p*theta./sec-2*pi/3+theta_i)+emf_c*I.*cos(p*theta./sec-4*pi/3+theta_i));

T_avg = -p/omega*3/2*E_n(1)*I;

% E_n_prime = (((R_m.^(n).*(-R_m.^(2.*n+1).*(mu_r-1).*(n.*A_zm_n-mu_0.*M_theta_n-A_zm_n)-R_m.*R_r.^(2.*n).*(mu_r+1).*(n.*A_zm_n+mu_0.*M_theta_n+A_zm_n)+2.*R_m.^(n).*R_r.^(1+n).*(A_zm_n+mu_r.*n.*A_zm_n+mu_0.*M_theta_n)))./(n.*(R_m.^(4.*n).*(mu_r-1)^2-R_m.^(2.*n).*R_r.^(2.*n).*(mu_r+1)^2-R_s.^(2.*n).*R_m.^(2.*n).*(mu_r+1).*(mu_r-1)+R_s.^(2.*n).*R_r.^(2.*n).*(mu_r-1).*(mu_r+1))))./(n+2).*(R_w.^(n+2)-R_ws.^(n+2))+R_s.^(2.*n).*((R_m.^(n).*(-R_m.^(2.*n+1).*(mu_r-1).*(n.*A_zm_n-mu_0.*M_theta_n-A_zm_n)-R_m.*R_r.^(2.*n).*(mu_r+1).*(n.*A_zm_n+mu_0.*M_theta_n+A_zm_n)+2.*R_m.^(n).*R_r.^(1+n).*(A_zm_n+mu_r.*n.*A_zm_n+mu_0.*M_theta_n)))./(n.*(R_m.^(4.*n).*(mu_r-1)^2-R_m.^(2.*n).*R_r.^(2.*n).*(mu_r+1)^2-R_s.^(2.*n).*R_m.^(2.*n).*(mu_r+1).*(mu_r-1)+R_s.^(2.*n).*R_r.^(2.*n).*(mu_r-1).*(mu_r+1))))./(2-n).*(R_w.^(2-n)-R_ws.^(2-n)));
% E_n_prime = -n_cs*l_a/S_Ph*omega*4*E_n_prime.*sin(n/p*pi/6);
% E_n_prime(isnan(E_n_prime))=0; % higher order harmonics are not represented in Matlab double precision floating point representation (NaN numbers are therefore zeroed out)
% E_n_prime(isinf(E_n_prime))=0; % higher order harmonics are not represented in Matlab double precision floating point representation (Inf numbers are therefore zeroed out)
% emf_a_RFF = E_n_prime*cos(n'.*theta);
% emf_b_RFF = E_n_prime*cos(n'.*theta-2*pi/3*(n/p)');
% emf_c_RFF = E_n_prime*cos(n'.*theta-4*pi/3*(n/p)');
% 
% T_RFF = -p/omega*(emf_a_RFF*I.*cos(p*theta+theta_i)+emf_b_RFF*I.*cos(p*theta-2*pi/3+theta_i)+emf_c_RFF*I.*cos(p*theta-4*pi/3+theta_i));

if contains(back_emf, 'yes')
    figure
    hold on
    plot(theta./sec,emf_a)	
    plot(theta./sec,emf_b)
    plot(theta./sec,emf_c)
    xlim([-pi/(2*p) pi/(2*p)]*sec)
    xticks(linspace(-pi/(2*p)*sec,pi/(2*p)*sec,3*sec+1))
    xticklabels({'0','\pi/3','2\pi/3','\pi'})
    title('E_{ph}')
    grid on
end

if contains(torque, 'yes')
figure
hold on
plot(theta./sec,T)
plot(theta./sec,linspace(T_avg,T_avg,m_th))
title('T_{em}')
grid on
end

end

%% MAGNETS REGION
if R_s>R_m % INRUNNER magnets region field compuatation
% FLUX DENSITY IN THE MAGNETS REGION

A_m_pl = (-((mu_r+1)-(R_i/R_r).^(2*n).*(mu_r-1)).*((mu_0*M_theta_n+A_zm_n+mu_r*n.*A_zm_n)+(R_m/R_s).^(2*n).*(mu_0*M_theta_n+A_zm_n-mu_r*n.*A_zm_n))+(R_r/R_m).^(n+1).*((R_m/R_s).^(2*n).*(mu_r+1)-(mu_r-1)).*((mu_0*M_theta_n+A_zm_n-mu_r*n.*A_zm_n)+(R_i/R_r).^(2*n).*(mu_0*M_theta_n+A_zm_n+mu_r*n.*A_zm_n)))./DEN_I;

A_m_mi = (-(R_r/R_m).^(n-1).*(-(mu_r-1)+(R_i/R_r).^(2*n).*(mu_r+1)).*((mu_0*M_theta_n+A_zm_n+mu_r*n.*A_zm_n)+(R_m/R_s).^(2*n).*(mu_0*M_theta_n+A_zm_n-mu_r*n.*A_zm_n))-((R_m/R_s).^(2*n).*(mu_r-1)-(mu_r+1)).*((mu_0*M_theta_n+A_zm_n-mu_r*n.*A_zm_n)+(R_i/R_r).^(2*n).*(mu_0*M_theta_n+A_zm_n+mu_r*n.*A_zm_n)))./DEN_I;


Amp_Az_m = ((R_m*A_m_pl.*(r_m./R_m).^n+R_r*A_m_mi.*(R_r./r_m).^n)./n+A_zm_n.*r_m);

% FLUX DENSITY IN THE MAGNETS REGION
% radial component

Amp_m_r_m = (A_m_pl.*(r_m./R_m).^(n-1)+A_m_mi.*(R_r./r_m).^(n+1)+n.*A_zm_n);

% circumferential component

Amp_m_theta_m = (A_m_pl.*(r_m./R_m).^(n-1)-A_m_mi.*(R_r./r_m).^(n+1)+A_zm_n);


% the case pn=1 is adjusted in the following, as the general solution
% exhibits a singularity in that very case
if p==1
    A_z_1_pl = ((A_zm_1+mu_0.*M_theta_n(1)).*(-R_m.^2.*R_r.^2.*(R_r.^2+R_i.^2).*(1+mu_r)+R_m.^4.*(R_r.^2.*(mu_r+1)-R_i.^2.*(mu_r-1))+R_s.^2.*(R_r.^2.*(R_r.^2+R_i.^2).*(mu_r-1)+R_m.^2.*(R_r.^2.*(mu_r+1)-R_i.^2.*(mu_r-1))))+R_m.^2.*A_zm_1.*(-R_m.^2.*(mu_r-1)+R_s.^2.*(mu_r+1)).*(-R_i.^2.*(mu_r-1)+R_r.^2.*(mu_r+1)).*log(R_m)-R_r.^2.*A_zm_1.*(R_s.^2.*(mu_r-1)-R_m.^2.*(mu_r+1)).*(R_r.^2.*(mu_r-1)-R_i.^2.*(mu_r+1)).*log(R_r))./DEN_I(1);
    A_z_1_mi = (-R_r^2*R_m^2*(2*(R_s^2*R_r^2-R_m^2*R_i^2)*(A_zm_1+mu_0*M_theta_n(1))*mu_r+A_zm_1*(R_s^2*(mu_r+1)-R_m^2*(mu_r-1))*(R_r^2*(mu_r-1)-R_i^2*(mu_r+1))*(log(R_m/R_r))))/DEN_I(1);
    Amp_Az_m(:,1) = (A_z_1_pl.*r_m+A_z_1_mi./r_m+A_zm_1.*r_m.*log(r_m));
    Amp_m_r_m(:,1) = (A_z_1_pl+A_z_1_mi./r_m.^2+A_zm_1.*log(r_m));
    Amp_m_theta_m(:,1) = (A_z_1_pl-A_z_1_mi./r_m.^2+A_zm_1.*(log(r_m)+1));
end

Az_m = sigma.*(Amp_Az_m)*sin(n'.*theta);

B_m_r_m = sigma.*(Amp_m_r_m)*cos(n'.*theta);

B_m_theta_m = -sigma.*(Amp_m_theta_m)*sin(n'.*theta);

if R_i == R_r
   B_sat = 1.8; % iron saturation flux density 
   t_bi = abs(Az_m(1,1)/B_sat); % suggested magnetic sleeve thickness
end

else % OUTRUNNER magnets region field compuatation
% FLUX DENSITY IN THE MAGNETS REGION

A_m_pl = (-(-(mu_r+1)+(R_s/R_m).^(2*n).*(mu_r-1)).*((mu_0*M_theta_n+A_zm_n+mu_r*n.*A_zm_n)+(R_r/R_i).^(2*n).*(mu_0*M_theta_n+A_zm_n-mu_r*n.*A_zm_n))-(R_m/R_r).^(n+1).*((R_r/R_i).^(2*n).*(mu_r+1)-(mu_r-1)).*((mu_0*M_theta_n+A_zm_n-mu_r*n.*A_zm_n)+(R_s/R_m).^(2*n).*(mu_0*M_theta_n+A_zm_n+mu_r*n.*A_zm_n)))./DEN_O;

A_m_mi = ((R_m/R_r).^(n-1).*(-(mu_r-1)+(R_s/R_m).^(2*n).*(mu_r+1)).*((mu_0*M_theta_n+A_zm_n+mu_r*n.*A_zm_n)+(R_r/R_i).^(2*n).*(mu_0*M_theta_n+A_zm_n-mu_r*n.*A_zm_n))+((R_r/R_i).^(2*n).*(mu_r-1)-(mu_r+1)).*((mu_0*M_theta_n+A_zm_n-mu_r*n.*A_zm_n)+(R_s/R_m).^(2*n).*(mu_0*M_theta_n+A_zm_n+mu_r*n.*A_zm_n)))./DEN_O;


Amp_Az_m = ((R_r*A_m_pl.*(r_m./R_r).^n+R_m*A_m_mi.*(R_m./r_m).^n)./n+A_zm_n.*r_m);

% FLUX DENSITY IN THE MAGNETS REGION
% radial component
Amp_m_r_m = (A_m_pl.*(r_m./R_r).^(n-1)+A_m_mi.*(R_m./r_m).^(n+1)+n.*A_zm_n);

% circumferential component
Amp_m_theta_m = (A_m_pl.*(r_m./R_r).^(n-1)-A_m_mi.*(R_m./r_m).^(n+1)+A_zm_n);


% the case pn=1 is adjusted in the following, as the general solution
% exhibits a singularity in that very case
if p==1
%     A_z_1_pl = ((A_zm_1+mu_0.*M_theta_n(1)).*(-R_m.^2.*R_r.^2.*(R_r.^2+R_i.^2).*(1+mu_r)+R_m.^4.*(R_r.^2.*(mu_r+1)-R_i.^2.*(mu_r-1))+R_s.^2.*(R_r.^2.*(R_r.^2+R_i.^2).*(mu_r-1)+R_m.^2.*(R_r.^2.*(mu_r+1)-R_i.^2.*(mu_r-1))))+R_m.^2.*A_zm_1.*(-R_m.^2.*(mu_r-1)+R_s.^2.*(mu_r+1)).*(-R_i.^2.*(mu_r-1)+R_r.^2.*(mu_r+1)).*log(R_m)-R_r.^2.*A_zm_1.*(R_s.^2.*(mu_r-1)-R_m.^2.*(mu_r+1)).*(R_r.^2.*(mu_r-1)-R_i.^2.*(mu_r+1)).*log(R_r))./DEN_I(1);
    A_z_1_pl = ((A_zm_1+mu_0.*M_theta_n(1)).*(-R_m.^2.*R_r.^2.*(R_r.^2+R_i.^2).*(1+mu_r)+R_m.^4.*(R_r.^2.*(mu_r+1)-R_i.^2.*(mu_r-1))+R_s.^2.*(R_r.^2.*(R_r.^2+R_i.^2).*(mu_r-1)+R_m.^2.*(R_r.^2.*(mu_r+1)-R_i.^2.*(mu_r-1))))+R_m.^2.*A_zm_1.*(-R_m.^2.*(mu_r-1)+R_s.^2.*(mu_r+1)).*(-R_i.^2.*(mu_r-1)+R_r.^2.*(mu_r+1)).*log(R_m)-R_r.^2.*A_zm_1.*(R_s.^2.*(mu_r-1)-R_m.^2.*(mu_r+1)).*(R_r.^2.*(mu_r-1)-R_i.^2.*(mu_r+1)).*log(R_r))./DEN_O(1);
    A_z_1_mi = -(R_m^2*R_r.^2.*(2.*(R_s.^2.*R_r.^2-R_m.^2.*R_i.^2).*(A_zm_1+mu_0.*M_theta_n(1)).*mu_r+A_zm_1.*(-R_m.^2.*(mu_r-1)+R_s.^2.*(mu_r+1)).*(R_r.^2.*(mu_r-1)-R_i.^2.*(mu_r+1)).*log(R_m/R_r)))/DEN_O(1);
    if R_i == Inf
        A_z_1_pl = ((A_zm_1+mu_0.*M_theta_n(1)).*(-R_m.^2.*R_r.^2.*(1+mu_r)-R_m.^4.*(mu_r-1)+R_s.^2.*(R_r.^2.*(mu_r-1)-R_m.^2.*(mu_r-1)))-R_m.^2.*A_zm_1.*(-R_m.^2.*(mu_r-1)+R_s.^2.*(mu_r+1)).*(mu_r-1).*log(R_m)+R_r.^2.*A_zm_1.*(R_s.^2.*(mu_r-1)-R_m.^2.*(mu_r+1)).*(mu_r+1).*log(R_r))./DEN_O(1);
        A_z_1_mi = -(R_m^2*R_r.^2.*(2.*-R_m.^2.*(A_zm_1+mu_0.*M_theta_n(1)).*mu_r-A_zm_1.*(-R_m.^2.*(mu_r-1)+R_s.^2.*(mu_r+1)).*(mu_r+1).*log(R_m/R_r)))/DEN_O(1);
    end
    Amp_Az_m(:,1) = (A_z_1_pl.*r_m+A_z_1_mi./r_m+A_zm_1.*r_m.*log(r_m));
    Amp_m_r_m(:,1) = (A_z_1_pl+A_z_1_mi./r_m.^2+A_zm_1.*log(r_m));
    Amp_m_theta_m(:,1) = (A_z_1_pl-A_z_1_mi./r_m.^2+A_zm_1.*(log(r_m)+1));
end

Az_m = sigma.*(Amp_Az_m)*sin(n'.*theta);

B_m_r_m = sigma.*(Amp_m_r_m)*cos(n'.*theta);

B_m_theta_m = -sigma.*(Amp_m_theta_m)*sin(n'.*theta);

if R_i == R_r
   B_sat = 1.8; % iron saturation flux density 
   t_bi = abs(Az_m(1,1)/B_sat); % suggested magnetic sleeve thickness
end

end

%% STATOR CORE REGION
if R_s>R_m % INRUNNER stator core region field computation
    
    Amp_AzS_m = (R_m*G_Sn./n.*((r_s./R_se).^(n).*(R_m/R_se).^(n)-(R_m./r_s).^(n)));
    Az_S_m = sigma.*Amp_AzS_m*sin(n'.*theta);
    
    Amp_Sr_m = (G_Sn.*((r_s./R_se).^(n-1).*(R_m/R_se).^(n+1)-(R_m./r_s).^(n+1)));
    B_Sr_m = sigma.*Amp_Sr_m*cos(n'.*theta);

    Amp_Stheta_m = -(G_Sn.*((r_s./R_se).^(n-1).*(R_m/R_se).^(n+1)+(R_m./r_s).^(n+1)));
    B_Stheta_m = sigma.*Amp_Stheta_m*sin(n'.*theta);
    
    
else % OUTRUNNER stator core region field computation
    
    Amp_AzS_m = (R_m*G_Sn_out./n.*((r_s./R_m).^(n)-(R_se/R_m).^(n).*(R_se./r_s).^(n)));
    Az_S_m = sigma.*Amp_AzS_m*sin(n'.*theta);
    
    Amp_Sr_m = (G_Sn_out.*((r_s./R_m).^(n-1)-(R_se/R_m).^(n-1).*(R_se./r_s).^(n+1)));
    B_Sr_m = sigma.*Amp_Sr_m*cos(n'.*theta);

    Amp_Stheta_m = -(G_Sn_out.*((r_s./R_m).^(n-1)+(R_se/R_m).^(n-1).*(R_se./r_s).^(n+1)));
    B_Stheta_m = sigma.*Amp_Stheta_m*sin(n'.*theta);
end

%% BACK-MAGNETS REGION
if R_s>R_m % INRUNNER back-magnets region field computation 
    if R_i == R_r % iron backing

        I_in = (A_m_pl.*(R_r/R_m).^(n-1)+A_m_mi+n.*A_zm_n)./(1-(R_ie/R_r).^(2*n));
        
        if p == 1
            I_in(1) = (A_z_1_pl+A_z_1_mi/R_r^2+A_zm_1*log(R_r))./(1-(R_ie/R_r).^2);
        end

        Amp_AzI_m = R_r*I_in./n.*((r_ext./R_r).^(n)-(R_ie/R_r).^(n).*(R_ie./r_ext).^(n));
        Az_I_m = sigma.*Amp_AzI_m*sin(n'.*theta);

        Amp_Ir_m = I_in.*((r_ext./R_r).^(n-1)-(R_ie/R_r).^(n-1).*(R_ie./r_ext).^(n+1));
        B_Ir_m = sigma.*Amp_Ir_m*cos(n'.*theta);

        Amp_Itheta_m = -I_in.*((r_ext./R_r).^(n-1)+(R_ie/R_r).^(n-1).*(R_ie./r_ext).^(n+1));
        B_Itheta_m = sigma.*Amp_Itheta_m*sin(n'.*theta);
    else % no iron backing

        I_in = (((A_zm_n+n.*A_zm_n+mu_0*M_theta_n).*((mu_r+1)-(R_m/R_s).^(2*n)*(mu_r-1))+(A_zm_n-n.*A_zm_n+mu_0*M_theta_n).*(R_r/R_m).^(2*n).*(-(mu_r-1)+(R_m/R_s).^(2*n)*(mu_r+1))-2*(R_r/R_m).^(n-1).*((A_zm_n+mu_0*M_theta_n+mu_r*n.*A_zm_n)+(R_m/R_s).^(2*n).*(A_zm_n+mu_0*M_theta_n-mu_r*n.*A_zm_n)))./DEN_I);

        if p == 1
            I_in(1) = R_r^2*((R_m^2-R_r^2)*(A_zm_1+mu_0*M_theta_n(1))*(R_m^2*(mu_r+1)-R_s^2*(mu_r-1))+2*R_m^2*A_zm_1*(-R_m^2*(mu_r-1)+R_s^2*(mu_r+1))*log(R_m/R_r))/DEN_I(1);
        end
        
        Amp_AzI_m = R_r*I_in./n.*((r_ext./R_r).^(n)+(R_i/R_r).^(n).*(R_i./r_ext).^(n));
        Az_I_m = sigma.*Amp_AzI_m*sin(n'.*theta);

        Amp_Ir_m = I_in.*((r_ext./R_r).^(n-1)+(R_i/R_r).^(n-1).*(R_i./r_ext).^(n+1));
        B_Ir_m = sigma.*Amp_Ir_m*cos(n'.*theta);

        Amp_Itheta_m = -I_in.*((r_ext./R_r).^(n-1)-(R_i/R_r).^(n-1).*(R_i./r_ext).^(n+1));
        B_Itheta_m = sigma.*Amp_Itheta_m*sin(n'.*theta);
    end


else %  OUTRUNNER back-magnets region field computation
    if R_i == R_r % iron backing

        I_out = (A_m_pl+A_m_mi.*(R_m/R_r).^(n+1)+n.*A_zm_n)./(1-(R_r/R_ie).^(2*n));
        
        if p == 1
            I_out(1) = (A_z_1_pl+A_z_1_mi./R_r.^2+A_zm_1*log(R_r))./(1-(R_r/R_ie).^2);
        end

        Amp_AzI_m = R_r*I_out./n.*(-(R_r/R_ie).^(n).*(r_ext./R_ie).^(n)+(R_r./r_ext).^(n));
        Az_I_m = sigma.*Amp_AzI_m*sin(n'.*theta);

        Amp_Ir_m = I_out.*(-(R_r/R_ie).^(n+1).*(r_ext./R_ie).^(n-1)+(R_r./r_ext).^(n+1));
        B_Ir_m = sigma.*Amp_Ir_m*cos(n'.*theta);

        Amp_Itheta_m = I_out.*((R_r/R_ie).^(n+1).*(r_ext./R_ie).^(n-1)+(R_r./r_ext).^(n+1));
        B_Itheta_m = sigma.*Amp_Itheta_m*sin(n'.*theta);
    else % no iron backing
        I_out = (((A_zm_n+n.*A_zm_n+mu_0*M_theta_n).*(R_m/R_r).^(2*n).*((R_s/R_m).^(2*n)*(mu_r+1)-(mu_r-1))+(A_zm_n-n.*A_zm_n+mu_0*M_theta_n).*((mu_r+1)-(R_s/R_m).^(2*n)*(mu_r-1))-2*(R_m/R_r).^(n+1).*((R_s/R_m).^(2*n).*(A_zm_n+mu_0*M_theta_n+mu_r*n.*A_zm_n)+(A_zm_n+mu_0*M_theta_n-mu_r*n.*A_zm_n)))./DEN_O);

        if p == 1
            I_out(1) = ((R_m^2-R_r^2)*(A_zm_1+mu_0*M_theta_n(1))*(R_m^2*(mu_r+1)-R_s^2*(mu_r-1))+2*R_m^2*A_zm_1*(-R_m^2*(mu_r-1)+R_s^2*(mu_r+1))*log(R_m/R_r))/DEN_O(1);
        end
        
        Amp_AzI_m = R_r*I_out./n.*(-(R_r/R_i).^(n).*(r_ext./R_i).^(n)+(R_r./r_ext).^(n));
        Az_I_m = sigma.*Amp_AzI_m*sin(n'.*theta);

        Amp_Ir_m = I_out.*(-(R_r/R_i).^(n+1).*(r_ext./R_i).^(n-1)+(R_r./r_ext).^(n+1));
        B_Ir_m = sigma.*Amp_Ir_m*cos(n'.*theta);

        Amp_Itheta_m = -I_out.*((R_r/R_i).^(n+1).*(r_ext./R_i).^(n-1)-(R_r./r_ext).^(n+1));
        B_Itheta_m = sigma.*Amp_Itheta_m*sin(n'.*theta);
    end
end


if contains(mapping, 'yes')
% Create polar data
[r,t] = meshgrid(r,theta);
[r_m,t_m] = meshgrid(r_m,theta);
[r_s,t_s] = meshgrid(r_s,theta);
[r_ext,t_ext] = meshgrid(r_ext,theta);
x = r.*sin(t);
y = r.*cos(t);
% h = polar(x,y);
x_m = r_m.*sin(t_m);
y_m = r_m.*cos(t_m);
x_s = r_s.*sin(t_s);
y_s = r_s.*cos(t_s);
x_ext = r_ext.*sin(t_ext);
y_ext = r_ext.*cos(t_ext);


%%

Norm_Bgm = sqrt(B_g_r_m'.^2+B_g_theta_m'.^2);
Norm_Bm_m = sqrt(B_m_r_m'.^2+B_m_theta_m'.^2);
Norm_BS_m = sqrt(B_Sr_m'.^2+B_Stheta_m'.^2);
Norm_BI_m = sqrt(B_Ir_m'.^2+B_Itheta_m'.^2);
Levels = linspace(min([min(min(Az_m)) min(min(A_z)) min(min(Az_S_m)) min(min(Az_I_m))]),max([max(max(Az_m)) max(max(A_z)) max(max(Az_S_m)) max(max(Az_I_m))]),16);
B_min = 0;
B_max = 1.7;
B_levels = linspace(B_min,B_max,100);


if Halbach ==0
figure;
hold on;
contourf(x,y,Norm_Bgm,B_levels,'edgecolor','none');
contourf(x_m,y_m,Norm_Bm_m,B_levels,'edgecolor','none');
contourf(x_s,y_s,Norm_BS_m,B_levels,'edgecolor','none');
contourf(x_ext,y_ext,Norm_BI_m,B_levels,'edgecolor','none');
contour(x, y,A_z',Levels,'LineColor','k','linewidth',2)
contour(x_m, y_m,Az_m',Levels,'LineColor','k','linewidth',2) %Create contour plots in polar coordinates onto polar chart
contour(x_s, y_s,Az_S_m',Levels,'LineColor','k','linewidth',2) %Create contour plots in polar coordinates onto polar chart
contour(x_ext, y_ext,Az_I_m',Levels,'LineColor','k','linewidth',2)
% Circumferential boundaries
plot(R_ie*sin(linspace(theta(1),theta(end),100)),R_ie*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot(R_r*sin(linspace(theta(1),theta(end),100)),R_r*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot(R_m*sin(linspace(theta(1),theta(end),100)),R_m*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot(R_s*sin(linspace(theta(1),theta(end),100)),R_s*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot(R_se*sin(linspace(theta(1),theta(end),100)),R_se*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot([R_r*sin(theta(1)+(1-alpha_p)*pi/(2*p):pi/(p):theta(end));R_m*sin(theta(1)+(1-alpha_p)*pi/(2*p):pi/(p):theta(end))],[R_r*cos(theta(1)+(1-alpha_p)*pi/(2*p):pi/(p):theta(end));R_m*cos(theta(1)+(1-alpha_p)*pi/(2*p):pi/(p):theta(end))],'linewidth',0.8,'color','k');
plot([R_r*sin(theta(1)+(1+alpha_p)*pi/(2*p):pi/(p):theta(end));R_m*sin(theta(1)+(1+alpha_p)*pi/(2*p):pi/(p):theta(end))],[R_r*cos(theta(1)+(1+alpha_p)*pi/(2*p):pi/(p):theta(end));R_m*cos(theta(1)+(1+alpha_p)*pi/(2*p):pi/(p):theta(end))],'linewidth',0.8,'color','k');
plot([R_r*sin(theta(1)+(1-alpha_p1)*pi/(2*p):pi/(p):theta(end));R_m*sin(theta(1)+(1-alpha_p1)*pi/(2*p):pi/(p):theta(end))],[R_r*cos(theta(1)+(1-alpha_p1)*pi/(2*p):pi/(p):theta(end));R_m*cos(theta(1)+(1-alpha_p1)*pi/(2*p):pi/(p):theta(end))],'linewidth',0.8,'color','k');
plot([R_r*sin(theta(1)+(1+alpha_p1)*pi/(2*p):pi/(p):theta(end));R_m*sin(theta(1)+(1+alpha_p1)*pi/(2*p):pi/(p):theta(end))],[R_r*cos(theta(1)+(1+alpha_p1)*pi/(2*p):pi/(p):theta(end));R_m*cos(theta(1)+(1+alpha_p1)*pi/(2*p):pi/(p):theta(end))],'linewidth',0.8,'color','k');
% Winding slots
plot(R_w*sin(linspace(theta(1),theta(end),100)),R_w*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot(R_ws*sin(linspace(theta(1),theta(end),100)),R_ws*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot([R_ws*sin(theta(1)+pi/(6*p):pi/(3*p):theta(end)-pi/(6*p));R_w*sin(theta(1)+pi/(6*p):pi/(3*p):theta(end)-pi/(6*p))],[R_ws*cos(theta(1)+pi/(6*p):pi/(3*p):theta(end)-pi/(6*p));R_w*cos(theta(1)+pi/(6*p):pi/(3*p):theta(end)-pi/(6*p))],'linewidth',0.8,'color','k');
set(gca,'visible','off');
colormap(jet)
c = colorbar;
c.Label.String = 'Flux density norm [T]';
caxis([B_min B_max])
% Hide the POLAR function data and leave annotations
% set(h,'Visible','off')
% set(h_m,'Visible','off')
% Turn off axes and set square aspect ratio
axis off
axis image
% camroll(-90)
title({'Flux density norm map and isopotential lines',' in the whole 2D domain (optimal formulation)'})


elseif Halbach == 1
    
figure;
hold on;
contourf(x,y,Norm_Bgm,B_levels,'LineStyle','none');
contourf(x_m,y_m,Norm_Bm_m,B_levels,'LineStyle','none');
contourf(x_s,y_s,Norm_BS_m,B_levels,'LineStyle','none');
contourf(x_ext,y_ext,Norm_BI_m,B_levels,'LineStyle','none');
contour(x, y,A_z',Levels,'LineColor','k','linewidth',1.2)
contour(x_m, y_m,Az_m',Levels,'LineColor','k','linewidth',1.2) 
contour(x_s, y_s,Az_S_m',Levels,'LineColor','k','linewidth',1.2) 
contour(x_ext, y_ext,Az_I_m',Levels,'LineColor','k','linewidth',1.2) 
% Circumferential boundaries
plot(R_ie*sin(linspace(theta(1),theta(end),100)),R_ie*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot(R_r*sin(linspace(theta(1),theta(end),100)),R_r*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot(R_m*sin(linspace(theta(1),theta(end),100)),R_m*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot(R_s*sin(linspace(theta(1),theta(end),100)),R_s*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot(R_se*sin(linspace(theta(1),theta(end),100)),R_se*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
% Magnets segments
plot([R_r*sin(theta(1)+(1-alpha_p)*pi/(2*p):pi/(p):theta(end));R_m*sin(theta(1)+(1-alpha_p)*pi/(2*p):pi/(p):theta(end))],[R_r*cos(theta(1)+(1-alpha_p)*pi/(2*p):pi/(p):theta(end));R_m*cos(theta(1)+(1-alpha_p)*pi/(2*p):pi/(p):theta(end))],'linewidth',0.8,'color','k');
plot([R_r*sin(theta(1)+(1+alpha_p)*pi/(2*p):pi/(p):theta(end));R_m*sin(theta(1)+(1+alpha_p)*pi/(2*p):pi/(p):theta(end))],[R_r*cos(theta(1)+(1+alpha_p)*pi/(2*p):pi/(p):theta(end));R_m*cos(theta(1)+(1+alpha_p)*pi/(2*p):pi/(p):theta(end))],'linewidth',0.8,'color','k');
plot([R_r*sin(theta(1)+(1-alpha_p1)*pi/(2*p):pi/(p):theta(end));R_m*sin(theta(1)+(1-alpha_p1)*pi/(2*p):pi/(p):theta(end))],[R_r*cos(theta(1)+(1-alpha_p1)*pi/(2*p):pi/(p):theta(end));R_m*cos(theta(1)+(1-alpha_p1)*pi/(2*p):pi/(p):theta(end))],'linewidth',0.8,'color','k');
plot([R_r*sin(theta(1)+(1+alpha_p1)*pi/(2*p):pi/(p):theta(end));R_m*sin(theta(1)+(1+alpha_p1)*pi/(2*p):pi/(p):theta(end))],[R_r*cos(theta(1)+(1+alpha_p1)*pi/(2*p):pi/(p):theta(end));R_m*cos(theta(1)+(1+alpha_p1)*pi/(2*p):pi/(p):theta(end))],'linewidth',0.8,'color','k');
plot([R_r*sin(theta(1):pi/(p):theta(end));R_m*sin(theta(1):pi/(p):theta(end))],[R_r*cos(theta(1):pi/(p):theta(end));R_m*cos(theta(1):pi/(p):theta(end))],'linewidth',0.8,'color','k');
% Winding slots
plot(R_w*sin(linspace(theta(1),theta(end),100)),R_w*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot(R_ws*sin(linspace(theta(1),theta(end),100)),R_ws*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot([R_ws*sin(theta(1)+pi/(6*p):pi/(3*p):theta(end)-pi/(6*p));R_w*sin(theta(1)+pi/(6*p):pi/(3*p):theta(end)-pi/(6*p))],[R_ws*cos(theta(1)+pi/(6*p):pi/(3*p):theta(end)-pi/(6*p));R_w*cos(theta(1)+pi/(6*p):pi/(3*p):theta(end)-pi/(6*p))],'linewidth',0.8,'color','k');
set(gca,'visible','off');
colormap(jet)
c = colorbar;
c.Label.String = 'Flux density norm [T]';
caxis([B_min B_max])
% Hide the POLAR function data and leave annotations
% set(h,'Visible','off')
% set(h_m,'Visible','off')
% Turn off axes and set square aspect ratio
axis off
axis image
title({'Flux density norm map and isopotential lines',' in the whole 2D domain (optimal formulation)'})
    
end
end

