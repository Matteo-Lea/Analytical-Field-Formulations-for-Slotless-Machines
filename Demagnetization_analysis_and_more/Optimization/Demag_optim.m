%% ANALYTICAL FIELD SOLUTION FROM PM (SCRIPT)

% Matteo Leandro
% matteo.leandro@ntnu.no
% VERSION: 24Oct2022

% This code is used to check the demagnetization phenomenon in optimizad
% slotless motors stored in RESULTS folder. The optimization was run with a
% specific magneet grade. Do not change the PM material properties to make
% hav a consistent final result.
%% Machine parameters
clearvars
clc
close all

mu_0 = 4*pi*1e-7; % air permeability

polepairs = 'p=17'; % available [3 5 7 9 11 13 15 17 19]
mode = 'NO_DEMAG'; % 'NO_DEMAG --> optimization without demagnetization
                   % 'W_DEMAG --> optimization with demagnetization
load(['RESULTS\' polepairs '\' mode '\' polepairs ])

% ex-opt-1
R_s = REP.pos(3); % stator radius [m]
l_a = R_s*1.2; % active lenght [m]
R_wi = R_s; % winding radius [m]
R_w = R_s-REP.pos(4); % winding radius [m]
R_m = R_w-REP.pos(5); % magnets array outer radius [m]
R_r = R_m-REP.pos(6); % outer rotor radius (inner magnets radius) [m]
R_se = R_s + 0.5*(R_m-R_r); % outer stator radius [m]
R_i = 0; % Iron boundary radius facing the magnets
R_ie = R_i;
g = R_w-R_m; % air-gap thickness [m]
R_1 = R_m +g/2; % mid-air-gap radius [m]
R_2 = R_s/2+R_se/2; % mid-stator radius
p = REP.pos(9);
alpha_p = REP.pos(1); % mid-magnet to pole ratio [-]
%% test-geometry
% l_a = 0.2673; % active lenght [m]
% R_s = 0.1447; % stator radius [m]
% R_wi = R_s; % winding radius [m]
% R_w = R_s-0.0096; % winding radius [m]
% R_m = R_w-0.0015; % magnets array outer radius [m]
% R_r = R_m-0.0163; % outer rotor radius (inner magnets radius) [m]
% R_se = 1.01*R_s; % outer stator radius [m]
% R_i = R_r; % Iron boundary radius facing the magnets
% R_ie = R_i-sign(R_w-R_m)*0.1*R_i;
% g = R_w-R_m; % air-gap thickness [m]
% R_1 = R_m +g/2; % mid-air-gap radius [m]
% R_2 = R_s/2+R_se/2; % mid-stator radius
% p = 19;
% alpha_p = 0.5321; % mid-magnet to pole ratio [-]

%% PM material properties used in the optimization
B_r = 1.3756; % remanent flux density [T] (mid-magnet)
B_rs = 1.3756; % remanent flux density [T] (side-magnets)
mu_r = 1.04; % PM recoil permeability [-]
mu_rc = 7.969; % PM coercive permeability (from knee-point to Hc)
B_knee = 0.39429;
%% ------------------------------------------------------------------------

% machine/magnets configuration parameters
Halbach = 0; % set 0 for even segments and 1 for odd segmens arrays 
Halbach_1 = 2; % set number of segments per pole

if mod(Halbach_1,2)~=0 && Halbach==0
    commandwindow
    error([' -)You have set an odd number of segments per pole', ...
          ' but \n selected an even-segment Halbach array.  ', ... 
          ' \n  change either of the following values:', ...
          ' \n Halbach or Halbach_1'])
elseif mod(Halbach_1,2)==0 && Halbach==1
    commandwindow
    error([' -)"You have set an even number of segments per pole \n',...
           ' but selected an odd-segment Halbach array.  \n',...
           ' change either of the following values:\n',...
           ' Halbach or Halbach_1"']);
end


if Halbach_1 == 2 && alpha_p ==0
    alpha_p1 =1;
elseif  Halbach_1 == 2 
    alpha_p1 = alpha_p; % side-magnets + mid-magnet to pole ratio [-]
elseif Halbach_1 == 3
    alpha_p1 = 1; 
else
    alpha_p1 = 0.7;
end
    

theta_m_side = 0; % orientation angle side-magnet [deg]
theta_m_end = 0; % orientation angle end-side-magnet [deg]
theta_m_side = theta_m_side*pi/180;
theta_m_end = theta_m_end*pi/180;
delta_m = (alpha_p+alpha_p1)*pi/(4*p); % mid-point angle of the side-magnet
delta_m1 = (alpha_p1+1)*pi/(4*p); % mid-point angle of the end-side-magnet


% B_r = 1.35; % remanent flux density [T] (mid-magnet)
% B_rs = B_r; % remanent flux density [T] (side-magnets)
% mu_r = 1.0; % PM recoil permeability [-]
% alpha_p = 0.5; % mid-magnet to pole ratio [-]
% alpha_p1 = 2-alpha_p; % side-magnets + mid-magnet to pole ratio [-]  
% p =3; % pole pairs

% % outrunner example
% R_r = (198/2)*1e-3; % outer rotor radius (inner magnets radius) [m]
% R_m = (187/2)*1e-3; % magnets array outer radius [m]
% R_w = (185/2)*1e-3; % winding radius [m]
% R_wi = (182/2)*1e-3; % winding radius [m]
% R_s = (180/2)*1e-3; % stator radius [m]
% R_se = (177/2)*1e-3; % outer stator radius [m]
% R_i = Inf;
% g = abs(R_w-R_m); % air-gap thickness [m]
% R_1 = R_m -g/2; % mid-air-gap radius [m]
% R_2 = R_s/2+R_se/2; % mid-stator radius
% % machine/magnets configuration parameters
% Halbach = 0; % set 0 for even segments and 1 for odd segmens arrays 
% Halbach_1 = 2; % set number of segments per pole
% 
% if mod(Halbach_1,2)~=0 && Halbach==0
%     commandwindow
%     error(['\n -)"You have set an odd number of segments per pole '...
%            '\n but selected an even-segment Halbach array.  '...
%            '\n change either of the following values:'...
%            '\n Halbach or Halbach_1"'])
% elseif mod(Halbach_1,2)==0 && Halbach==1
%     commandwindow
%     error(['\n -)"You have set an even number of segments per pole '...
%            '\n but selected an odd-segment Halbach array.  '...
%            '\n change either of the following values:'...
%            '\n Halbach or Halbach_1"'])
% end
% B_r = 1.4; % remanent flux density [T] (mid-magnet)
% B_rs = B_r; % remanent flux density [T] (side-magnets)
% mu_r = 1.0; % PM recoil permeability [-]
% alpha_p = 0; % mid-magnet to pole ratio [-]
% if Halbach_1 == 2
% alpha_p1 = alpha_p; % side-magnets + mid-magnet to pole ratio [-]
% elseif Halbach_1 == 3
%     alpha_p1 = 1; 
% elseif Halbach_1 == 2 && alpha_p ==0
%     alpha_p1 =1;
% end
% 
% if Halbach_1 == 2 && alpha_p ==0
%     alpha_p1 =1;   
% end
% theta_m_side = 45; % orientation angle side-magnet [deg]
% theta_m_end = 0; % orientation angle end-side-magnet [deg]
% % theta_m_end = 30; % orientation angle end-side-magnet [deg]
% theta_m_side = theta_m_side*pi/180;
% theta_m_end = theta_m_end*pi/180;
% p =26; % pole pairs
% delta_m = (alpha_p+alpha_p1)*pi/(4*p); % mid-point angle of the side-magnet
% delta_m1 = (alpha_p1+1)*pi/(4*p); % mid-point angle of the end-side-magnet



%% radial discretization         
if R_s>R_m %  inrunner
r = linspace(R_m,R_s,10)'; % radius at which the flux density is to be evaluated
r_m = linspace(R_r,R_m,55)'; % radius at which the flux density is to be evaluated
% r_m = r_m(2:end-1);
r_s = linspace(R_s,R_se,10)';
elseif R_s<R_m %  outrunner
r = linspace(R_s,R_m,10)'; % radius at which the flux density is to be evaluated
r_m = linspace(R_m,R_r,55)'; % radius at which the flux density is to be evaluated
r_s = linspace(R_se,R_s,10)';
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
m_PM = 700; % number of harmonics or series components tionfor the magnetization functions
x = 0:1:m_PM;
% x = linspace(0,40,m_PM+1);
n = p*(2*x+1); % existing harmonic series components (n=1,3,5,...)


sigma = (sin(pi*n./n(end))./(pi*n./n(end))).^1; % Lanczos sigma for Gibbs phenomenon reduction
% sigma = 1; % Lanczos sigma for Gibbs phenomenon reduction
%% circumferential discretization
% tic
m_th = 1000; % points along the circumference/arc
theta = linspace(0,2*pi/(2*p),m_th)-pi/(2*p); % circumferential discretization (leave the -pi/(2*p))
Theta = repmat(theta,m_PM+1,1);



%% Magnetization distribution

% PARALLEL MAGNETIZED PM
% Coefficients for defining the magnetization distribution
if Halbach ==0
% Even segments
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

else
% Mid-magnets coefficients
M_1_0 = sin((n-1)*alpha_p*pi/(2*p))./((n-1)*alpha_p*pi/(2*p));
M_2_0 = sin((n+1)*alpha_p*pi/(2*p))./((n+1)*alpha_p*pi/(2*p));
% Side-magnets coefficients
M_1_2 = cos((n-1)*alpha_p*pi/(2*p)+theta_m_side+delta_m)./((n-1)*alpha_p*pi/(2*p));
M_2_2 = cos((n+1)*alpha_p*pi/(2*p)-theta_m_side-delta_m)./((n+1)*alpha_p*pi/(2*p));
M_3_2 = cos((n-1)*alpha_p1*pi/(2*p)+theta_m_side+delta_m)./((n-1)*alpha_p1*pi/(2*p));
M_4_2 = cos((n+1)*alpha_p1*pi/(2*p)-theta_m_side-delta_m)./((n+1)*alpha_p1*pi/(2*p));
% End-side-magnets coefficients
M_1_3 = cos((n-1)*alpha_p1*pi/(2*p)+pi/(2*p))./((n-1)*alpha_p1*pi/(2*p));
M_2_3 = cos((n+1)*alpha_p1*pi/(2*p)-pi/(2*p))./((n+1)*alpha_p1*pi/(2*p));

end
    

if p==1
   if Halbach_1 == 2 && alpha_p ==0
   M_1_2(1) = 0;
%    M_1_2(1) = -sin(theta_m_side+delta_m);
   M_3_2(1) = -sin(theta_m_side+delta_m);
   else
   M_1_0(1) = 1;
   M_1_3(1) = sin(pi/2*p)*(1/alpha_p1-1); 
   M_1_2(1) = -sin(theta_m_side+delta_m);
   M_3_2(1) = -sin(theta_m_side+delta_m);
   end
end

% Magnetization harmonics amplitude
if R_w>R_m % inrunner case
    if Halbach_1 == 2 && alpha_p ==0
        M_r_n_par_mid = B_r/mu_0*alpha_p*(M_1_0+M_2_0);
        M_theta_n_par_mid = -B_r/mu_0*alpha_p*(M_1_0-M_2_0);
        M_r_n_par_end_side = B_rs/mu_0*alpha_p1*(M_1_3-M_2_3);
        M_theta_n_par_end_side = -B_rs/mu_0*alpha_p1*(M_1_3+M_2_3);
        M_r_n_par_side = B_rs/mu_0*((M_1_2-M_2_2)+alpha_p1*(M_4_2-M_3_2));
        M_theta_n_par_side = -B_rs/mu_0*((M_2_2+M_1_2)-alpha_p1*(M_3_2+M_4_2));
    else
        M_r_n_par_mid = B_r/mu_0*alpha_p*(M_1_0+M_2_0);
        M_theta_n_par_mid = -B_r/mu_0*alpha_p*(M_1_0-M_2_0);
        M_r_n_par_end_side = B_rs/mu_0*alpha_p1*(M_1_3-M_2_3);
        M_theta_n_par_end_side = -B_rs/mu_0*alpha_p1*(M_1_3+M_2_3);
        M_r_n_par_side = B_rs/mu_0*(alpha_p*(M_1_2-M_2_2)+alpha_p1*(M_4_2-M_3_2));
        M_theta_n_par_side = -B_rs/mu_0*(alpha_p*(M_2_2+M_1_2)-alpha_p1*(M_3_2+M_4_2));
    end
elseif R_w<R_m % outrunner case
    if p ==1 && Halbach_1 == 2 && alpha_p ==0
        M_r_n_par_mid = B_r/mu_0*alpha_p*(M_1_0+M_2_0);
        M_theta_n_par_mid = -B_r/mu_0*alpha_p*(M_1_0-M_2_0);
        M_r_n_par_end_side = -B_rs/mu_0*alpha_p1*(M_1_3-M_2_3);
        M_theta_n_par_end_side = B_rs/mu_0*alpha_p1*(M_1_3+M_2_3);
        M_r_n_par_side = -B_rs/mu_0*((M_1_2-M_2_2)+alpha_p1*(M_4_2-M_3_2));
        M_theta_n_par_side = B_rs/mu_0*((M_2_2+M_1_2)-alpha_p1*(M_3_2+M_4_2));
    elseif p == 1 && Halbach_1 == 2 && alpha_p ~=0
        M_r_n_par_mid = B_r/mu_0*alpha_p*(M_1_0+M_2_0);
        M_theta_n_par_mid = -B_r/mu_0*alpha_p*(M_1_0-M_2_0);
        M_r_n_par_end_side = -B_rs/mu_0*alpha_p1*(M_1_3-M_2_3);
        M_theta_n_par_end_side = B_rs/mu_0*alpha_p1*(M_1_3+M_2_3);
        M_r_n_par_side = -B_rs/mu_0*(alpha_p*(M_1_2-M_2_2)+alpha_p1*(M_4_2-M_3_2));
        M_theta_n_par_side = B_rs/mu_0*(alpha_p*(M_2_2+M_1_2)-alpha_p1*(M_3_2+M_4_2));
    elseif Halbach_1 == 2 && alpha_p ==0
        M_r_n_par_mid = B_r/mu_0*alpha_p*(M_1_0+M_2_0);
        M_theta_n_par_mid = B_r/mu_0*alpha_p*(M_1_0-M_2_0);
        M_r_n_par_end_side = B_rs/mu_0*alpha_p1*(M_1_3-M_2_3);
        M_theta_n_par_end_side = B_rs/mu_0*alpha_p1*(M_1_3+M_2_3);
        M_r_n_par_side = B_rs/mu_0*((M_1_2-M_2_2)+alpha_p1*(M_4_2-M_3_2));
        M_theta_n_par_side = B_rs/mu_0*((M_2_2+M_1_2)-alpha_p1*(M_3_2+M_4_2));
    elseif Halbach_1 == 2 && alpha_p ~=0
        M_r_n_par_mid = B_r/mu_0*alpha_p*(M_1_0+M_2_0);
        M_theta_n_par_mid = B_r/mu_0*alpha_p*(M_1_0-M_2_0);
        M_r_n_par_end_side = B_rs/mu_0*alpha_p1*(M_1_3-M_2_3);
        M_theta_n_par_end_side = B_rs/mu_0*alpha_p1*(M_1_3+M_2_3);
        M_r_n_par_side = B_rs/mu_0*(alpha_p*(M_1_2-M_2_2)+alpha_p1*(M_4_2-M_3_2));
        M_theta_n_par_side = B_rs/mu_0*(alpha_p*(M_2_2+M_1_2)-alpha_p1*(M_3_2+M_4_2));
    end
end

elseif Halbach ==1
% Odd-segments
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

if p==1
   M_1_0(1) = 1; 
   M_1_2(1) = -sin(theta_m_side+delta_m);
   M_3_2(1) = -sin(theta_m_side+delta_m);
   M_1_3(1) = -sin(theta_m_end+delta_m1);
   M_3_3(1) = -sin(theta_m_end+delta_m1);
   
end

% Magnetization harmonics amplitude
if R_w>R_m % inrunner case
    M_r_n_par_mid = B_r/mu_0*alpha_p*(M_1_0+M_2_0);
    M_theta_n_par_mid = -B_r/mu_0*alpha_p*(M_1_0-M_2_0);
    M_r_n_par_end_side = B_rs/mu_0*(alpha_p1*(M_1_3-M_2_3)+(M_4_3-M_3_3));
    M_theta_n_par_end_side = -B_rs/mu_0*(alpha_p1*(M_2_3+M_1_3)-(M_3_3+M_4_3));
    M_r_n_par_side = B_rs/mu_0*(alpha_p*(M_1_2-M_2_2)+alpha_p1*(M_4_2-M_3_2));
    M_theta_n_par_side = -B_rs/mu_0*(alpha_p*(M_2_2+M_1_2)-alpha_p1*(M_3_2+M_4_2));
elseif p == 1 && R_w<R_m % outrunner case
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
M_theta_par_mid = M_theta_n_par_mid*sin(n'.*Theta);
M_r_par_mid = M_r_n_par_mid*cos(n'.*Theta);


M_r_par_end_side = M_r_n_par_end_side*cos(n'.*Theta);
M_theta_par_end_side = M_theta_n_par_end_side*sin(n'.*Theta);

M_r_par_side = M_r_n_par_side*cos(n'.*Theta);
M_theta_par_side = M_theta_n_par_side*sin(n'.*Theta);

M_r_par = (M_r_n_par_mid+M_r_n_par_end_side+M_r_n_par_side)*cos(n'.*Theta);
M_theta_par = (M_theta_n_par_mid+M_theta_n_par_end_side+M_theta_n_par_side)*sin(n'.*Theta);

    
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

%% WINDING-AIRGAP REGION
if R_s>R_m % INRUNNER winding/air-gap region field compuatation
% FLUX DENSITY IN THE WINDINGS/AIR-GAP REGIONS

% radial component
Amp_r_m = (K_Bn.*((r./R_s).^(n-1).*(R_m/R_s).^(n+1)+(R_m./r).^(n+1)));
B_g_r_m = sigma.*Amp_r_m*cos(n'.*Theta);

Amp_Az = ((K_Bn./(n)*R_m).*((r./R_s).^(n).*(R_m/R_s).^(n)+(R_m./r).^(n)));
A_z = sigma.*Amp_Az*sin(n'.*Theta);

% circumferential component
Amp_theta_m = (K_Bn.*(-(r./R_s).^(n-1).*(R_m/R_s).^(n+1)+(R_m./r).^(n+1)));
B_g_theta_m = sigma.*Amp_theta_m*sin(n'.*Theta);

else % OUTRUNNER winding/air-gap region field compuatation
% FLUX DENSITY IN THE WINDINGS/AIR-GAP REGIONS

% radial component
Amp_r_m = K_Bn_out.*((r./R_m).^(n-1)+(R_s/R_m).^(n-1).*(R_s./r).^(n+1));
B_g_r_m = sigma.*Amp_r_m*cos(n'.*Theta);

% circumferential component
Amp_theta_m = K_Bn_out.*(-(r./R_m).^(n-1)+(R_s/R_m).^(n-1).*(R_s./r).^(n+1));
B_g_theta_m = sigma.*Amp_theta_m*sin(n'.*Theta);

Amp_Az = ((K_Bn_out./(n)*R_m).*((r./R_m).^(n)+(R_s/R_m).^(n).*(R_s./r).^(n)));
A_z = sigma.*Amp_Az*sin(n'.*Theta);

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

Az_m = sigma.*(Amp_Az_m)*sin(n'.*Theta);

B_m_r_m = sigma.*(Amp_m_r_m)*cos(n'.*Theta);

B_m_theta_m = -sigma.*(Amp_m_theta_m)*sin(n'.*Theta);

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
        A_z_1_pl = ((A_zm_1+mu_0.*M_theta_n(1)).*(-R_m.^2.*R_r.^2.*(1+mu_r)-R_m.^4.*(mu_r-1)-R_s.^2.*(R_r.^2.*(mu_r-1)+R_m.^2.*(mu_r-1)))-R_m.^2.*A_zm_1.*(-R_m.^2.*(mu_r-1)+R_s.^2.*(mu_r+1)).*(mu_r-1).*log(R_m)+R_r.^2.*A_zm_1.*(R_s.^2.*(mu_r-1)-R_m.^2.*(mu_r+1)).*(mu_r+1).*log(R_r))./DEN_O(1);
        A_z_1_mi = -(R_m^2*R_r.^2.*(2.*-R_m.^2.*(A_zm_1+mu_0.*M_theta_n(1)).*mu_r-A_zm_1.*(-R_m.^2.*(mu_r-1)+R_s.^2.*(mu_r+1)).*(mu_r+1).*log(R_m/R_r)))/DEN_O(1);
    end
    Amp_Az_m(:,1) = (A_z_1_pl.*r_m+A_z_1_mi./r_m+A_zm_1.*r_m.*log(r_m));
    Amp_m_r_m(:,1) = (A_z_1_pl+A_z_1_mi./r_m.^2+A_zm_1.*log(r_m));
    Amp_m_theta_m(:,1) = (A_z_1_pl-A_z_1_mi./r_m.^2+A_zm_1.*(log(r_m)+1));
end

Az_m = sigma.*(Amp_Az_m)*sin(n'.*Theta);

B_m_r_m = sigma.*(Amp_m_r_m)*cos(n'.*Theta);

B_m_theta_m = -sigma.*(Amp_m_theta_m)*sin(n'.*Theta);


end
% B_m_r_m = sum(B_m_r_m,1);
% B_m_theta_m = sum(B_m_theta_m,1);

%% STATOR CORE REGION
if R_s>R_m % INRUNNER stator core region field computation
    
    Amp_AzS_m = (R_m*G_Sn./n.*((r_s./R_se).^(n).*(R_m/R_se).^(n)-(R_m./r_s).^(n)));
    Az_S_m = sigma.*Amp_AzS_m*sin(n'.*Theta);
    
    Amp_Sr_m = (G_Sn.*((r_s./R_se).^(n-1).*(R_m/R_se).^(n+1)-(R_m./r_s).^(n+1)));
    B_Sr_m = sigma.*Amp_Sr_m*cos(n'.*Theta);

    Amp_Stheta_m = -(G_Sn.*((r_s./R_se).^(n-1).*(R_m/R_se).^(n+1)+(R_m./r_s).^(n+1)));
    B_Stheta_m = sigma.*Amp_Stheta_m*sin(n'.*Theta);
    
else % OUTRUNNER stator core region field computation
    
    Amp_AzS_m = (R_m*G_Sn_out./n.*((r_s./R_m).^(n)-(R_se/R_m).^(n).*(R_se./r_s).^(n)));
    Az_S_m = sigma.*Amp_AzS_m*sin(n'.*Theta);
    
    Amp_Sr_m = (G_Sn_out.*((r_s./R_m).^(n-1)-(R_se/R_m).^(n-1).*(R_se./r_s).^(n+1)));
    B_Sr_m = sigma.*Amp_Sr_m*cos(n'.*Theta);

    Amp_Stheta_m = -(G_Sn_out.*((r_s./R_m).^(n-1)+(R_se/R_m).^(n-1).*(R_se./r_s).^(n+1)));
    B_Stheta_m = sigma.*Amp_Stheta_m*sin(n'.*Theta);
end

%% BACK-MAGNETS REGION
if R_s>R_m % INRUNNER back-magnets region field computation 
    if R_i == R_r % iron backing

        I_in = (A_m_pl.*(R_r/R_m).^(n-1)+A_m_mi+n.*A_zm_n)./(1-(R_ie/R_r).^(2*n));
        
        if p == 1
            I_in(1) = (A_z_1_pl+A_z_1_mi/R_r^2+A_zm_1*log(R_r))./(1-(R_ie/R_r).^2);
        end

        Amp_AzI_m = R_r*I_in./n.*((r_ext./R_r).^(n)-(R_ie/R_r).^(n).*(R_ie./r_ext).^(n));
        Az_I_m = sigma.*Amp_AzI_m*sin(n'.*Theta);

        Amp_Ir_m = I_in.*((r_ext./R_r).^(n-1)-(R_ie/R_r).^(n-1).*(R_ie./r_ext).^(n+1));
        B_Ir_m = sigma.*Amp_Ir_m*cos(n'.*Theta);

        Amp_Itheta_m = -I_in.*((r_ext./R_r).^(n-1)+(R_ie/R_r).^(n-1).*(R_ie./r_ext).^(n+1));
        B_Itheta_m = sigma.*Amp_Itheta_m*sin(n'.*Theta);
    else % no iron backing

        I_in = (((A_zm_n+n.*A_zm_n+mu_0*M_theta_n).*((mu_r+1)-(R_m/R_s).^(2*n)*(mu_r-1))+(A_zm_n-n.*A_zm_n+mu_0*M_theta_n).*(R_r/R_m).^(2*n).*(-(mu_r-1)+(R_m/R_s).^(2*n)*(mu_r+1))-2*(R_r/R_m).^(n-1).*((A_zm_n+mu_0*M_theta_n+mu_r*n.*A_zm_n)+(R_m/R_s).^(2*n).*(A_zm_n+mu_0*M_theta_n-mu_r*n.*A_zm_n)))./DEN_I);

        if p == 1
            I_in(1) = R_r^2*((R_m^2-R_r^2)*(A_zm_1+mu_0*M_theta_n(1))*(R_m^2*(mu_r+1)-R_s^2*(mu_r-1))+2*R_m^2*A_zm_1*(-R_m^2*(mu_r-1)+R_s^2*(mu_r+1))*log(R_m/R_r))/DEN_I(1);
        end
        
        Amp_AzI_m = R_r*I_in./n.*((r_ext./R_r).^(n)+(R_i/R_r).^(n).*(R_i./r_ext).^(n));
        Az_I_m = sigma.*Amp_AzI_m*sin(n'.*Theta);

        Amp_Ir_m = I_in.*((r_ext./R_r).^(n-1)+(R_i/R_r).^(n-1).*(R_i./r_ext).^(n+1));
        B_Ir_m = sigma.*Amp_Ir_m*cos(n'.*Theta);

        Amp_Itheta_m = -I_in.*((r_ext./R_r).^(n-1)-(R_i/R_r).^(n-1).*(R_i./r_ext).^(n+1));
        B_Itheta_m = sigma.*Amp_Itheta_m*sin(n'.*Theta);
    end


else %  OUTRUNNER back-magnets region field computation
    if R_i == R_r % iron backing

        I_out = (A_m_pl+A_m_mi.*(R_m/R_r).^(n+1)+n.*A_zm_n)./(1-(R_r/R_ie).^(2*n));
        
        if p == 1
            I_out(1) = (A_z_1_pl+A_z_1_mi./R_r.^2+A_zm_1*log(R_r))./(1-(R_r/R_ie).^2);
        end

        Amp_AzI_m = R_r*I_out./n.*(-(R_r/R_ie).^(n).*(r_ext./R_ie).^(n)+(R_r./r_ext).^(n));
        Az_I_m = sigma.*Amp_AzI_m*sin(n'.*Theta);

        Amp_Ir_m = I_out.*(-(R_r/R_ie).^(n+1).*(r_ext./R_ie).^(n-1)+(R_r./r_ext).^(n+1));
        B_Ir_m = sigma.*Amp_Ir_m*cos(n'.*Theta);

        Amp_Itheta_m = I_out.*((R_r/R_ie).^(n+1).*(r_ext./R_ie).^(n-1)+(R_r./r_ext).^(n+1));
        B_Itheta_m = sigma.*Amp_Itheta_m*sin(n'.*Theta);
    else % no iron backing
        I_out = (((A_zm_n+n.*A_zm_n+mu_0*M_theta_n).*(R_m/R_r).^(2*n).*((R_s/R_m).^(2*n)*(mu_r+1)-(mu_r-1))+(A_zm_n-n.*A_zm_n+mu_0*M_theta_n).*((mu_r+1)-(R_s/R_m).^(2*n)*(mu_r-1))-2*(R_m/R_r).^(n+1).*((R_s/R_m).^(2*n).*(A_zm_n+mu_0*M_theta_n+mu_r*n.*A_zm_n)+(A_zm_n+mu_0*M_theta_n-mu_r*n.*A_zm_n)))./DEN_O);

        if p == 1
            I_out(1) = ((R_m^2-R_r^2)*(A_zm_1+mu_0*M_theta_n(1))*(R_m^2*(mu_r+1)-R_s^2*(mu_r-1))+2*R_m^2*A_zm_1*(-R_m^2*(mu_r-1)+R_s^2*(mu_r+1))*log(R_m/R_r))/DEN_O(1);
        end
        
        Amp_AzI_m = R_r*I_out./n.*(-(R_r/R_i).^(n).*(r_ext./R_i).^(n)+(R_r./r_ext).^(n));
        Az_I_m = sigma.*Amp_AzI_m*sin(n'.*Theta);

        Amp_Ir_m = I_out.*(-(R_r/R_i).^(n+1).*(r_ext./R_i).^(n-1)+(R_r./r_ext).^(n+1));
        B_Ir_m = sigma.*Amp_Ir_m*cos(n'.*Theta);

        Amp_Itheta_m = -I_out.*((R_r/R_i).^(n+1).*(r_ext./R_i).^(n-1)-(R_r./r_ext).^(n+1));
        B_Itheta_m = sigma.*Amp_Itheta_m*sin(n'.*Theta);
    end
end



%% Demagnetization analysis

if alpha_p~=0
    
    Vol_mag = 0;
    Vol_dem = 0;
    Vol_mag1 = alpha_p*pi*abs(R_m^2-R_r^2)/(2*p);
    dis = 400;
    theta_1 = linspace(-alpha_p*pi/(2*p),alpha_p*pi/(2*p),dis);
%     theta_1 = theta_1(2:end-1);
    Theta_1mag = repmat(theta_1,size(r_m,2),1);
    B_r1 = sigma.*Amp_m_r_m*cos(n'.*theta_1);
    B_theta1 = -sigma.*(Amp_m_theta_m)*sin(n'.*theta_1);
    B_mag_1 = B_r1.*cos(Theta_1mag)+B_theta1.*sin(Theta_1mag);
    B_dem_1 = B_knee - B_mag_1;
%     B_dem_1 = B_knee - B_mag_1(2:end-1,2:end-1);
    B_dem_1(B_dem_1 <= 0) = 0;
    Vol_el_1 = Vol_mag1/numel(B_mag_1);
    
    
    Vol_demag1 = nnz(B_dem_1)*Vol_mag1/numel(B_dem_1);
    Vol_dem = Vol_dem + Vol_demag1;
    Vol_mag = Vol_mag + Vol_mag1;
%     Weighted_vol_demag_perc_1 = nnz(B_dem_1)*Vol_mag1/numel(B_dem_1)*sum(sum(B_dem_1./B_r*mu_r/mu_rc))*100; 
%     Weighted_vol_demag_perc_1 = sum(sum(B_dem_1 > 0.*B_dem_1./B_r*mu_r/mu_rc))/sum(sum(B_dem_1./B_r*mu_r/mu_rc))/numel(B_dem_1)*100;
    Weight = B_dem_1./B_knee*Vol_el_1;
%     Weighted_vol_demag_perc_1 = sum(sum((B_dem_1 > 0).*(B_dem_1*mu_r/mu_rc)./B_r))/numel(B_dem_1)*100;
    Weighted_vol_demag_perc_1 = sum(sum(Weight))/numel(B_dem_1)*100;
    
    [r_m1,t_m1] = meshgrid(r_m,theta_1);
    x_m1 = r_m1.*sin(t_m1);
    y_m1 = r_m1.*cos(t_m1); 
    
    if  mod(Halbach_1,2)==0
        
        Vol_mag2 = 1/2*(1-alpha_p1)*pi*abs(R_m^2-R_r^2)/(2*p);
        Vol_mag3 = Vol_mag2;
        theta_2 = linspace(-pi/(2*p),-alpha_p1*pi/(2*p),dis);
        theta_2 = theta_2(2:end-1);
        Theta_2mag = repmat(theta_2,size(r_m,2),1);
%         theta_3 = linspace(alpha_p1*pi/(2*p),pi/(p)-alpha_p1*pi/(2*p),1000);
        theta_3 = linspace(alpha_p1*pi/(2*p),pi/(2*p),dis);
        theta_3 = theta_3(2:end-1);
        Theta_3mag = repmat(theta_3,size(r_m,2),1);
        B_r2 = sigma.*Amp_m_r_m*cos(n'.*theta_2);
        B_theta2 = sign(R_m-R_w)*sigma.*(Amp_m_theta_m)*sin(n'.*theta_2);
        B_r3 = sigma.*Amp_m_r_m*cos(n'.*theta_3);
        B_theta3 = sign(R_m-R_w)*sigma.*(Amp_m_theta_m)*sin(n'.*theta_3);

        B_mag_2 = B_r2.*cos(Theta_2mag+pi/(2*p)+pi/2)+B_theta2.*sin(Theta_2mag+pi/(2*p)+pi/2);
        B_mag_3 = B_r3.*cos(Theta_3mag-pi/(2*p)-pi/2)+B_theta3.*sin(Theta_3mag-pi/(2*p)-pi/2);
        
        B_dem_2 = B_knee - B_mag_2;
%         B_dem_2 = B_knee - B_mag_2(2:end-1,2:end-1);
        B_dem_2(B_dem_2 <= 0) = 0;
        B_dem_3 = B_knee - B_mag_3;
%         B_dem_3 = B_knee - B_mag_3(2:end-1,2:end-1);
        B_dem_3(B_dem_3 <= 0) = 0;
        
        Vol_demag2 = nnz(B_dem_2)*Vol_mag2/numel(B_dem_2);
        Vol_demag3 = nnz(B_dem_3)*Vol_mag3/numel(B_dem_3);
        Vol_dem = Vol_dem + Vol_demag2 + Vol_demag3;
        Vol_mag = Vol_mag + Vol_mag2 + Vol_mag3;

        Vol_el_2 = Vol_mag2/numel(B_mag_2);
%         Weighted_vol_demag_perc_2 = nnz(B_dem_2)/numel(B_dem_2)*sum(sum(B_dem_2./B_r*mu_r/mu_rc))*100;
%         Weighted_vol_demag_perc_3 = nnz(B_dem_3)/numel(B_dem_3)*sum(sum(B_dem_3./B_r*mu_r/mu_rc))*100;
        Weighted_vol_demag_perc_2 = sum(sum((B_dem_2 > 0).*(B_dem_2*mu_r/mu_rc)./B_r))/numel(B_dem_2)*100;
        Weighted_vol_demag_perc_3 = sum(sum((B_dem_3 > 0).*(B_dem_3*mu_r/mu_rc)./B_r))/numel(B_dem_3)*100;
        
        Weight = [Weight B_dem_2./B_knee*Vol_el_2 B_dem_3./B_knee*Vol_el_2];
        
        [r_m2,t_m2] = meshgrid(r_m,theta_2);
        [r_m3,t_m3] = meshgrid(r_m,theta_3);
        x_m2 = r_m2.*sin(t_m2);
        y_m2 = r_m2.*cos(t_m2);
        x_m3 = r_m3.*sin(t_m3);
        y_m3 = r_m3.*cos(t_m3);
    end
    
    if Halbach_1>=3
        
    
        Vol_mag4 = 1/2*(alpha_p1-alpha_p)*pi*abs(R_m^2-R_r^2)/(2*p);
        Vol_mag5 = Vol_mag4;
        theta_4 = linspace(-alpha_p1*pi/(2*p),-alpha_p*pi/(2*p),1000);
        Theta_4 = repmat(theta_4,m_PM+1,1);
        Theta_4mag = repmat(theta_4,size(r_m,2),1);
        theta_5 = linspace(alpha_p*pi/(2*p),alpha_p1*pi/(2*p),1000);
        Theta_5= repmat(theta_5,m_PM+1,1);
        Theta_5mag = repmat(theta_5,size(r_m,2),1);
        B_r4 = sigma.*Amp_m_r_m*cos(n'.*Theta_4);
        B_theta4 = -sigma.*(Amp_m_theta_m)*sin(n'.*Theta_4);
        B_r5 = sigma.*Amp_m_r_m*cos(n'.*Theta_5);
        B_theta5 = -sigma.*(Amp_m_theta_m)*sin(n'.*Theta_5);

        B_mag_4 = B_r4.*cos(Theta_4mag+pi/(2*p)+pi/2-theta_m_side-delta_m)+B_theta4.*sin(Theta_4mag+pi/(2*p)+pi/2-theta_m_side-delta_m);
        B_mag_5 = B_r5.*cos(Theta_5mag-pi/(2*p)-pi/2+theta_m_side+delta_m)+B_theta5.*sin(Theta_5mag-pi/(2*p)-pi/2+theta_m_side+delta_m);
        
        B_dem_4 = B_knee - B_mag_4;
        B_dem_4(B_dem_4 <= 0) = 0;
        B_dem_5 = B_knee - B_mag_5;
        B_dem_5(B_dem_5 <= 0) = 0;
        
        Vol_demag4 = nnz(B_dem_4)*Vol_mag4/numel(B_dem_4);
        Vol_demag5 = nnz(B_dem_5)*Vol_mag5/numel(B_dem_5);
        Vol_dem = Vol_dem + Vol_demag4 + Vol_demag5;
        Vol_mag = Vol_mag + Vol_mag4 + Vol_mag5;
        
        Weight = [Weight B_dem_4./B_knee B_dem_5./B_knee];
        
        [r_m4,t_m4] = meshgrid(r_m,theta_4);
        [r_m5,t_m5] = meshgrid(r_m,theta_5);
        x_m4 = r_m4.*sin(t_m4);
        y_m4 = r_m4.*cos(t_m4);
        x_m5 = r_m5.*sin(t_m5);
        y_m5 = r_m5.*cos(t_m5);

    end
    
    if Halbach_1==5
        
        Vol_mag6 = 1/2*(1-alpha_p1)*pi*abs(R_m^2-R_r^2)/(2*p);
        Vol_mag7 = Vol_mag6;
        theta_6 = linspace(-pi/(2*p),-alpha_p1*pi/(2*p),1000);
        Theta_6 = repmat(theta_6,m_PM+1,1);
        Theta_6mag = repmat(theta_6,size(r_m,2),1);
        theta_7 = linspace(alpha_p1*pi/(2*p),pi/(2*p),1000);
        Theta_7 = repmat(theta_7,m_PM+1,1);
        Theta_7mag = repmat(theta_7,size(r_m,2),1);
        B_r6 = sigma.*Amp_m_r_m*cos(n'.*Theta_6);
        B_theta6 = -sigma.*(Amp_m_theta_m)*sin(n'.*Theta_6);
        B_r7 = sigma.*Amp_m_r_m*cos(n'.*Theta_7);
        B_theta7 = -sigma.*(Amp_m_theta_m)*sin(n'.*Theta_7);

        B_mag_6 = B_r6.*sin(Theta_6mag+theta_m_end+delta_m1)+B_theta6.*cos(Theta_6mag+theta_m_end+delta_m1);
        B_mag_7 = -B_r7.*sin(Theta_7mag-theta_m_end-delta_m1)-B_theta7.*cos(Theta_7mag-theta_m_end-delta_m1);
        
        B_dem_6 = B_knee - B_mag_6;
        B_dem_6(B_dem_6 <= 0) = 0;
        B_dem_7 = B_knee - B_mag_7;
        B_dem_7(B_dem_7 <= 0) = 0;
        
        Vol_demag6 = nnz(B_dem_6)*Vol_mag6/numel(B_dem_6);
        Vol_demag7 = nnz(B_dem_7)*Vol_mag7/numel(B_dem_7);
        Vol_dem = Vol_dem + Vol_demag6 + Vol_demag7;
        Vol_mag = Vol_mag + Vol_mag6 + Vol_mag7;
        
        Weight = [Weight B_dem_6./B_knee B_dem_7./B_knee];
        
        
        [r_m6,t_m6] = meshgrid(r_m,theta_6);
        [r_m7,t_m7] = meshgrid(r_m,theta_7);
        x_m6 = r_m6.*sin(t_m6);
        y_m6 = r_m6.*cos(t_m6);
        x_m7 = r_m7.*sin(t_m7);
        y_m7 = r_m7.*cos(t_m7);
        
    end
    
    
Perc_demag = Vol_dem/Vol_mag*100;
Weighted_Perc_demag = sum(sum(Weight))/Vol_mag*100;
        
        
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

Levels = linspace(min([min(min(Az_m)) min(min(A_z)) min(min(Az_S_m)) min(min(Az_I_m))]),max([max(max(Az_m)) max(max(A_z)) max(max(Az_S_m)) max(max(Az_I_m))]),10);


figure('Renderer','Painters');
hold on;
if alpha_p~=0
contourf(x_m1,y_m1,B_mag_1',linspace(-4,4,100),'LineColor','none');
end
if  mod(Halbach_1,2)==0
contourf(x_m2,y_m2,B_mag_2',linspace(-4,4,100),'LineColor','none');
contourf(x_m3,y_m3,B_mag_3',linspace(-4,4,100),'LineColor','none');
end
if Halbach_1>=3
contourf(x_m4,y_m4,B_mag_4',linspace(-4,4,100),'LineColor','none');
contourf(x_m5,y_m5,B_mag_5',linspace(-4,4,100),'LineColor','none');
end
if Halbach_1==5
contourf(x_m6,y_m6,B_mag_6',linspace(-4,4,100),'LineColor','none');
contourf(x_m7,y_m7,B_mag_7',linspace(-4,4,100),'LineColor','none');
end
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
plot(R_s*sin(linspace(theta(1),theta(end),100)),R_s*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot(R_w*sin(linspace(theta(1),theta(end),100)),R_w*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot(R_wi*sin(linspace(theta(1),theta(end),100)),R_wi*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot([R_wi*sin(theta(1)+pi/(6*p):pi/(3*p):theta(end)-pi/(6*p));R_w*sin(theta(1)+pi/(6*p):pi/(3*p):theta(end)-pi/(6*p))],[R_wi*cos(theta(1)+pi/(6*p):pi/(3*p):theta(end)-pi/(6*p));R_w*cos(theta(1)+pi/(6*p):pi/(3*p):theta(end)-pi/(6*p))],'linewidth',0.8,'color','k');
set(gca,'visible','off');
colormap(jet)
c = colorbar;
c.Label.String = 'Flux density norm [T]';
caxis([0 1])
axis off
axis image
title('Resulting flux density map along the magnetization direction')


figure('Renderer','Painters');
hold on;
if alpha_p~=0
contourf(x_m1,y_m1,B_dem_1',linspace(0,B_knee,100),'LineColor','none');
end
if  mod(Halbach_1,2)==0
contourf(x_m2,y_m2,B_dem_2',linspace(0,B_knee,100),'LineColor','none');
contourf(x_m3,y_m3,B_dem_3',linspace(0,B_knee,100),'LineColor','none');
end
if Halbach_1>=3
contourf(x_m4,y_m4,B_dem_4',linspace(0,B_knee,100),'LineColor','none');
contourf(x_m5,y_m5,B_dem_5',linspace(0,B_knee,100),'LineColor','none');
end
if Halbach_1==5
contourf(x_m6,y_m6,B_dem_6',linspace(0,B_knee,100),'LineColor','none');
contourf(x_m7,y_m7,B_dem_7',linspace(0,B_knee,100),'LineColor','none');
end
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
plot(R_s*sin(linspace(theta(1),theta(end),100)),R_s*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot(R_w*sin(linspace(theta(1),theta(end),100)),R_w*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot(R_wi*sin(linspace(theta(1),theta(end),100)),R_wi*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot([R_wi*sin(theta(1)+pi/(6*p):pi/(3*p):theta(end)-pi/(6*p));R_w*sin(theta(1)+pi/(6*p):pi/(3*p):theta(end)-pi/(6*p))],[R_wi*cos(theta(1)+pi/(6*p):pi/(3*p):theta(end)-pi/(6*p));R_w*cos(theta(1)+pi/(6*p):pi/(3*p):theta(end)-pi/(6*p))],'linewidth',0.8,'color','k');
set(gca,'visible','off');
colormap(flipud(gray))
% c = colorbar;
% c.Label.String = 'Flux density norm [T]';
caxis([0 B_knee])
axis off
axis image
title(['Demagnetization map (',num2str(Perc_demag),'% of the PM volume)'])
end

if alpha_p == 0
    
    Vol_mag1 = pi/2*abs(R_m^2-R_r^2)/(2*p);
    Vol_mag2 = Vol_mag1; 
    theta_1 = linspace(-pi/(2*p),0,1000);
    theta_1 = theta_1(1:end-1);
    Theta_1 = repmat(theta_1,m_PM+1,1);
    Theta_1mag = repmat(theta_1,size(r_m,2),1);
    theta_2 = linspace(0,pi/(2*p),1000);
    theta_2 = theta_2(1:end-1);
    Theta_2 = repmat(theta_2,m_PM+1,1);
    Theta_2mag = repmat(theta_2,size(r_m,2),1);
    B_r1 = sigma.*Amp_m_r_m*cos(n'.*Theta_1);
    B_theta1 = -sigma.*(Amp_m_theta_m)*sin(n'.*Theta_1);
    B_r2 = sigma.*Amp_m_r_m*cos(n'.*Theta_2);
    B_theta2 = -sigma.*(Amp_m_theta_m)*sin(n'.*Theta_2);
    
    B_mag_1 = B_r1.*sin(Theta_1mag+theta_m_side+delta_m)+B_theta1.*cos(Theta_1mag+theta_m_side+delta_m);
    B_mag_2 = -B_r2.*sin(Theta_2mag-theta_m_side-delta_m)-B_theta2.*cos(Theta_2mag-theta_m_side-delta_m);
    
    B_dem_1 = B_knee - B_mag_1;
    B_dem_1(B_dem_1 <= 0) = 0;
    
    B_dem_2 = B_knee - B_mag_2;
    B_dem_2(B_dem_2 <= 0) = 0;
    
    Vol_demag1 = nnz(B_dem_1)*Vol_mag1/numel(B_dem_1);
    Vol_demag2 = nnz(B_dem_2)*Vol_mag2/numel(B_dem_2);
    Vol_dem =  Vol_demag1 + Vol_demag2;
    Vol_mag = Vol_mag1 + Vol_mag2;
    Perc_demag = (nnz(B_dem_1)+nnz(B_dem_2))/(numel(B_dem_1)+numel(B_dem_1))*100;
    
    % Create polar data
        
[r,t] = meshgrid(r,theta);
[r_m1,t_m1] = meshgrid(r_m,theta_1);
[r_m2,t_m2] = meshgrid(r_m,theta_2);
[r_m,t_m] = meshgrid(r_m,theta);
[r_s,t_s] = meshgrid(r_s,theta);
[r_ext,t_ext] = meshgrid(r_ext,theta);
x = r.*sin(t);
y = r.*cos(t);
h = polar(x,y);
x_m = r_m.*sin(t_m);
y_m = r_m.*cos(t_m);
x_m1 = r_m1.*sin(t_m1);
y_m1 = r_m1.*cos(t_m1);
x_m2 = r_m2.*sin(t_m2);
y_m2 = r_m2.*cos(t_m2);
x_s = r_s.*sin(t_s);
y_s = r_s.*cos(t_s);
x_ext = r_ext.*sin(t_ext);
y_ext = r_ext.*cos(t_ext);
Levels = linspace(min([min(min(Az_m)) min(min(A_z)) min(min(Az_S_m)) min(min(Az_I_m))]),max([max(max(Az_m)) max(max(A_z)) max(max(Az_S_m)) max(max(Az_I_m))]),30);


figure('Renderer','Painters');
hold on;
contourf(x_m1,y_m1,B_mag_1',linspace(-4,4,1000),'LineColor','none');
contourf(x_m2,y_m2,B_mag_2',linspace(-4,4,1000),'LineColor','none');
contour(x, y,A_z',Levels,'LineColor','k','linewidth',2)  
contour(x_m, y_m,Az_m',Levels,'LineColor','k','linewidth',2) %Create contour plots in polar coordinates onto polar chart
plot(R_r*sin(linspace(theta(1),theta(end),100)),R_r*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot(R_m*sin(linspace(theta(1),theta(end),100)),R_m*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot([R_r*sin(theta(1)+(1-alpha_p)*pi/(2*p):pi/(p):theta(end));R_m*sin(theta(1)+(1-alpha_p)*pi/(2*p):pi/(p):theta(end))],[R_r*cos(theta(1)+(1-alpha_p)*pi/(2*p):pi/(p):theta(end));R_m*cos(theta(1)+(1-alpha_p)*pi/(2*p):pi/(p):theta(end))],'linewidth',0.8,'color','k');
plot([R_r*sin(theta(1)+(1+alpha_p)*pi/(2*p):pi/(p):theta(end));R_m*sin(theta(1)+(1+alpha_p)*pi/(2*p):pi/(p):theta(end))],[R_r*cos(theta(1)+(1+alpha_p)*pi/(2*p):pi/(p):theta(end));R_m*cos(theta(1)+(1+alpha_p)*pi/(2*p):pi/(p):theta(end))],'linewidth',0.8,'color','k');
plot([R_r*sin(theta(1)+(1-alpha_p1)*pi/(2*p):pi/(p):theta(end));R_m*sin(theta(1)+(1-alpha_p1)*pi/(2*p):pi/(p):theta(end))],[R_r*cos(theta(1)+(1-alpha_p1)*pi/(2*p):pi/(p):theta(end));R_m*cos(theta(1)+(1-alpha_p1)*pi/(2*p):pi/(p):theta(end))],'linewidth',0.8,'color','k');
plot([R_r*sin(theta(1)+(1+alpha_p1)*pi/(2*p):pi/(p):theta(end));R_m*sin(theta(1)+(1+alpha_p1)*pi/(2*p):pi/(p):theta(end))],[R_r*cos(theta(1)+(1+alpha_p1)*pi/(2*p):pi/(p):theta(end));R_m*cos(theta(1)+(1+alpha_p1)*pi/(2*p):pi/(p):theta(end))],'linewidth',0.8,'color','k');
plot(R_s*sin(linspace(theta(1),theta(end),100)),R_s*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot(R_w*sin(linspace(theta(1),theta(end),100)),R_w*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot(R_wi*sin(linspace(theta(1),theta(end),100)),R_wi*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot([R_wi*sin(theta(1)+pi/(6*p):pi/(3*p):theta(end)-pi/(6*p));R_w*sin(theta(1)+pi/(6*p):pi/(3*p):theta(end)-pi/(6*p))],[R_wi*cos(theta(1)+pi/(6*p):pi/(3*p):theta(end)-pi/(6*p));R_w*cos(theta(1)+pi/(6*p):pi/(3*p):theta(end)-pi/(6*p))],'linewidth',0.8,'color','k');
set(gca,'visible','off');
colormap(jet)
c = colorbar;
c.Label.String = 'Flux density norm [T]';
caxis([0 1])
axis off
axis image
title('Resulting flux density map along the magnetization direction')

% saveas(gcf,'Magnetization_oriented_field_special_Halbach.pdf','pdf')



figure;
hold on;
contourf(x_m1,y_m1,B_dem_1',linspace(0,B_knee,200),'LineColor','none');
contourf(x_m2,y_m2,B_dem_2',linspace(0,B_knee,200),'LineColor','none');
contour(x, y,A_z',Levels,'LineColor','k','linewidth',2)  
contour(x_m, y_m,Az_m',Levels,'LineColor','k','linewidth',2) %Create contour plots in polar coordinates onto polar chart
plot(R_r*sin(linspace(theta(1),theta(end),100)),R_r*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot(R_m*sin(linspace(theta(1),theta(end),100)),R_m*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot([R_r*sin(theta(1)+(1-alpha_p)*pi/(2*p):pi/(p):theta(end));R_m*sin(theta(1)+(1-alpha_p)*pi/(2*p):pi/(p):theta(end))],[R_r*cos(theta(1)+(1-alpha_p)*pi/(2*p):pi/(p):theta(end));R_m*cos(theta(1)+(1-alpha_p)*pi/(2*p):pi/(p):theta(end))],'linewidth',0.8,'color','k');
plot([R_r*sin(theta(1)+(1+alpha_p)*pi/(2*p):pi/(p):theta(end));R_m*sin(theta(1)+(1+alpha_p)*pi/(2*p):pi/(p):theta(end))],[R_r*cos(theta(1)+(1+alpha_p)*pi/(2*p):pi/(p):theta(end));R_m*cos(theta(1)+(1+alpha_p)*pi/(2*p):pi/(p):theta(end))],'linewidth',0.8,'color','k');
plot([R_r*sin(theta(1)+(1-alpha_p1)*pi/(2*p):pi/(p):theta(end));R_m*sin(theta(1)+(1-alpha_p1)*pi/(2*p):pi/(p):theta(end))],[R_r*cos(theta(1)+(1-alpha_p1)*pi/(2*p):pi/(p):theta(end));R_m*cos(theta(1)+(1-alpha_p1)*pi/(2*p):pi/(p):theta(end))],'linewidth',0.8,'color','k');
plot([R_r*sin(theta(1)+(1+alpha_p1)*pi/(2*p):pi/(p):theta(end));R_m*sin(theta(1)+(1+alpha_p1)*pi/(2*p):pi/(p):theta(end))],[R_r*cos(theta(1)+(1+alpha_p1)*pi/(2*p):pi/(p):theta(end));R_m*cos(theta(1)+(1+alpha_p1)*pi/(2*p):pi/(p):theta(end))],'linewidth',0.8,'color','k');
plot(R_s*sin(linspace(theta(1),theta(end),100)),R_s*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot(R_w*sin(linspace(theta(1),theta(end),100)),R_w*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot(R_wi*sin(linspace(theta(1),theta(end),100)),R_wi*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot([R_wi*sin(theta(1)+pi/(6*p):pi/(3*p):theta(end)-pi/(6*p));R_w*sin(theta(1)+pi/(6*p):pi/(3*p):theta(end)-pi/(6*p))],[R_wi*cos(theta(1)+pi/(6*p):pi/(3*p):theta(end)-pi/(6*p));R_w*cos(theta(1)+pi/(6*p):pi/(3*p):theta(end)-pi/(6*p))],'linewidth',0.8,'color','k');
set(gca,'visible','off');
colormap(flipud(gray))
% c = colorbar;
% c.Label.String = 'Flux density norm [T]';
caxis([0 B_knee])
axis off
axis image
title(['Demagnetization map (',num2str(Perc_demag),'% of the PM volume)'])

% saveas(gcf,'Demagnetization_special_Halbach.pdf','pdf') 



    
end

