%% outrunner example
% l_a = 0.1; % active lenght [m]
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
% B_r = 1.4; % remanent flux density [T] (mid-magnet)
% B_rs = 1.4; % remanent flux density [T] (side-magnets)
% mu_r = 1.0; % PM recoil permeability [-]
% alpha_p = 0.5; % mid-magnet to pole ratio [-]
% alpha_p1 = 0.7; % side-magnets + mid-magnet to pole ratio [-]  
% p = 26; % pole pairs


% X100
l_a =   0.03; % active lenght [m]
R_r = (0.097/2); % outer rotor radius (inner magnets radius) [m]
R_m = (0.0882/2); % magnets array outer radius [m]
R_w = (0.087/2); % winding radius [m]
R_s = (0.0804/2); % stator radius [m]
R_wi = (0.0816/2); % winding radius [m]
% R_s = R_wi; % winding radius [m]
R_se = (0.078/2); % outer stator radius [m]
R_i = Inf; % Iron boundary radius facing the magnets
% R_i = R_r; % Iron boundary radius facing the magnets
% R_ie = R_i-sign(R_w-R_m)*0.1*R_i;
g = R_w-R_m; % air-gap thickness [m]
R_1 = R_m +g/2; % mid-air-gap radius [m]
R_2 = R_s/2+R_se/2; % mid-stator radius
p = 14;

% machine/magnets configuration parameters
B_r = 1.4; % remanent flux density [T] (mid-magnet)
B_rs = B_r; % remanent flux density [T] (side-magnets)
mu_r = 1.05; % PM recoil permeability [-]
alpha_p = 0.5; % mid-magnet to pole ratio [-]
alpha_p1 = 0.7; % side-magnets + mid-magnet to pole ratio [-] 

% machine/magnets configuration parameters
Halbach = 0; % set 0 for even segments and 1 for odd segmens arrays 
Halbach_1 = 2; % set number of segments per pole

if mod(Halbach_1,2)~=0 && Halbach==0
    commandwindow
    error(['\n -)"You have set an odd number of segments per pole '...
           '\n but selected an even-segment Halbach array.  '...
           '\n change either of the following values:'...
           '\n Halbach or Halbach_1"'])
elseif mod(Halbach_1,2)==0 && Halbach==1
    commandwindow
    error(['\n -)"You have set an even number of segments per pole '...
           '\n but selected an odd-segment Halbach array.  '...
           '\n change either of the following values:'...
           '\n Halbach or Halbach_1"'])
end

theta_m_end = 30; % orientation angle end-side-magnet [deg]
theta_m_side = 60; % orientation angle side-magnet [deg]
theta_m_side = theta_m_side*pi/180;
theta_m_end = theta_m_end*pi/180;
if Halbach_1 == 2 && alpha_p ==0
    alpha_p1 = 1;
    delta_m = (alpha_p+alpha_p1)*pi/(4*p); % mid-point angle of the side-magnet
    delta_m1 = (alpha_p1+1)*pi/(4*p); % mid-point angle of the end-side-magnet
    theta_m_side = theta_m_side; % orientation angle side-magnet [deg]
    theta_m_end = 0; % orientation angle end-side-magnet [deg]
elseif Halbach_1 == 2
    alpha_p1 = alpha_p; % side-magnets + mid-magnet to pole ratio [-]
    delta_m = (alpha_p+alpha_p1)*pi/(4*p); % mid-point angle of the side-magnet
    delta_m1 = pi/(2*p); % mid-point angle of the end-side-magnet
    theta_m_side = 0; % orientation angle side-magnet [deg]
    theta_m_end = 0; % orientation angle end-side-magnet [deg] IT MUST BE ZERO!!!
elseif Halbach_1 == 3
    alpha_p1 = 1; 
    delta_m = (alpha_p+alpha_p1)*pi/(4*p); % mid-point angle of the side-magnet
    delta_m1 = (alpha_p1+1)*pi/(4*p); % mid-point angle of the end-side-magnet
    theta_m_side = theta_m_side; % orientation angle side-magnet [deg]
    theta_m_end = theta_m_end; % orientation angle end-side-magnet [deg]
else
    alpha_p1 = alpha_p1;
    delta_m = (alpha_p+alpha_p1)*pi/(4*p); % mid-point angle of the side-magnet
    delta_m1 = (alpha_p1+1)*pi/(4*p); % mid-point angle of the end-side-magnet
    theta_m_side = theta_m_side; % orientation angle side-magnet [deg]
    theta_m_end = theta_m_end; % orientation angle end-side-magnet [deg]
end



S_Ph = abs(pi*(R_wi^2-R_w^2)/(6*p)); % phase belt cross section
N_tc = 6; % number of conductors per coil
q = 1; %coils per pole and per phase
b = 3; % parallel ways per phase
I = 1; % phase current peak value [A]
n_cs = N_tc*q/b;
omega = 1000*2*pi/60*p; % mechanical angular frequency [rad/s]

var = 0; % variable for setting the type of magnetization: cartesian (var=0) 
         % or cylindrical (var=1)