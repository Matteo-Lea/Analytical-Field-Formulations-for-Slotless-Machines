%% MACHINE DATA
clearvars;
close all;
clc;
filename = 'Plutus.fem';
% Initial data
p = 16; % pole pairs
f = 1900; % [Hz] rated frequency
omegam_nom = 2*pi*f/p; % rated angular speed
Qs = 3*2*p; % number of slots (total in the stator)
np = 1; % number of poles considered in the model (works with either 1 or 2)
if np>2
    commandwindow
    error([' -)You may consider to set a number of poles in the model', ...
          ' to be equal to either 1 or 2.'])
end
Qsim = np*Qs/(2*p); % number of stator slot in the simulation
Qs_p = Qs/(2*p); % number of slots over a pole
t = gcd(Qs,p); % greatest common divisor of (Q,p)
alpha = np*360/(2*p); % [deg] portion of the machine to be studied (2 poles, 6 slots)
alpha_s = 360/Qs; % slot/tooth angle [deg]
theta_p = 2*180/(2*p); % pole angle [deg]
nmag = 2*2*p; % total number of magnets
q = Qs/(3*2*p); % number of slots per pole and per phase
alpha_se = p*alpha_s; % electrical slot/tooth angle [deg]
mu_0 = 4*pi*1e-7;

% machine/magnets configuration parameters
Halbach = 0; % set 0 for even segments and 1 for odd segmens arrays 
Halbach_1 = 2; % set number of segments per pole
if Halbach == 0
    n_mag = Halbach_1*np+1; % number of magnets (a partial magnet counts as one)
else
    n_mag = Halbach_1*np; % number of magnets (a partial magnet counts as one)
end

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

B_r = 1.3756; % remanent flux density [T] (mid-magnet)
B_rs = B_r; % remanent flux density [T] (side-magnets)
H_c = 760*1e3; % PM hysteresis curve coercivity [A/m] 1 [Oe] = 79.57747 [A/m]
mu_r = 1.04; % PM recoil permeability [-]
mu_rc = 7.969; % PM coercive permeability (from knee-point to Hc)
H_int = (B_r-mu_0*mu_rc*H_c)/(mu_0*(mu_rc-mu_r));
B_int = mu_0*mu_rc*(H_int+H_c);
alpha_p = 0.7175; % mid-magnet to pole ratio [-]
% if Halbach_1 == 2 
% alpha_p1 = alpha_p; % side-magnets + mid-magnet to pole ratio [-]
% elseif Halbach_1 == 3
%     alpha_p1 = 1; 
% elseif Halbach_1 == 2 && alpha_p ==0
%     alpha_p1 =1;
% else
%     alpha_p1 = 0.8;
% end

if Halbach_1 == 2 && alpha_p == 0
    alpha_p1 =1;
elseif  Halbach_1 == 2 
    alpha_p1 = alpha_p; % side-magnets + mid-magnet to pole ratio [-]
elseif Halbach_1 == 3
    alpha_p1 = 1; 
else
    alpha_p1 = 0.7;
end
    

theta_m_side = 45; % orientation angle side-magnet [deg]
theta_m_end = 25; % orientation angle end-side-magnet [deg]
theta_m_side = theta_m_side*pi/180;
theta_m_end = theta_m_end*pi/180;
delta_m = (alpha_p+alpha_p1)*pi/(4*p); % mid-point angle of the side-magnet
delta_m1 = (alpha_p1+1)*pi/(4*p); % mid-point angle of the end-side-magnet


% Define the double air-gap geometry 
% ex-opt 1
% l_a = 0.2764; % active lenght [m]
% R_s = 0.15; % stator radius [m]
% R_wi = R_s; % winding radius [m]
% R_w = R_s-0.0078; % winding radius [m]
% R_m = R_w-0.0015; % magnets array outer radius [m]
% R_r = R_m-0.0087; % outer rotor radius (inner magnets radius) [m]
% R_se = 1.01*R_s; % outer stator radius [m]
% R_i = R_r; % Iron boundary radius facing the magnets
% R_ie = R_i-sign(R_w-R_m)*0.1*R_i;
% g = R_w-R_m; % air-gap thickness [m]
% R_1 = R_m +g/2; % mid-air-gap radius [m]
% R_2 = R_s/2+R_se/2; % mid-stator radius

% ex-opt 2
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

% ex-opt 3
% l_a = 0.2673; % active lenght [m]
% R_s = 0.15; % stator radius [m]
% R_wi = R_s; % winding radius [m]
% R_w = R_s-0.0083; % winding radius [m]
% R_m = R_w-0.0015; % magnets array outer radius [m]
% R_r = R_m-0.0087; % outer rotor radius (inner magnets radius) [m]
% R_se = 1.01*R_s; % outer stator radius [m]
% R_i = R_r; % Iron boundary radius facing the magnets
% R_ie = R_i-sign(R_w-R_m)*0.1*R_i;
% g = R_w-R_m; % air-gap thickness [m]
% R_1 = R_m +g/2; % mid-air-gap radius [m]
% R_2 = R_s/2+R_se/2; % mid-stator radius

% ex-opt 4
l_a = 0.1204; % active lenght [m]
R_s = 0.15; % stator radius [m]
R_wi = R_s; % winding radius [m]
R_w = R_s-0.0082; % winding radius [m]
R_m = R_w-0.0015; % magnets array outer radius [m]
R_r = R_m-0.0089; % outer rotor radius (inner magnets radius) [m]
R_se = 1.01*R_s; % outer stator radius [m]
R_i = R_r; % Iron boundary radius facing the magnets
R_ie = R_i-sign(R_w-R_m)*0.1*R_i;
g = R_w-R_m; % air-gap thickness [m]
R_1 = R_m +g/2; % mid-air-gap radius [m]
R_2 = R_s/2+R_se/2; % mid-stator radius

% l_a =   0.015155218524355094; % active lenght [m]
% R_r = (  0.07792470053917462/2); % outer rotor radius (inner magnets radius) [m]
% R_m = (  0.08356746598268504/2); % magnets array outer radius [m]
% R_w = (  0.08476746598268504/2); % winding radius [m]
% R_s = (  0.08765845176028497/2); % stator radius [m]
% R_wi = (  0.08735845176028498/2); % winding radius [m]
% R_se = ( 0.08968372987341358/2); % outer stator radius [m]
% R_i = R_r; % Iron boundary radius facing the magnets
% R_ie = R_i-sign(R_w-R_m)*0.1*R_i;

gap_w = abs(R_wi-R_s); % [mm] winding air-gap
gap_m = abs(R_m-R_w); % [mm] main air-gap
s_t = abs(R_w-R_wi); % [mm] stator thickness
t_m = abs(R_m-R_r); % [mm] magnet thickness
% tang = 0.4; % pole portion occupied by the tangential magnet 
% rad = 0.6; % pole portion occupied by the radial magnet 
tau_p = pi*2*R_w/(2*p); % stator pole pitch
% lew = 0.5*pi*tau_p/2; % end-winding length
% lew=10.8;

a = tau_p/2;
b = 0.5*a;
h_e = (a-b)^2/(a+b)^2;
lew = pi*(a+b)*(1+1/4*h_e+1/64*h_e^2+1/256*h_e^3)/2;


%% MAGNETS DOMAIN DISCRETIZATION (over one pole)

% Mid-magnet is intended as the magnet radially magnetized, side-magnets as
% the eventual magnet with the angled magnetization, while end-side-magnet
% the one cosing the array (being it tangentially magnetized or angled

dis_ang_mid = 10; % number of angle partitions mid-magnet
if Halbach == 0
    dis_ang_mid = dis_ang_mid*2;
end

alpha_dis_mid = linspace((1-alpha_p)*180/(2*p),(1+alpha_p)*180/(2*p),dis_ang_mid); % vector containing the mid-magnet discretized angle
dis_ang_side = dis_ang_mid/2; % number of angle partitions side-magnet
alpha_dis_side = linspace((1-alpha_p1)*180/(2*p),(1-alpha_p)*180/(2*p),dis_ang_side); % vector containing the side-magnet discretized angle
dis_ang_end = dis_ang_mid/2; % number of angle partitions end-side-magnet
alpha_dis_end = linspace(0,(1-alpha_p1)*180/(2*p),dis_ang_end); % vector containing the end-side-magnet discretized angle
rad_dis = 10; % divisions along the radial direction
r_m_dis = linspace(R_r,R_m,rad_dis); % magnets discretization along the radial direction

% mid-points vectors for domains definition
alpha_span_mid = abs(alpha_dis_mid(2)-alpha_dis_mid(1)); % one division angle
alpha_span_side = abs(alpha_dis_side(2)-alpha_dis_side(1)); % one division angle
alpha_span_end = abs(alpha_dis_end(2)-alpha_dis_end(1)); % one division angle
rad_thick = abs(r_m_dis(2)-r_m_dis(1)); % one division radial thickness





%% SLOT MATRIX
% these three sets of elements must be changed according to the winding
% configuration!
% phase A: 9 elements for each pole
ka_1 = [1 0 0 ]' ; ka_2 = [-1 0 0]';
% phase B: 9 elements for each pole
kb_1 = [0 0 1]' ; kb_2 = [0 0 -1]';
% phase C: 9 elements for each pole
kc_1 = [0 -1 0]'; kc_2 = [0 1 0]';

ka = zeros(Qsim,1);
kb = zeros(Qsim,1);
kc = zeros(Qsim,1);
for ii = 1:np
    if mod(ii,2)==1
        ka((ii-1)*Qs_p+1:ii*Qs_p) = ka_1;
        kb((ii-1)*Qs_p+1:ii*Qs_p) = kb_1;
        kc((ii-1)*Qs_p+1:ii*Qs_p) = kc_1;
    elseif mod(ii,2)==0
        ka((ii-1)*Qs_p+1:ii*Qs_p) = ka_2;
        kb((ii-1)*Qs_p+1:ii*Qs_p) = kb_2;
        kc((ii-1)*Qs_p+1:ii*Qs_p) = kc_2;
    end    
end


%% WINDING PARAMETERS
% J = 6.3; % [A/mm^2] current density
npp = 1; % parallel paths
In = 128; % [A] nominal current (rms value)
E = 40; % [V] maximum available voltage (rms phase voltage)
d_litz = 0.14; % [mm] litz wire diameter
lpc = 150; % litz per conductor
K_fill = 0.43; % fill factor


Slitz = pi*(d_litz*1e-3)^2/4;
nc = 2;
ncs = nc/npp; % equivalent-series conductors per slot
Scond = Slitz*lpc*1e6; % [mm2]conductor section
Sc_eq = npp*Scond; % [mm2]conductor equivalent section
Jcu = (In/npp)/Scond; % [A/mm2] copper current density
Ns = ncs*Qs/3; % equivalent number of series conductor per phase
T_rise = 0; % [°C]temperature rise 
rho = 0.018*(1+0.004*T_rise); % copper resistivity [ohm*mm^2/m]
Rs =rho*(l_a+lew)*1e-3*Ns/Sc_eq;

K_d = sin(q*alpha_se/2*pi/180)/(q*sin(alpha_se/2*pi/180)); % distribution factor
Ls_a = 3/pi*4*pi*1e-7*(K_d*Ns/(2*p))^2*((2*R_w-sign(R_w-R_m)*gap_m))*l_a/(gap_m+gap_w+s_t+t_m)*1e-3;
% Ls_ew = 4*pi*1e-7*tau_p/2*(log(8*tau_p/s_t/2)-1.75)*(nc*q)^2*2*p*1e-3/2;
Ls_ew = 4*pi*1e-7*lew*1e-3*2*p*q^2*ncs^2*0.35;
L_s = Ls_a+Ls_ew;

%% FEMM Analysis

openfemm(1)
newdocument(0)

%problem definition
mi_probdef(0, 'millimeters', 'planar', 1e-8, l_a, 10,0)

%% geometry definition

Az_zero_1 = 0.3*R_ie;
Az_zero_2 = 0.1*R_se;
% magnets back-iron boundaries
mi_drawarc(R_i,0,R_i*cos(alpha*pi/180),R_i*sin(alpha*pi/180),alpha,1)
mi_drawarc(R_ie,0,R_ie*cos(alpha*pi/180),R_ie*sin(alpha*pi/180),alpha,1)
% stator winding boundaries
mi_drawarc(R_w,0,R_w*cos(alpha*pi/180),R_w*sin(alpha*pi/180),alpha,1)
mi_drawarc(R_wi,0,R_wi*cos(alpha*pi/180),R_wi*sin(alpha*pi/180),alpha,1)
% winding back-iron boundaries
mi_drawarc(R_s,0,R_s*cos(alpha*pi/180),R_s*sin(alpha*pi/180),alpha,1)
mi_drawarc(R_se,0,R_se*cos(alpha*pi/180),R_se*sin(alpha*pi/180),alpha,1)
% external/internal boundaries
% mi_drawarc((R_se+sign(R_w-R_m)*Az_zero_2),0,(R_se+sign(R_w-R_m)*Az_zero_2)*cos(alpha*pi/180),(R_se+sign(R_w-R_m)*Az_zero_2)*sin(alpha*pi/180),alpha,1)
mi_drawarc(R_ie-sign(R_w-R_m)*Az_zero_1,0,(R_ie-sign(R_w-R_m)*Az_zero_1)*cos(alpha*pi/180),(R_ie-sign(R_w-R_m)*Az_zero_1)*sin(alpha*pi/180),alpha,1)
% slicing boundaries non-linear back-iron
% mi_drawline(R_ie-sign(R_w-R_m)*Az_zero_1,0,(R_se+sign(R_w-R_m)*Az_zero_2),0)
% mi_drawline((R_se+sign(R_w-R_m)*Az_zero_2)*cos(alpha*pi/180),(R_se+sign(R_w-R_m)*Az_zero_2)*sin(alpha*pi/180),(R_ie-sign(R_w-R_m)*Az_zero_1)*cos(alpha*pi/180),(R_ie-sign(R_w-R_m)*Az_zero_1)*sin(alpha*pi/180))
% slicing boundaries
mi_drawline(R_ie-sign(R_w-R_m)*Az_zero_1,0,R_se,0)
mi_drawline((R_se)*cos(alpha*pi/180),(R_se)*sin(alpha*pi/180),(R_ie-sign(R_w-R_m)*Az_zero_1)*cos(alpha*pi/180),(R_ie-sign(R_w-R_m)*Az_zero_1)*sin(alpha*pi/180))

% coils boundaries
 mi_selectsegment((R_w+sign(R_w-R_m)*s_t/2),0);
 mi_copyrotate2(0,0,alpha_s,Qsim-1,1)

% Magnets array drawing
mi_drawarc(R_m,0,R_m*cos(alpha*pi/180),R_m*sin(alpha*pi/180),alpha,1)

% Comment these three lines if a good demag map is not wanted
mi_selectarcsegment(R_m*cos(alpha/2*pi/180),R_m*sin(alpha/2*pi/180));
mi_selectarcsegment(R_r*cos(alpha/2*pi/180),R_r*sin(alpha/2*pi/180));
mi_setarcsegmentprop(alpha_span_end/20, 0, 0, 0)




mi_drawline(R_m*cos((1-alpha_p1)*pi/(2*p)),R_m*sin((1-alpha_p1)*pi/(2*p)),R_r*cos((1-alpha_p1)*pi/(2*p)),R_r*sin((1-alpha_p1)*pi/(2*p)))
mi_drawline(R_m*cos((1-alpha_p)*pi/(2*p)),R_m*sin((1-alpha_p)*pi/(2*p)),R_r*cos((1-alpha_p)*pi/(2*p)),R_r*sin((1-alpha_p)*pi/(2*p)))
mi_drawline(R_m*cos((1+alpha_p)*pi/(2*p)),R_m*sin((1+alpha_p)*pi/(2*p)),R_r*cos((1+alpha_p)*pi/(2*p)),R_r*sin((1+alpha_p)*pi/(2*p)))
mi_drawline(R_m*cos((1+alpha_p1)*pi/(2*p)),R_m*sin((1+alpha_p1)*pi/(2*p)),R_r*cos((1+alpha_p1)*pi/(2*p)),R_r*sin((1+alpha_p1)*pi/(2*p)))

if alpha_p == 0 % two-segments special Halbach
    mi_selectsegment((R_m-sign(R_w-R_m)*t_m/2)*cos(pi/(2*p)),(R_m-sign(R_w-R_m)*t_m/2)*sin(pi/(2*p)));
    mi_copyrotate2(0,0, alpha_span_side, dis_ang_side-2,1)
    mi_selectsegment((R_m-sign(R_w-R_m)*t_m/2)*cos(pi/(2*p)),(R_m-sign(R_w-R_m)*t_m/2)*sin(pi/(2*p)));
    mi_copyrotate2(0,0, -alpha_span_side, dis_ang_side-2,1)
    for ii = 1:2*dis_ang_side-2
        mi_selectsegment((R_m-sign(R_w-R_m)*t_m/2)*cos((ii-1)*alpha_span_side*pi/180),(R_m-sign(R_w-R_m)*t_m/2)*sin((ii-1)*alpha_span_side*pi/180));
    end
elseif alpha_p ~= alpha_p1 && alpha_p ~= 0 && Halbach == 0 % Even-segments Halbach (4-segments)
    mi_selectsegment((R_m-sign(R_w-R_m)*t_m/2)*cos((1-alpha_p)*pi/(2*p)),(R_m-sign(R_w-R_m)*t_m/2)*sin((1-alpha_p)*pi/(2*p)));
    mi_copyrotate2(0,0, -alpha_span_side, dis_ang_side-1,1)
    mi_selectsegment((R_m-sign(R_w-R_m)*t_m/2)*cos((1+alpha_p)*pi/(2*p)),(R_m-sign(R_w-R_m)*t_m/2)*sin((1+alpha_p)*pi/(2*p)));
    mi_copyrotate2(0,0, +alpha_span_side, dis_ang_side-1,1)
    mi_selectsegment((R_m-sign(R_w-R_m)*t_m/2)*cos((1-alpha_p1)*pi/(2*p)),(R_m-sign(R_w-R_m)*t_m/2)*sin((1-alpha_p1)*pi/(2*p)));
    mi_copyrotate2(0,0, -alpha_span_end, dis_ang_end-1,1)
    mi_selectsegment((R_m-sign(R_w-R_m)*t_m/2)*cos((1+alpha_p1)*pi/(2*p)),(R_m-sign(R_w-R_m)*t_m/2)*sin((1+alpha_p1)*pi/(2*p)));
    mi_copyrotate2(0,0, alpha_span_end, dis_ang_end-1,1)
    mi_selectsegment((R_m-sign(R_w-R_m)*t_m/2)*cos((1+alpha_p)*pi/(2*p)),(R_m-sign(R_w-R_m)*t_m/2)*sin((1+alpha_p)*pi/(2*p)));
    mi_copyrotate2(0,0, -alpha_span_mid, dis_ang_mid-2,1)
    for ii = 1:dis_ang_end-2
        mi_selectsegment((R_m-sign(R_w-R_m)*t_m/2)*cos((1+alpha_p1)*pi/(2*p)+ii*alpha_span_end*pi/180),(R_m-sign(R_w-R_m)*t_m/2)*sin((1+alpha_p1)*pi/(2*p)+ii*alpha_span_end*pi/180));
        mi_selectsegment((R_m-sign(R_w-R_m)*t_m/2)*cos((1-alpha_p1)*pi/(2*p)-ii*alpha_span_end*pi/180),(R_m-sign(R_w-R_m)*t_m/2)*sin((1-alpha_p1)*pi/(2*p)-ii*alpha_span_end*pi/180));
    end
    for ii = 1:dis_ang_side-1
        mi_selectsegment((R_m-sign(R_w-R_m)*t_m/2)*cos((1+alpha_p)*pi/(2*p)+ii*alpha_span_side*pi/180),(R_m-sign(R_w-R_m)*t_m/2)*sin((1+alpha_p)*pi/(2*p)+ii*alpha_span_side*pi/180));
        mi_selectsegment((R_m-sign(R_w-R_m)*t_m/2)*cos((1-alpha_p)*pi/(2*p)-ii*alpha_span_side*pi/180),(R_m-sign(R_w-R_m)*t_m/2)*sin((1-alpha_p)*pi/(2*p)-ii*alpha_span_side*pi/180));
    end
    for ii = 1:dis_ang_mid
        mi_selectsegment((R_m-sign(R_w-R_m)*t_m/2)*cos((1-alpha_p)*pi/(2*p)+(ii-1)*alpha_span_mid*pi/180),(R_m-sign(R_w-R_m)*t_m/2)*sin((1-alpha_p)*pi/(2*p)+(ii-1)*alpha_span_mid*pi/180));
    end
    
elseif alpha_p ~= alpha_p1 && alpha_p ~= 0 && Halbach == 1
    mi_selectsegment((R_m-sign(R_w-R_m)*t_m/2)*cos((1-alpha_p)*pi/(2*p)),(R_m-sign(R_w-R_m)*t_m/2)*sin((1-alpha_p)*pi/(2*p)));
    mi_copyrotate2(0,0, -alpha_span_side, dis_ang_side-1,1)
    mi_selectsegment((R_m-sign(R_w-R_m)*t_m/2)*cos((1+alpha_p)*pi/(2*p)),(R_m-sign(R_w-R_m)*t_m/2)*sin((1+alpha_p)*pi/(2*p)));
    mi_copyrotate2(0,0, +alpha_span_side, dis_ang_side-1,1)
    mi_selectsegment((R_m-sign(R_w-R_m)*t_m/2)*cos((1-alpha_p1)*pi/(2*p)),(R_m-sign(R_w-R_m)*t_m/2)*sin((1-alpha_p1)*pi/(2*p)));
    mi_copyrotate2(0,0, -alpha_span_end, dis_ang_end-1,1)
    mi_selectsegment((R_m-sign(R_w-R_m)*t_m/2)*cos((1+alpha_p1)*pi/(2*p)),(R_m-sign(R_w-R_m)*t_m/2)*sin((1+alpha_p1)*pi/(2*p)));
    mi_copyrotate2(0,0, alpha_span_end, dis_ang_end-1,1)
    mi_selectsegment((R_m-sign(R_w-R_m)*t_m/2)*cos((1+alpha_p)*pi/(2*p)),(R_m-sign(R_w-R_m)*t_m/2)*sin((1+alpha_p)*pi/(2*p)));
    mi_copyrotate2(0,0, -alpha_span_mid, dis_ang_mid-2,1)
    for ii = 1:dis_ang_end-2
        mi_selectsegment((R_m-sign(R_w-R_m)*t_m/2)*cos((1+alpha_p1)*pi/(2*p)+ii*alpha_span_end*pi/180),(R_m-sign(R_w-R_m)*t_m/2)*sin((1+alpha_p1)*pi/(2*p)+ii*alpha_span_end*pi/180));
        mi_selectsegment((R_m-sign(R_w-R_m)*t_m/2)*cos((1-alpha_p1)*pi/(2*p)-ii*alpha_span_end*pi/180),(R_m-sign(R_w-R_m)*t_m/2)*sin((1-alpha_p1)*pi/(2*p)-ii*alpha_span_end*pi/180));
    end
    for ii = 1:dis_ang_side-1
        mi_selectsegment((R_m-sign(R_w-R_m)*t_m/2)*cos((1+alpha_p)*pi/(2*p)+ii*alpha_span_side*pi/180),(R_m-sign(R_w-R_m)*t_m/2)*sin((1+alpha_p)*pi/(2*p)+ii*alpha_span_side*pi/180));
        mi_selectsegment((R_m-sign(R_w-R_m)*t_m/2)*cos((1-alpha_p)*pi/(2*p)-ii*alpha_span_side*pi/180),(R_m-sign(R_w-R_m)*t_m/2)*sin((1-alpha_p)*pi/(2*p)-ii*alpha_span_side*pi/180));
    end
    for ii = 1:dis_ang_mid
        mi_selectsegment((R_m-sign(R_w-R_m)*t_m/2)*cos((1-alpha_p)*pi/(2*p)+(ii-1)*alpha_span_mid*pi/180),(R_m-sign(R_w-R_m)*t_m/2)*sin((1-alpha_p)*pi/(2*p)+(ii-1)*alpha_span_mid*pi/180));
    end




elseif alpha_p == alpha_p1 % 2-segments Halbach
    mi_selectsegment((R_m-sign(R_w-R_m)*t_m/2)*cos((1-alpha_p)*pi/(2*p)),(R_m-sign(R_w-R_m)*t_m/2)*sin((1-alpha_p)*pi/(2*p)));
    mi_copyrotate2(0,0, -alpha_span_end, dis_ang_end-1,1)
    mi_selectsegment((R_m-sign(R_w-R_m)*t_m/2)*cos((1-alpha_p)*pi/(2*p)),(R_m-sign(R_w-R_m)*t_m/2)*sin((1-alpha_p)*pi/(2*p)));
    mi_copyrotate2(0,0, alpha_span_mid, dis_ang_mid-2,1)
    mi_selectsegment((R_m-sign(R_w-R_m)*t_m/2)*cos((1+alpha_p)*pi/(2*p)),(R_m-sign(R_w-R_m)*t_m/2)*sin((1+alpha_p)*pi/(2*p)));
    mi_copyrotate2(0,0, alpha_span_end, dis_ang_end-1,1)
    for ii = 1:dis_ang_end-2
        mi_selectsegment((R_m-sign(R_w-R_m)*t_m/2)*cos((1+alpha_p)*pi/(2*p)+ii*alpha_span_end*pi/180),(R_m-sign(R_w-R_m)*t_m/2)*sin((1+alpha_p)*pi/(2*p)+ii*alpha_span_end*pi/180));
        mi_selectsegment((R_m-sign(R_w-R_m)*t_m/2)*cos((1-alpha_p)*pi/(2*p)-ii*alpha_span_end*pi/180),(R_m-sign(R_w-R_m)*t_m/2)*sin((1-alpha_p)*pi/(2*p)-ii*alpha_span_end*pi/180));
    end
    for ii = 1:dis_ang_mid
        mi_selectsegment((R_m-sign(R_w-R_m)*t_m/2)*cos((1-alpha_p)*pi/(2*p)+(ii-1)*alpha_span_mid*pi/180),(R_m-sign(R_w-R_m)*t_m/2)*sin((1-alpha_p)*pi/(2*p)+(ii-1)*alpha_span_mid*pi/180));
    end
end

mi_setsegmentprop(0, rad_thick/20, 0, 0, 0)
mi_copyrotate2(0,0, theta_p, np-1,1)

for ii = 1:rad_dis
   mi_drawarc(r_m_dis(ii),0,r_m_dis(ii)*cos(alpha*pi/180),r_m_dis(ii)*sin(alpha*pi/180),alpha,1) 
end



%% domains definition

% Magnets domains
if alpha_p == 0
    for ii = 1:rad_dis-1
        mi_addblocklabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos(pi/(2*p)-alpha_span_side/2*pi/180),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin(pi/(2*p)-alpha_span_side/2*pi/180))
    end
    for ii = 1:rad_dis-1
        mi_selectlabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos(pi/(2*p)-alpha_span_side/2*pi/180),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin(pi/(2*p)-alpha_span_side/2*pi/180));
    end
    mi_copyrotate2(0,0, -alpha_span_side, dis_ang_side-2,2)
    
    for ii = 1:rad_dis-1
        mi_addblocklabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos(pi/(2*p)+alpha_span_side/2*pi/180),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin(pi/(2*p)+alpha_span_side/2*pi/180))
    end
    for ii = 1:rad_dis-1
        mi_selectlabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos(pi/(2*p)+alpha_span_side/2*pi/180),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin(pi/(2*p)+alpha_span_side/2*pi/180));
    end
    mi_copyrotate2(0,0, alpha_span_side, dis_ang_side-2,2)
   
    mi_selectcircle(0,0,max(R_r,R_m),2);
    mi_copyrotate2(0,0, theta_p, np-1,2)
    
elseif alpha_p ~= alpha_p1 && alpha_p ~= 0 && Halbach == 1
    for ii = 1:rad_dis-1
        mi_addblocklabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1-alpha_p1)*pi/(2*p)-alpha_span_end/2*pi/180),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1-alpha_p1)*pi/(2*p)-alpha_span_end/2*pi/180))
        mi_addblocklabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1+alpha_p1)*pi/(2*p)+alpha_span_end/2*pi/180),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1+alpha_p1)*pi/(2*p)+alpha_span_end/2*pi/180))
    end
    for ii = 1:rad_dis-1
        mi_selectlabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1-alpha_p1)*pi/(2*p)-alpha_span_end/2*pi/180),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1-alpha_p1)*pi/(2*p)-alpha_span_end/2*pi/180));
    end
    mi_copyrotate2(0,0, -alpha_span_end, dis_ang_end-2,2)
    for ii = 1:rad_dis-1
        mi_selectlabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1+alpha_p1)*pi/(2*p)+alpha_span_end/2*pi/180),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1+alpha_p1)*pi/(2*p)+alpha_span_end/2*pi/180));
    end
    mi_copyrotate2(0,0, +alpha_span_end, dis_ang_end-2,2)
    
    for ii = 1:rad_dis-1
        mi_addblocklabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1-alpha_p)*pi/(2*p)-alpha_span_side/2*pi/180),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1-alpha_p)*pi/(2*p)-alpha_span_side/2*pi/180))
        mi_addblocklabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1+alpha_p)*pi/(2*p)+alpha_span_side/2*pi/180),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1+alpha_p)*pi/(2*p)+alpha_span_side/2*pi/180))
    end
    for ii = 1:rad_dis-1
        mi_selectlabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1-alpha_p)*pi/(2*p)-alpha_span_side/2*pi/180),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1-alpha_p)*pi/(2*p)-alpha_span_side/2*pi/180));
    end
    mi_copyrotate2(0,0, -alpha_span_side, dis_ang_side-2,2)
    for ii = 1:rad_dis-1
        mi_selectlabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1+alpha_p)*pi/(2*p)+alpha_span_side/2*pi/180),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1+alpha_p)*pi/(2*p)+alpha_span_side/2*pi/180));
    end
    mi_copyrotate2(0,0, +alpha_span_side, dis_ang_side-2,2)
    
    for ii = 1:rad_dis-1
        mi_addblocklabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1-alpha_p)*pi/(2*p)+alpha_span_mid/2*pi/180),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1-alpha_p)*pi/(2*p)+alpha_span_mid/2*pi/180))
    end
    for ii = 1:rad_dis-1
        mi_selectlabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1-alpha_p)*pi/(2*p)+alpha_span_mid/2*pi/180),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1-alpha_p)*pi/(2*p)+alpha_span_mid/2*pi/180));
    end
    mi_copyrotate2(0,0, +alpha_span_mid, dis_ang_mid-2,2)
    
    mi_selectcircle(0,0,max(R_r,R_m),2);
    mi_copyrotate2(0,0, theta_p, np-1,2)
    
elseif alpha_p ~= alpha_p1 && alpha_p ~= 0 && Halbach == 0
   for ii = 1:rad_dis-1
        mi_addblocklabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1-alpha_p1)*pi/(2*p)-alpha_span_end/2*pi/180),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1-alpha_p1)*pi/(2*p)-alpha_span_end/2*pi/180))
        mi_addblocklabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1+alpha_p1)*pi/(2*p)+alpha_span_end/2*pi/180),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1+alpha_p1)*pi/(2*p)+alpha_span_end/2*pi/180))
    end
    for ii = 1:rad_dis-1
        mi_selectlabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1-alpha_p1)*pi/(2*p)-alpha_span_end/2*pi/180),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1-alpha_p1)*pi/(2*p)-alpha_span_end/2*pi/180));
    end
    mi_copyrotate2(0,0, -alpha_span_end, dis_ang_end-2,2)
    for ii = 1:rad_dis-1
        mi_selectlabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1+alpha_p1)*pi/(2*p)+alpha_span_end/2*pi/180),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1+alpha_p1)*pi/(2*p)+alpha_span_end/2*pi/180));
    end
    mi_copyrotate2(0,0, +alpha_span_end, dis_ang_end-2,2)
    
    for ii = 1:rad_dis-1
        mi_addblocklabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1-alpha_p)*pi/(2*p)-alpha_span_side/2*pi/180),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1-alpha_p)*pi/(2*p)-alpha_span_side/2*pi/180))
        mi_addblocklabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1+alpha_p)*pi/(2*p)+alpha_span_side/2*pi/180),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1+alpha_p)*pi/(2*p)+alpha_span_side/2*pi/180))
    end
    for ii = 1:rad_dis-1
        mi_selectlabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1-alpha_p)*pi/(2*p)-alpha_span_side/2*pi/180),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1-alpha_p)*pi/(2*p)-alpha_span_side/2*pi/180));
    end
    mi_copyrotate2(0,0, -alpha_span_side, dis_ang_side-2,2)
    for ii = 1:rad_dis-1
        mi_selectlabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1+alpha_p)*pi/(2*p)+alpha_span_side/2*pi/180),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1+alpha_p)*pi/(2*p)+alpha_span_side/2*pi/180));
    end
    mi_copyrotate2(0,0, +alpha_span_side, dis_ang_side-2,2)
    
    for ii = 1:rad_dis-1
        mi_addblocklabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1-alpha_p)*pi/(2*p)+alpha_span_mid/2*pi/180),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1-alpha_p)*pi/(2*p)+alpha_span_mid/2*pi/180))
    end
    for ii = 1:rad_dis-1
        mi_selectlabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1-alpha_p)*pi/(2*p)+alpha_span_mid/2*pi/180),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1-alpha_p)*pi/(2*p)+alpha_span_mid/2*pi/180));
    end
    mi_copyrotate2(0,0, +alpha_span_mid, dis_ang_mid-2,2)
    
    mi_selectcircle(0,0,max(R_r,R_m),2);
    mi_copyrotate2(0,0, theta_p, np-1,2)
   
elseif alpha_p == alpha_p1 % 2-segments Halbach array (mid+end magnets)
    for ii = 1:rad_dis-1
        mi_addblocklabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1-alpha_p)*pi/(2*p)-alpha_span_end/2*pi/180),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1-alpha_p)*pi/(2*p)-alpha_span_end/2*pi/180))
        mi_addblocklabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1+alpha_p)*pi/(2*p)+alpha_span_end/2*pi/180),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1+alpha_p)*pi/(2*p)+alpha_span_end/2*pi/180))
    end
    for ii = 1:rad_dis-1
        mi_selectlabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1-alpha_p)*pi/(2*p)-alpha_span_end/2*pi/180),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1-alpha_p)*pi/(2*p)-alpha_span_end/2*pi/180));
    end
    mi_copyrotate2(0,0, -alpha_span_end, dis_ang_end-2,2)
    for ii = 1:rad_dis-1
        mi_selectlabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1+alpha_p)*pi/(2*p)+alpha_span_end/2*pi/180),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1+alpha_p)*pi/(2*p)+alpha_span_end/2*pi/180));
    end
    mi_copyrotate2(0,0, +alpha_span_end, dis_ang_end-2,2)
    
    for ii = 1:rad_dis-1
        mi_addblocklabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1-alpha_p)*pi/(2*p)+alpha_span_mid/2*pi/180),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1-alpha_p)*pi/(2*p)+alpha_span_mid/2*pi/180))
    end
    for ii = 1:rad_dis-1
        mi_selectlabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1-alpha_p)*pi/(2*p)+alpha_span_mid/2*pi/180),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1-alpha_p)*pi/(2*p)+alpha_span_mid/2*pi/180));
    end
    mi_copyrotate2(0,0, +alpha_span_mid, dis_ang_mid-2,2)
    
    mi_selectcircle(0,0,max(R_r,R_m),2);
    mi_copyrotate2(0,0, theta_p, np-1,2)
   
end


% Inner air
mi_addblocklabel((R_ie-sign(R_w-R_m)*0.05*R_ie)*cos(alpha/2*pi/180),(R_ie-sign(R_w-R_m)*0.05*R_ie)*sin(alpha/2*pi/180))
% Outer air if non-linear iron
% mi_addblocklabel((R_se+sign(R_w-R_m)*0.05*R_se)*cos(alpha/2*pi/180),(R_se+sign(R_w-R_m)*0.05*R_se)*sin(alpha/2*pi/180))
% Inner rotor back-iron
mi_addblocklabel((R_s-(R_s-R_se)/2)*cos(alpha/2*pi/180),(R_s-(R_s-R_se)/2)*sin(alpha/2*pi/180))
% Outer rotor back-iron domain
mi_addblocklabel((R_i-(R_i-R_ie)/2)*cos(alpha/2*pi/180),(R_i-(R_i-R_ie)/2)*sin(alpha/2*pi/180))
% Inner air-gaps
if R_wi~=R_s
mi_addblocklabel(((R_wi+sign(R_w-R_m)*gap_w/2))*cos(alpha/2*pi/180),((R_wi+sign(R_w-R_m)*gap_w/2))*sin(alpha/2*pi/180))
end
mi_addblocklabel(((R_w-sign(R_w-R_m)*gap_m/2))*cos(alpha/2*pi/180),((R_w-sign(R_w-R_m)*gap_m/2))*sin(alpha/2*pi/180))
% Coils domains
 for iq = 1:Qsim
        mi_addblocklabel((R_w+sign(R_w-R_m)*s_t/2)*cos((2*iq-1)*alpha_s/2*pi/180),(R_w+sign(R_w-R_m)*s_t/2)*sin((2*iq-1)*alpha_s/2*pi/180))
end



%% materials definition
mi_addmaterial('AIR')
mi_addmaterial('COPPER', 1, 1, 0, 0, 58)
% iron BH curve points definition
mi_addmaterial('IRON',100000,100000,0 ,0,0)
mi_addmaterial('M330-50A')
mi_addbhpoint('M330-50A',0.501378,53.000000)
mi_addbhpoint('M330-50A',0.596626,60.000000)
mi_addbhpoint('M330-50A',0.700687,69.000000)
mi_addbhpoint('M330-50A',0.795684,79.000000)
mi_addbhpoint('M330-50A',0.899644,93.000000)
mi_addbhpoint('M330-50A',1.000826,112.000000)
mi_addbhpoint('M330-50A',1.104045,141.000000)
mi_addbhpoint('M330-50A',1.200822,192.000000)
mi_addbhpoint('M330-50A',1.302307,313.000000)
mi_addbhpoint('M330-50A',1.400781,698.000000)
mi_addbhpoint('M330-50A',1.502822,1932.000000)
mi_addbhpoint('M330-50A',1.604263,4284.000000)
mi_addbhpoint('M330-50A',1.714053,7750.000000)
mi_addbhpoint('M330-50A',1.824212,13318.000000)
mi_addbhpoint('M330-50A',1.922127,19610.000000)
mi_addbhpoint('M330-50A',2,50000.000000)
mi_addbhpoint('M330-50A',2.05,200000.000000)
mi_addbhpoint('M330-50A',2.1,1000000.000000)
mi_addmaterial('ALUMINUM', 1, 1, 0, 0, 34.45)
% mi_addmaterial('N30EH@100C', 1, 1, B_r/(mu_r*mu_0), 0, 0.625)

%% definition of boundary conditions
mi_addboundprop('Az_zero')
if mod(np,2)==1
    mi_addboundprop('Apbc_1', 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0)
    mi_addboundprop('Apbc_2', 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0)
    mi_addboundprop('Apbc_3', 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0)
    mi_addboundprop('Apbc_4', 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0)
    mi_addboundprop('Apbc_5', 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0)
    mi_addboundprop('Apbc_6', 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0)
    mi_addboundprop('Apbc_7', 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0)
    mi_addboundprop('Apbc_8', 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0)
    for ii = 1:rad_dis
        mi_addboundprop(['Apbc_', num2str(8+ii)], 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0)       
    end
%     mi_addboundprop('Apbc_9', 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0)
    mi_addboundprop('Apagc_o', 0, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0)
elseif mod(np,2)==0
    mi_addboundprop('Pbc_1', 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0)
    mi_addboundprop('Pbc_2', 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0)
    mi_addboundprop('Pbc_3', 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0)
    mi_addboundprop('Pbc_4', 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0)
    mi_addboundprop('Pbc_5', 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0)
    mi_addboundprop('Pbc_6', 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0)
    mi_addboundprop('Pbc_7', 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0)
    mi_addboundprop('Pbc_8', 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0)
    for ii = 1:rad_dis-1
        mi_addboundprop(['Pbc_', num2str(8+ii)], 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0)       
    end
    mi_addboundprop('Pagc_o', 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0)
end



%% electric sources definition 
 for iq = 1:Qsim
        mi_addcircprop(['Is', num2str(iq)], 0, 1)       
end


%% assign boundary condition
% non-linear iron case
% mi_selectarcsegment((R_se+sign(R_w-R_m)*Az_zero_2)*cos(alpha/2*pi/180),(R_se+sign(R_w-R_m)*Az_zero_2)*sin(alpha/2*pi/180));
% linear iron case
mi_selectarcsegment((R_se)*cos(alpha/2*pi/180),(R_se)*sin(alpha/2*pi/180));
mi_selectarcsegment((R_ie-sign(R_w-R_m)*Az_zero_1)*cos(alpha/2*pi/180),(R_ie-sign(R_w-R_m)*Az_zero_1)*sin(alpha/2*pi/180));
mi_setarcsegmentprop(1, 'Az_zero', 0, 0)
mi_clearselected

if mod(np,2)==1 
mi_selectsegment((R_ie-sign(R_w-R_m)*0.05*R_ie),0);
mi_selectsegment((R_ie-sign(R_w-R_m)*0.05*R_ie)*cos(alpha*pi/180),(R_ie-sign(R_w-R_m)*0.05*R_ie)*sin(alpha*pi/180));
mi_setsegmentprop('Apbc_1', 0, 1, 0, 0)
mi_clearselected
mi_selectsegment((R_s-(R_s-R_se)/2),0);
mi_selectsegment((R_s-(R_s-R_se)/2)*cos(alpha*pi/180),(R_s-(R_s-R_se)/2)*sin(alpha*pi/180));
mi_setsegmentprop('Apbc_2', 0, 1, 0, 0)
mi_clearselected
mi_selectsegment((R_i-(R_i-R_ie)/2),0);
mi_selectsegment((R_i-(R_i-R_ie)/2)*cos(alpha*pi/180),(R_i-(R_i-R_ie)/2)*sin(alpha*pi/180));
mi_setsegmentprop('Apbc_3', 0, 1, 0, 0)
mi_clearselected
mi_selectsegment((R_wi+sign(R_w-R_m)*gap_w/2),0);
mi_selectsegment((R_wi+sign(R_w-R_m)*gap_w/2)*cos(alpha*pi/180),(R_wi+sign(R_w-R_m)*gap_w/2)*sin(alpha*pi/180));
mi_setsegmentprop('Apbc_4', 0, 1, 0, 0)
mi_clearselected
mi_selectsegment(((R_w-sign(R_w-R_m)*gap_m/2)-gap_m/4),0);
mi_selectsegment(((R_w-sign(R_w-R_m)*gap_m/2)-gap_m/4)*cos(alpha*pi/180),((R_w-sign(R_w-R_m)*gap_m/2)-gap_m/4)*sin(alpha*pi/180));
mi_setsegmentprop('Apbc_5', 0, 1, 0, 0)
mi_clearselected
mi_selectsegment((R_w+sign(R_w-R_m)*s_t/2),0);
mi_selectsegment((R_w+sign(R_w-R_m)*s_t/2)*cos(alpha*pi/180),(R_w+sign(R_w-R_m)*s_t/2)*sin(alpha*pi/180));
mi_setsegmentprop('Apbc_7', 0, 1, 0, 0)
mi_clearselected
mi_selectsegment((R_se+sign(R_w-R_m)*0.05*R_se),0);
mi_selectsegment((R_se+sign(R_w-R_m)*0.05*R_se)*cos(alpha*pi/180),(R_se+sign(R_w-R_m)*0.05*R_se)*sin(alpha*pi/180));
mi_setsegmentprop('Apbc_8', 0, 1, 0, 0)
mi_clearselected
for ii=1:rad_dis-1
mi_selectsegment((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick),0);
mi_selectsegment((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos(alpha*pi/180),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin(alpha*pi/180));
mi_setsegmentprop(['Apbc_', num2str(8+ii)], rad_thick/20, 0, 0, 0) 
mi_clearselected
end

elseif mod(np,2)==0 
mi_selectarcsegment(((R_w-sign(R_w-R_m)*gap_m/2)-gap_m/8)*cos(alpha/2*pi/180),((R_w-sign(R_w-R_m)*gap_m/2)-gap_m/8)*sin(alpha/2*pi/180));
mi_selectarcsegment(((R_w-sign(R_w-R_m)*gap_m/2)+gap_m/8)*cos(alpha/2*pi/180),((R_w-sign(R_w-R_m)*gap_m/2)+gap_m/8)*sin(alpha/2*pi/180));
mi_setarcsegmentprop(1, 'Pagc_o', 0, 0)
mi_clearselected
mi_selectsegment((R_ie-sign(R_w-R_m)*0.05*R_ie),0);
mi_selectsegment((R_ie-sign(R_w-R_m)*0.05*R_ie)*cos(alpha*pi/180),(R_ie-sign(R_w-R_m)*0.05*R_ie)*sin(alpha*pi/180));
mi_setsegmentprop('Pbc_1', 0, 1, 0, 0)
mi_clearselected
mi_selectsegment((R_s-(R_s-R_se)/2),0);
mi_selectsegment((R_s-(R_s-R_se)/2)*cos(alpha*pi/180),(R_s-(R_s-R_se)/2)*sin(alpha*pi/180));
mi_setsegmentprop('Pbc_2', 0, 1, 0, 0)
mi_clearselected
mi_selectsegment((R_i-(R_i-R_ie)/2),0);
mi_selectsegment((R_i-(R_i-R_ie)/2)*cos(alpha*pi/180),(R_i-(R_i-R_ie)/2)*sin(alpha*pi/180));
mi_setsegmentprop('Pbc_3', 0, 1, 0, 0)
mi_clearselected
mi_selectsegment((R_wi+sign(R_w-R_m)*gap_w/2),0);
mi_selectsegment((R_wi+sign(R_w-R_m)*gap_w/2)*cos(alpha*pi/180),(R_wi+sign(R_w-R_m)*gap_w/2)*sin(alpha*pi/180));
mi_setsegmentprop('Pbc_4', 0, 1, 0, 0)
mi_clearselected
mi_selectsegment(((R_w-sign(R_w-R_m)*gap_m/2)-gap_m/4),0);
mi_selectsegment(((R_w-sign(R_w-R_m)*gap_m/2)-gap_m/4)*cos(alpha*pi/180),((R_w-sign(R_w-R_m)*gap_m/2)-gap_m/4)*sin(alpha*pi/180));
mi_setsegmentprop('Pbc_5', 0, 1, 0, 0)
mi_clearselected
mi_selectsegment(((R_w-sign(R_w-R_m)*gap_m/2)+gap_m/4),0);
mi_selectsegment(((R_w-sign(R_w-R_m)*gap_m/2)+gap_m/4)*cos(alpha*pi/180),((R_w-sign(R_w-R_m)*gap_m/2)+gap_m/4)*sin(alpha*pi/180));
mi_setsegmentprop('Pbc_6', 0, 1, 0, 0)
mi_clearselected
mi_selectsegment((R_w+sign(R_w-R_m)*s_t/2),0);
mi_selectsegment((R_w+sign(R_w-R_m)*s_t/2)*cos(alpha*pi/180),(R_w+sign(R_w-R_m)*s_t/2)*sin(alpha*pi/180));
mi_setsegmentprop('Pbc_7', 0, 1, 0, 0)
mi_clearselected
mi_selectsegment((R_se+sign(R_w-R_m)*0.05*R_se),0);
mi_selectsegment((R_se+sign(R_w-R_m)*0.05*R_se)*cos(alpha*pi/180),(R_se+sign(R_w-R_m)*0.05*R_se)*sin(alpha*pi/180));
mi_setsegmentprop('Pbc_8', 0, 1, 0, 0)
mi_clearselected
for ii=1:rad_dis-1
mi_selectsegment((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick),0);
mi_selectsegment((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos(alpha*pi/180),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin(alpha*pi/180));
mi_setsegmentprop(['Pbc_', num2str(8+ii)], rad_thick/20, 0, 0, 0) 
mi_clearselected
end

end


%% assign material properties and sources
mi_selectlabel((R_ie-sign(R_w-R_m)*0.05*R_ie)*cos(alpha/2*pi/180),(R_ie-sign(R_w-R_m)*0.05*R_ie)*sin(alpha/2*pi/180));
mi_setblockprop('AIR', 1, rad_thick/20, 0,0, 0,1)
mi_clearselected
% inner rotor
mi_selectlabel((R_s-(R_s-R_se)/2)*cos(alpha/2*pi/180),(R_s-(R_s-R_se)/2)*sin(alpha/2*pi/180));
mi_setblockprop('IRON', 1, 0, 0,0, 0,1)
mi_clearselected
% double air-gap
if R_wi~=R_s
mi_selectlabel(((R_wi+sign(R_w-R_m)*gap_w/2))*cos(alpha/2*pi/180),((R_wi+sign(R_w-R_m)*gap_w/2))*sin(alpha/2*pi/180));
mi_setblockprop('AIR', 1, 0, 0,0, 0,1)
mi_clearselected
end
mi_selectlabel(((R_w-sign(R_w-R_m)*gap_m/2))*cos(alpha/2*pi/180),((R_w-sign(R_w-R_m)*gap_m/2))*sin(alpha/2*pi/180));
mi_setblockprop('AIR', 1, rad_thick/20, 0,0, 0,1)
mi_clearselected
 
 c = 0;
 
first_side (1,np) = 0; 
second_side (1,np) = 0;
mid (1,np) = 0; 
first_end (1,np) = 0;
second_end (1,np) = 0;
for im = 1:np
    if alpha_p == 0
        for ii = 1:rad_dis-1
            for jj = 1:dis_ang_side-1
                c = c+1;
                mi_addmaterial(['PM',num2str(c)], mu_r, mu_r, B_r/(mu_r*mu_0), 0, 0.625)
                mi_selectlabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos(pi/(2*p)-alpha_span_side/2*pi/180-(jj-1)*alpha_span_side*pi/180 + (im-1)*pi/p),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin(pi/(2*p)-alpha_span_side/2*pi/180-(jj-1)*alpha_span_side*pi/180 + (im-1)*pi/p));
                mi_setblockprop(['PM',num2str(c)], 1, rad_thick/20, 0 , num2str((90/p-delta_m*180/pi+180/p*(im-1))+sign(R_w-R_m)*(90-theta_m_side*180/pi)+sign(R_w-R_m)*180*(im-1)), c ,1)
                mi_clearselected
            end
        end
        
        first_side(1,im) = c; % maximum group number of the first side-magnet over each pole 
        
        for ii = 1:rad_dis-1
            for jj = 1:dis_ang_side-1
                c = c+1;
                mi_addmaterial(['PM',num2str(c)], mu_r, mu_r, B_r/(mu_r*mu_0), 0, 0.625)
                mi_selectlabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos(pi/(2*p)+alpha_span_side/2*pi/180+(jj-1)*alpha_span_side*pi/180 + (im-1)*pi/p),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin(pi/(2*p)+alpha_span_side/2*pi/180+(jj-1)*alpha_span_side*pi/180 + (im-1)*pi/p));
                mi_setblockprop(['PM',num2str(c)], 1, rad_thick/20, 0 , num2str((90/p+delta_m*180/pi+180/p*(im-1))-sign(R_w-R_m)*(90-theta_m_side*180/pi)+sign(R_w-R_m)*180*(im-1)), c,1)
                mi_clearselected
            end
        end
        
        second_side(1,im) = c; % maximum group number of the first second-magnet over each pole 
    end
        
        if Halbach == 0 && alpha_p ~= 0
            
            for ii = 1:rad_dis-1
                for jj = 1:dis_ang_end-1
                    c = c+1;
                    mi_addmaterial(['PM',num2str(c)], mu_r, mu_r, B_r/(mu_r*mu_0), 0, 0.625)
                    mi_selectlabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1-alpha_p1)*pi/(2*p)-alpha_span_end/2*pi/180-(jj-1)*alpha_span_end*pi/180 + (im-1)*pi/p),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1-alpha_p1)*pi/(2*p)-alpha_span_end/2*pi/180-(jj-1)*alpha_span_end*pi/180 + (im-1)*pi/p));
                    mi_setblockprop(['PM',num2str(c)], 1, rad_thick/20, 0 , num2str(180/p*(im-1)+sign(R_w-R_m)*90+180*(im-1)), c,1)
                    mi_clearselected
                end
            end
            
            first_end(1,im) = c; % maximum group number of the first side-magnet over each pole
            
            if Halbach_1 == 4
                for ii = 1:rad_dis-1
                    for jj = 1:dis_ang_side-1
                        c = c+1;
                        mi_addmaterial(['PM',num2str(c)], mu_r, mu_r, B_r/(mu_r*mu_0), 0, 0.625)
                        mi_selectlabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1-alpha_p)*pi/(2*p)-alpha_span_side/2*pi/180-(jj-1)*alpha_span_side*pi/180 + (im-1)*pi/p),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1-alpha_p)*pi/(2*p)-alpha_span_side/2*pi/180-(jj-1)*alpha_span_side*pi/180 + (im-1)*pi/p));
                        mi_setblockprop(['PM',num2str(c)], 1, rad_thick/20, 0 , num2str((90/p-delta_m*180/pi+180/p*(im-1))+sign(R_w-R_m)*(90-theta_m_side*180/pi)+sign(R_w-R_m)*180*(im-1)), c,1)
                        mi_clearselected
                    end
                end
                
                first_side(1,im) = c; % maximum group number of the first side-magnet over each pole
            end
            
            for ii = 1:rad_dis-1
                for jj = 1:dis_ang_mid-1
                    c = c+1;
                    mi_addmaterial(['PM',num2str(c)], mu_r, mu_r, B_r/(mu_r*mu_0), 0, 0.625)
                    mi_selectlabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1-alpha_p)*pi/(2*p)+alpha_span_mid/2*pi/180+(jj-1)*alpha_span_mid*pi/180 + (im-1)*pi/(p)),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1-alpha_p)*pi/(2*p)+alpha_span_mid/2*pi/180+(jj-1)*alpha_span_mid*pi/180 + (im-1)*pi/(p)));
                    mi_setblockprop(['PM',num2str(c)], 1, rad_thick/20, 0 , num2str(180/(2*p)+(im-1)*180/p+sign(R_w-R_m)*(im-1)*180), c,1)
                    mi_clearselected
                end
            end
            
            mid(1,im) = c; % maximum group number of the mid-magnet over each pole
            
            if Halbach_1 == 4
                for ii = 1:rad_dis-1
                    for jj = 1:dis_ang_side-1
                        c = c+1;
                        mi_addmaterial(['PM',num2str(c)], mu_r, mu_r, B_r/(mu_r*mu_0), 0, 0.625)
                        mi_selectlabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1+alpha_p)*pi/(2*p)+alpha_span_side/2*pi/180+(jj-1)*alpha_span_side*pi/180 + (im-1)*pi/p),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1+alpha_p)*pi/(2*p)+alpha_span_side/2*pi/180+(jj-1)*alpha_span_side*pi/180 + (im-1)*pi/p));
                        mi_setblockprop(['PM',num2str(c)], 1, rad_thick/20, 0 , num2str((90/p+delta_m*180/pi+180/p*(im-1))-sign(R_w-R_m)*(90-theta_m_side*180/pi)+sign(R_w-R_m)*180*(im-1)), c,1)
                        mi_clearselected
                    end
                end
                
                second_side(1,im) = c; % maximum group number of the first side-magnet over each pole
            end
            
            for ii = 1:rad_dis-1
                for jj = 1:dis_ang_end-1
                    c = c+1;
                    mi_addmaterial(['PM',num2str(c)], mu_r, mu_r, B_r/(mu_r*mu_0), 0, 0.625)
                    mi_selectlabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1+alpha_p1)*pi/(2*p)+alpha_span_end/2*pi/180+(jj-1)*alpha_span_end*pi/180 + (im-1)*pi/p),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1+alpha_p1)*pi/(2*p)+alpha_span_end/2*pi/180+(jj-1)*alpha_span_end*pi/180 + (im-1)*pi/p));
                    mi_setblockprop(['PM',num2str(c)], 1, rad_thick/20, 0 , num2str(180/p-sign(R_w-R_m)*90 +(im-1)*(180+180/p)), c,1)
                    mi_clearselected
                end
            end
            
            second_end(1,im) = c; % maximum group number of the first side-magnet over each pole
        end
        
        
            
    if Halbach == 1
        if Halbach_1 ==5
            for ii = 1:rad_dis-1
                    for jj = 1:dis_ang_end-1
                        c = c+1;
                        mi_addmaterial(['PM',num2str(c)], mu_r, mu_r, B_r/(mu_r*mu_0), 0, 0.625)
                        mi_selectlabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1-alpha_p1)*pi/(2*p)-alpha_span_end/2*pi/180-(jj-1)*alpha_span_end*pi/180 + (im-1)*(pi/p)),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1-alpha_p1)*pi/(2*p)-alpha_span_end/2*pi/180-(jj-1)*alpha_span_end*pi/180 + (im-1)*(pi/p)));
                        mi_setblockprop(['PM',num2str(c)], 1, 0, 0 , num2str((90/p-delta_m1*180/pi+180/p*(im-1))+sign(R_w-R_m)*(90-theta_m_end*180/pi)+sign(R_w-R_m)*180*(im-1)), c,1)
                        mi_clearselected
                    end
            end
            
            first_end(1,im) = c; % maximum group number of the first side-magnet over each pole
        end
            
                for ii = 1:rad_dis-1
                    for jj = 1:dis_ang_side-1
                        c = c+1;
                        mi_addmaterial(['PM',num2str(c)], mu_r, mu_r, B_r/(mu_r*mu_0), 0, 0.625)
                        mi_selectlabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1-alpha_p)*pi/(2*p)-alpha_span_side/2*pi/180-(jj-1)*alpha_span_side*pi/180 + (im-1)*pi/p),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1-alpha_p)*pi/(2*p)-alpha_span_side/2*pi/180-(jj-1)*alpha_span_side*pi/180 + (im-1)*pi/p));
                        mi_setblockprop(['PM',num2str(c)], 1, 0, 0 , num2str((90/p-delta_m*180/pi+180/p*(im-1))+sign(R_w-R_m)*(90-theta_m_side*180/pi)+sign(R_w-R_m)*180*(im-1)), c,1)
                        mi_clearselected
                    end
                end
                
                first_side(1,im) = c; % maximum group number of the first side-magnet over each pole
                
                
                for ii = 1:rad_dis-1
                    for jj = 1:dis_ang_mid-1
                        c = c+1;
                        mi_addmaterial(['PM',num2str(c)], mu_r, mu_r, B_r/(mu_r*mu_0), 0, 0.625)
                        mi_selectlabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1-alpha_p)*pi/(2*p)+alpha_span_mid/2*pi/180+(jj-1)*alpha_span_mid*pi/180 + (im-1)*pi/(p)),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1-alpha_p)*pi/(2*p)+alpha_span_mid/2*pi/180+(jj-1)*alpha_span_mid*pi/180 + (im-1)*pi/(p)));
                        mi_setblockprop(['PM',num2str(c)], 1, 0, 0 , num2str(180/(2*p)+(im-1)*180/p+sign(R_w-R_m)*(im-1)*180), c,1)
                        mi_clearselected
                    end
                end
            
                mid(1,im) = c;
                
                
                for ii = 1:rad_dis-1
                    for jj = 1:dis_ang_side-1
                        c = c+1;
                        mi_addmaterial(['PM',num2str(c)], mu_r, mu_r, B_r/(mu_r*mu_0), 0, 0.625)
                        mi_selectlabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1+alpha_p)*pi/(2*p)+alpha_span_side/2*pi/180+(jj-1)*alpha_span_side*pi/180 + (im-1)*pi/p),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1+alpha_p)*pi/(2*p)+alpha_span_side/2*pi/180+(jj-1)*alpha_span_side*pi/180 + (im-1)*pi/p));
                        mi_setblockprop(['PM',num2str(c)], 1, 0, 0 , num2str((90/p+delta_m*180/pi+180/p*(im-1))-sign(R_w-R_m)*(90-theta_m_side*180/pi)+sign(R_w-R_m)*180*(im-1)), c,1)
                        mi_clearselected
                    end
                end
                
                second_side(1,im) = c; % maximum group number of the first side-magnet over each pole
            if Halbach_1 == 5
                for ii = 1:rad_dis-1
                    for jj = 1:dis_ang_end-1
                        c = c+1;
                        mi_addmaterial(['PM',num2str(c)], mu_r, mu_r, B_r/(mu_r*mu_0), 0, 0.625)
                        mi_selectlabel((R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*cos((1+alpha_p1)*pi/(2*p)+alpha_span_end/2*pi/180+(jj-1)*alpha_span_end*pi/180 + (im-1)*pi/(p)),(R_m-sign(R_w-R_m)*rad_thick/2-(ii-1)*sign(R_w-R_m)*rad_thick)*sin((1+alpha_p1)*pi/(2*p)+alpha_span_end/2*pi/180+(jj-1)*alpha_span_end*pi/180 + (im-1)*pi/(p)));
                        mi_setblockprop(['PM',num2str(c)], 1, 0, 0 , num2str((90/p+delta_m1*180/pi+180/p*(im-1))-sign(R_w-R_m)*(90-theta_m_end*180/pi)+sign(R_w-R_m)*180*(im-1)), c,1)
                        mi_clearselected
                    end
                end

                second_end(1,im) = c; % maximum group number of the first side-magnet over each pole
            end
    end
                 
end

% coil 
 for iq = 1:Qsim
        mi_selectlabel((R_w+sign(R_w-R_m)*s_t/2)*cos((2*iq-1)*alpha_s/2*pi/180),(R_w+sign(R_w-R_m)*s_t/2)*sin((2*iq-1)*alpha_s/2*pi/180));
        mi_setblockprop('COPPER', 1, rad_thick/20, ['Is', num2str(iq)],0, c+iq,1)
        mi_clearselected
        
 end

% rotor backing
mi_selectlabel((R_i-(R_i-R_ie)/2)*cos(alpha/2*pi/180),(R_i-(R_i-R_ie)/2)*sin(alpha/2*pi/180));
mi_setblockprop('IRON', 1, rad_thick/20, 0,0, 0,1)
mi_clearselected
% if non-linear iron
% mi_selectlabel((R_se+sign(R_w-R_m)*0.05*R_se)*cos(alpha/2*pi/180),(R_se+sign(R_w-R_m)*0.05*R_se)*sin(alpha/2*pi/180));
% mi_setblockprop('AIR', 1, rad_thick/20, 0 ,0, 0,1)
% mi_clearselected

smartmesh(1)
% mi_saveas(['Inductor_',num2str(L),'_uH.fem'])
mi_saveas(filename)
%% post-processing 
iter_num = 20;
B_avg(c,iter_num) = 0; % matrix holding the computed average flux density 
                       % on every magnets group at every iteration
B_knee_iter(c,iter_num) = 0; % matrix holding the updated knee point flux 
                             % density for each magnets group at every
                             % iteration
B_rem_up(c,iter_num) = 0; % matrix holding the updated remanent flux 
                          % density for each magnets group at every
                          % iteration
H_knee_iter(c,iter_num) = 0; % matrix holding the updated knee point H 
                             % field for each magnets group at every
                             % iteration
k = 0;
for iter = 1:iter_num
 
if iter==1
    B_rem_up(:,iter) = B_r; % remanence at the first iteration
    B_knee_iter(:,iter) = B_int; % knee-point flux density at the first iteration
    H_knee_iter(:,iter) = H_int; % knee-point flux density at the first iteration
elseif iter>1
    B_knee_iter(:,iter) = comp_opp.*B_knee_iter(:,iter-1)+comp.*B_avg(:,iter-1); % updated knee-point flux density
    H_knee_iter(:,iter) = comp_opp.*H_knee_iter(:,iter-1)+comp.*(B_knee_iter(:,iter-1)./(mu_0*mu_rc)-H_c); % updated knee-point flux density
    B_rem_up(:,iter) = B_knee_iter(:,iter)+abs(mu_0*mu_r*H_knee_iter(:,iter)); % updated remanence
    
    for ii = 1:c
        mi_modifymaterial(['PM',num2str(ii)],3,B_rem_up(ii,iter)/(mu_0*mu_r)) 
    end

Var = max((B_rem_up(:,iter-1)-B_rem_up(:,iter))./B_rem_up(:,iter-1)*100); % maximum percentage remanence variation between consecutive iterations

if Var < 1
    iter_num = iter;
    
    if iter == 2
        disp('No demagnetization detected :)')
        return
    end
    k = k+1;
    if np==2
        mo_addcontour(R_w+sign(R_w-R_m)*s_t/2,0);
        mo_addcontour((R_w+sign(R_w-R_m)*s_t/2)*cos((alpha)*2*pi/360),(R_w+sign(R_w-R_m)*s_t/2)*sin((alpha)*2*pi/360));
        mo_bendcontour(alpha,1);

        mo_makeplot(2,2000,'myfile.txt',1);
        xf=-importdata('myfile.txt');
        x2=xf(:,2);
        N = length(x2);
        % plot(x(:,1),x2)


        [ah,bh]=fft_femm(xf);
        x1=linspace(0,2*pi,N/2); 
        figure
        hold on
        h = ((1:N/2) - 1)';
        Bgf = ah'*cos(h*x1) + bh'*sin(h*x1);
        plot(x1,Bgf)
        plot(x1,ah(np)*cos((np-1)*x1) + bh(np)*sin((np-1)*x1))

        figure
        bar(h,sqrt(ah.^2+bh.^2))
        axis([0 10 -0.5 1.5])

        Bg1(k) = sqrt(ah(2).^2+bh(2).^2);
    elseif np==1
        mo_addcontour(R_w+sign(R_w-R_m)*s_t/2,0);
        mo_addcontour((R_w+sign(R_w-R_m)*s_t/2)*cos((alpha)*2*pi/360),(R_w+sign(R_w-R_m)*s_t/2)*sin((alpha)*2*pi/360));
        mo_bendcontour(alpha,1);

        mo_makeplot(2,2000,'myfile.txt',1);
        xf=-importdata('myfile.txt');
%         xf = [xf ; -xf(2:size(xf,1),:)];
        xf = [xf ; -flip(xf(1:size(xf,1)-1,:))];
        x2= xf(:,2);
        N = length(x2);
        % plot(x(:,1),x2)


        [ah,bh]=fft_femm(xf);
        x1=linspace(0,2*pi,N/2); 
        figure
        hold on
        h = ((1:N/2) - 1)';
        Bgf = ah'*cos(h*x1) + bh'*sin(h*x1);
        plot(x1,Bgf)
        plot(x1,ah(np)*cos((np-1)*x1) + bh(np)*sin((np-1)*x1))

        figure
        bar(h,sqrt(ah.^2+bh.^2))
        axis([0 10 -0.5 1.5])

        Bg1(k) = sqrt(ah(2).^2+bh(2).^2);
    end
        break
end

end
    

mi_analyze(1); % analyze the model
mi_loadsolution;

%% POST PROCESSING (DEMAGNETIZATION ANALYSIS)
if (iter == 1)
        % Record the initial mesh elements for the first iteration 
        nn = mo_numelements; % number of mesh elements
        b = zeros(1,nn); % matrix that will hold the flux density info
        z = zeros(nn,1); % Location of the centroid of each element
        a = zeros(nn,1); % Area of each element
        g = zeros(nn,1); % Block label of each element
        for m = 1:nn
            elm = mo_getelement(m);
            % z is a vector of complex numbers that represents the location of
            % the centroid of each element.
            z(m) = elm(4) + 1j*elm(5);
            % element area in the length units used to draw the geometry
            a(m) = elm(6);
            % group number associated with the element
            g(m) = elm(7); 
        end
end

% Store element flux densities *)

% Special Halbach demagnetization analysis
if alpha_p == 0
    if np == 1
        b_tot(c,nn) = 0;
        
        for m = 1:nn
        
            if g(m) > 0 && g(m) <= c
                % Element is on the first or second side magnet (first pole)
                % Store flux density at the element centroid for these elements			
                b_tot(g(m),m) = (mo_getb(real(z(m)),imag(z(m)))*[1;1j]);
                
            end
                  
        end
        
        % with b_tot_ones we put ones on every non-zero entry of b_tot
        b_tot_ones = (b_tot~= 0);
        Areas = a'.*b_tot_ones;
        % Projection along the magnetization direction for each magnet
        b_first_side_1_mag = real(b_tot(1:first_side(1,1),:))*cos(pi/(2*p)-delta_m+sign(R_w-R_m)*(pi/2-theta_m_side))+imag(b_tot(1:first_side(1,1),:))*sin(pi/(2*p)-delta_m+sign(R_w-R_m)*(pi/2-theta_m_side));
        b_second_side_1_mag = real(b_tot(first_side(1,1)+1:second_side(1,1),:))*cos(pi/(2*p)+delta_m-sign(R_w-R_m)*(pi/2-theta_m_side))+imag(b_tot(first_side(1,1)+1:second_side(1,1),:))*sin(pi/(2*p)+delta_m-sign(R_w-R_m)*(pi/2-theta_m_side));
        b_mag = [b_first_side_1_mag ; b_second_side_1_mag];
        
        % elements within a single group have different size. The flux
        % density is weighted with the respective element size to evaluate
        % the average flux density level within each group
        weighted_B = sum(b_mag.*Areas,2)./(sum(Areas,2));
        
        if iter==1    
            Enabler1 = ([b_first_side_1_mag]<B_knee_iter(1,1)) .*([b_first_side_1_mag ]~= 0);
            Demag_portion1 = sum(sum(Areas(1:first_side(1,1),:).*Enabler1,1),2);
            Perc_demag1 = Demag_portion1/sum(sum(Areas(1:first_side(1,1),:),1),2)*100;
            Enabler2 = ([b_second_side_1_mag]<B_knee_iter(1,1)) .*([b_second_side_1_mag ]~= 0);
            Demag_portion2 = sum(sum(Areas(first_side(1,1)+1:c,:).*Enabler2,1),2);
            Perc_demag2 = Demag_portion2/sum(sum(Areas(first_side(1,1)+1:c,:),1),2)*100;
            Enabler = ([b_first_side_1_mag ; b_second_side_1_mag]<B_knee_iter(1,1)) .*([b_first_side_1_mag ; b_second_side_1_mag]~= 0);
            Demag_portion = sum(sum(Areas.*Enabler,1),2);
            Perc_demag = Demag_portion/sum(sum(Areas,1),2)*100;
            % Normalized distance to the knee-point flux density
            Weights = (B_knee_iter(1,1)-[b_first_side_1_mag ; b_second_side_1_mag]).*Enabler./B_knee_iter(1,1);
            Weighted_Perc_demag = sum(sum(Areas.*Weights,2),1)/sum(sum(Areas,1),2)*100;
            B_dem_M = (B_knee_iter(1,1) - b_mag).*(b_mag~= 0);
            B_dem_M(B_dem_M <= 0) = 0;
            B_dem = sum(B_dem_M,1);
        end
        
    elseif np == 2
        % b_tot has size (total number os magnets divisions,total number of
        % mesh elements) and it's a highly sparse matrix
        b_tot(c,nn) = 0;
        
        tic
        for m = 1:nn
        
            if g(m) > 0 && g(m) <= c
                % Element is on the first or second side magnet (first and second pole)
                % Store flux density at the element centroid for these elements			
                b_tot(g(m),m) = (mo_getb(real(z(m)),imag(z(m)))*[1;1j]);
                
            end
                  
        end
        toc
        % with b_tot_ones we put ones on every non-zero entry of b_tot
        b_tot_ones = (b_tot~= 0);
        Areas = a'.*b_tot_ones;
        % Projection along the magnetization direction for each magnet
        b_first_side_1_mag = real(b_tot(1:first_side(1,1),:))*cos(pi/(2*p)-delta_m+sign(R_w-R_m)*(pi/2-theta_m_side))+imag(b_tot(1:first_side(1,1),:))*sin(pi/(2*p)-delta_m+sign(R_w-R_m)*(pi/2-theta_m_side));
        b_second_side_1_mag = real(b_tot(first_side(1,1)+1:second_side(1,1),:))*cos(pi/(2*p)+delta_m-sign(R_w-R_m)*(pi/2-theta_m_side))+imag(b_tot(first_side(1,1)+1:second_side(1,1),:))*sin(pi/(2*p)+delta_m-sign(R_w-R_m)*(pi/2-theta_m_side));
        b_first_side_2_mag = real(b_tot(second_side(1,1)+1:first_side(1,2),:))*cos(pi/(2*p)-delta_m+pi/p+sign(R_w-R_m)*(pi/2-theta_m_side)+pi)+imag(b_tot(second_side(1,1)+1:first_side(1,2),:))*sin(pi/(2*p)-delta_m+pi/p+sign(R_w-R_m)*(pi/2-theta_m_side)+pi);
        b_second_side_2_mag = real(b_tot(first_side(1,2)+1:second_side(1,2),:))*cos(pi/(2*p)+delta_m+pi/p-sign(R_w-R_m)*(pi/2-theta_m_side)+pi)+imag(b_tot(first_side(1,1)+1:second_side(1,1),:))*sin(pi/(2*p)+delta_m+pi/p-sign(R_w-R_m)*(pi/2-theta_m_side)+pi);

        % elements within a single group have different size. The flux
        % density is weighted with the respective element size to evaluate
        % the average flux density level within each group
        weighted_B = sum([b_first_side_1_mag ; b_second_side_1_mag ; b_first_side_2_mag ; b_second_side_2_mag].*Areas,2)./(sum(Areas,2));
        
    end  
end


% Two-segments Halbach array    
if Halbach_1 == 2 && alpha_p~= 0 % Even-segments Halbach (2-segments)
    if np == 1
        b_tot(c,nn) = 0;
        
        for m = 1:nn
        
            if g(m) > 0 && g(m) <= c
                % Element is on the first or second side magnet (first pole)
                % Store flux density at the element centroid for these elements			
                b_tot(g(m),m) = (mo_getb(real(z(m)),imag(z(m)))*[1;1j]);
                
            end
                  
        end
        
        % with b_tot_ones we put ones on every non-zero entry of b_tot
        b_tot_ones = (b_tot~= 0);
        Areas = a'.*b_tot_ones;
        % Projection along the magnetization direction for each magnet
        b_first_end_1_mag = real(b_tot(1:first_end(1,1),:))*cos(pi/2)+imag(b_tot(1:first_end(1,1),:))*sin(pi/2);
        b_mid_1_mag = real(b_tot(first_end(1,1)+1:mid(1,1),:))*cos(pi/(2*p))+imag(b_tot(first_end(1,1)+1:mid(1,1),:))*sin(pi/(2*p));
        b_second_end_1_mag = real(b_tot(mid(1,1)+1:second_end(1,1),:))*cos(pi/p-pi/2)+imag(b_tot(mid(1,1)+1:second_end(1,1),:))*sin(pi/p-pi/2);
        b_mag = [b_first_end_1_mag ; b_mid_1_mag ; b_second_end_1_mag];
        
        % elements within a single group have different size. The flux
        % density is weighted with the respective element size to evaluate
        % the average flux density level within each group
        weighted_B = sum(b_mag.*Areas,2)./(sum(Areas,2));
        
        if iter==1     
            Enabler1 = ([b_first_end_1_mag]<B_knee_iter(1,1)) .*([b_first_end_1_mag]~= 0);
            Demag_portion1 = sum(sum(Areas(1:first_end(1,1),:).*Enabler1,1),2);
            Perc_demag1 = Demag_portion1/sum(sum(Areas(1:first_end(1,1),:),1),2)*100;
            Enabler2 = ([b_mid_1_mag]<B_knee_iter(1,1)) .*([b_mid_1_mag]~= 0);
            Demag_portion2 = sum(sum(Areas(first_end(1,1)+1:mid(1,1),:).*Enabler2,1),2);
            Perc_demag2 = Demag_portion2/sum(sum(Areas(first_end(1,1)+1:mid(1,1),:),1),2)*100;
            Enabler3 = ([b_second_end_1_mag]<B_knee_iter(1,1)) .*([b_second_end_1_mag]~= 0);
            Demag_portion3 = sum(sum(Areas(mid(1,1)+1:second_end(1,1),:).*Enabler3,1),2);
            Perc_demag3 = Demag_portion3/sum(sum(Areas(mid(1,1)+1:second_end(1,1),:),1),2)*100;
            
            Enabler = ([b_first_end_1_mag ; b_mid_1_mag ; b_second_end_1_mag]<B_knee_iter(1,1)) .*([b_first_end_1_mag ; b_mid_1_mag ; b_second_end_1_mag]~= 0);
            Demag_portion = sum(sum(Areas.*Enabler,1),2);
            Perc_demag = Demag_portion/sum(sum(Areas,1),2)*100;
           % Normalized distance to the knee-point flux density
            Weights = (B_knee_iter(1,1)-[b_first_end_1_mag ; b_mid_1_mag ; b_second_end_1_mag]).*Enabler./B_knee_iter(1,1);
            Weighted_Perc_demag = sum(sum(Areas.*Weights,2),1)/sum(sum(Areas,1),2)*100;
            B_dem_M = (B_knee_iter(1,1) - b_mag).*(b_mag~= 0);
            B_dem_M(B_dem_M <= 0) = 0;
            B_dem = sum(B_dem_M,1);
        end
        
    elseif np == 2
        % b_tot has size (total number os magnets divisions,total number of
        % mesh elements) and it's a highly sparse matrix
        b_tot(c,nn) = 0;
        
        
        for m = 1:nn
        
            if g(m) > 0 && g(m) <= c
                % Element is on the first or second side magnet (first pole)
                % Store flux density at the element centroid for these elements			
                b_tot(g(m),m) = (mo_getb(real(z(m)),imag(z(m)))*[1;1j]);
                
            end
                  
        end
        
        % with b_tot_ones we put ones on every non-zero entry of b_tot
        b_tot_ones = (b_tot~= 0);
        Areas = a'.*b_tot_ones;
        % Projection along the magnetization direction for each magnet
        b_first_end_1_mag = real(b_tot(1:first_end(1,1),:))*cos(pi/2)+imag(b_tot(1:first_end(1,1),:))*sin(pi/2);
        b_mid_1_mag = real(b_tot(first_end(1,1)+1:mid(1,1),:))*cos(pi/(2*p))+imag(b_tot(first_end(1,1)+1:mid(1,1),:))*sin(pi/(2*p));
        b_second_end_1_mag = real(b_tot(mid(1,1)+1:second_end(1,1),:))*cos(pi/p-pi/2)+imag(b_tot(mid(1,1)+1:second_end(1,1),:))*sin(pi/p-pi/2);
        
        b_first_end_2_mag = real(b_tot(second_end(1,1)+1:first_end(1,2),:))*cos(pi/p-pi/2)+imag(b_tot(second_end(1,1)+1:first_end(1,2),:))*sin(pi/p-pi/2);
        b_mid_2_mag = real(b_tot(first_end(1,2)+1:mid(1,2),:))*cos(pi/(2*p)+pi/p+pi)+imag(b_tot(first_end(1,2)+1:mid(1,2),:))*sin(pi/(2*p)+pi/p+pi);
        b_second_end_2_mag = real(b_tot(mid(1,2)+1:second_end(1,2),:))*cos(pi/p-pi/2+pi/p+pi)+imag(b_tot(mid(1,2)+1:second_end(1,2),:))*sin(pi/p-pi/2+pi/p+pi);
        
        % elements within a single group have different size. The flux
        % density is weighted with the respective element size to evaluate
        % the average flux density level within each group
        weighted_B = sum([b_first_end_1_mag ; b_mid_1_mag ; b_second_end_1_mag ; b_first_end_2_mag ; b_mid_2_mag ; b_second_end_2_mag].*Areas,2)./(sum(Areas,2));
        
    end  
   
    
end


% Three-segments Halbach array    
if Halbach_1 == 3 % Even-segments Halbach (4-segments)
    if np == 1
        b_tot(c,nn) = 0;
        
        for m = 1:nn
        
            if g(m) > 0 && g(m) <= c
                % Element is on the first or second side magnet (first pole)
                % Store flux density at the element centroid for these elements			
                b_tot(g(m),m) = (mo_getb(real(z(m)),imag(z(m)))*[1;1j]);
                
            end
                  
        end
        
        % with b_tot_ones we put ones on every non-zero entry of b_tot
        b_tot_ones = (b_tot~= 0);
        Areas = a'.*b_tot_ones;
        % Projection along the magnetization direction for each magnet
        b_first_side_1_mag = real(b_tot(1:first_side(1,1),:))*cos(pi/(2*p)-delta_m+sign(R_w-R_m)*(pi/2-theta_m_side))+imag(b_tot(1:first_side(1,1),:))*sin(pi/(2*p)-delta_m+sign(R_w-R_m)*(pi/2-theta_m_side));
        b_mid_1_mag = real(b_tot(first_side(1,1)+1:mid(1,1),:))*cos(pi/(2*p))+imag(b_tot(first_side(1,1)+1:mid(1,1),:))*sin(pi/(2*p));
        b_second_side_1_mag = real(b_tot(mid(1,1)+1:second_side(1,1),:))*cos(pi/(2*p)+delta_m-sign(R_w-R_m)*(pi/2-theta_m_side))+imag(b_tot(mid(1,1)+1:second_side(1,1),:))*sin(pi/(2*p)+delta_m-sign(R_w-R_m)*(pi/2-theta_m_side));
        b_mag = [b_first_side_1_mag ; b_mid_1_mag ; b_second_side_1_mag];
        
        % elements within a single group have different size. The flux
        % density is weighted with the respective element size to evaluate
        % the average flux density level within each group
        weighted_B = sum(b_mag.*Areas,2)./(sum(Areas,2));
        
        
        if iter==1         
            
            Enabler1 = ([b_first_side_1_mag]<B_knee_iter(1,1)) .*([b_first_side_1_mag]~= 0);
            Demag_portion1 = sum(sum(Areas(1:first_side(1,1),:).*Enabler1,1),2);
            Perc_demag1 = Demag_portion1/sum(sum(Areas(1:first_side(1,1),:),1),2)*100;
            Enabler2 = ([b_mid_1_mag]<B_knee_iter(1,1)) .*([b_mid_1_mag]~= 0);
            Demag_portion2 = sum(sum(Areas(first_side(1,1)+1:mid(1,1),:).*Enabler2,1),2);
            Perc_demag2 = Demag_portion2/sum(sum(Areas(first_side(1,1)+1:mid(1,1),:),1),2)*100;
            Enabler3 = ([b_second_side_1_mag]<B_knee_iter(1,1)) .*([b_second_side_1_mag]~= 0);
            Demag_portion3 = sum(sum(Areas(mid(1,1)+1:second_side(1,1),:).*Enabler3,1),2);
            Perc_demag3 = Demag_portion3/sum(sum(Areas(mid(1,1)+1:second_side(1,1),:),1),2)*100;
            
            
            Enabler = ([b_first_side_1_mag ; b_mid_1_mag ; b_second_side_1_mag]<B_knee_iter(1,1)) .*([b_first_side_1_mag ; b_mid_1_mag ; b_second_side_1_mag]~= 0);
            Demag_portion = sum(sum(Areas.*Enabler,1),2);
            Perc_demag = Demag_portion/sum(sum(Areas,1),2)*100;
            % Normalized distance to the knee-point flux density
            Weights = (B_knee_iter(1,1)-[b_first_side_1_mag ; b_mid_1_mag ; b_second_side_1_mag]).*Enabler./B_knee_iter(1,1);
            Weighted_Perc_demag = sum(sum(Areas.*Weights,2),1)/sum(sum(Areas,1),2)*100;
            B_dem_M = (B_knee_iter(1,1) - b_mag).*(b_mag~= 0);
            B_dem_M(B_dem_M <= 0) = 0;
            B_dem = sum(B_dem_M,1);
         end
        
    elseif np == 2
        % b_tot has size (total number os magnets divisions,total number of
        % mesh elements) and it's a highly sparse matrix
        b_tot(c,nn) = 0;
        
        
        for m = 1:nn
        
            if g(m) > 0 && g(m) <= c
                % Element is on the first or second side magnet (first pole)
                % Store flux density at the element centroid for these elements			
                b_tot(g(m),m) = (mo_getb(real(z(m)),imag(z(m)))*[1;1j]);
                
            end
                  
        end
        
        % with b_tot_ones we put ones on every non-zero entry of b_tot
        b_tot_ones = (b_tot~= 0);
        Areas = a'.*b_tot_ones;
        % Projection along the magnetization direction for each magnet
        b_first_side_1_mag = real(b_tot(1:first_side(1,1),:))*cos(pi/(2*p)-delta_m+sign(R_w-R_m)*(pi/2-theta_m_side))+imag(b_tot(1:first_side(1,1),:))*sin(pi/(2*p)-delta_m+sign(R_w-R_m)*(pi/2-theta_m_side));
        b_mid_1_mag = real(b_tot(first_side(1,1)+1:mid(1,1),:))*cos(pi/(2*p))+imag(b_tot(first_side(1,1)+1:mid(1,1),:))*sin(pi/(2*p));
        b_second_side_1_mag = real(b_tot(mid(1,1)+1:second_side(1,1),:))*cos(pi/(2*p)+delta_m-sign(R_w-R_m)*(pi/2-theta_m_side))+imag(b_tot(mid(1,1)+1:second_side(1,1),:))*sin(pi/(2*p)+delta_m-sign(R_w-R_m)*(pi/2-theta_m_side));
        
        b_first_side_2_mag = real(b_tot(second_side(1,1)+1:first_side(1,2),:))*cos(pi/(2*p)-delta_m+sign(R_w-R_m)*(pi/2-theta_m_side+pi/p+pi))+imag(b_tot(second_side(1,1)+1:first_side(1,2),:))*sin(pi/(2*p)-delta_m+sign(R_w-R_m)*(pi/2-theta_m_side)+pi/p+pi);
        b_mid_2_mag = real(b_tot(first_side(1,2)+1:mid(1,2),:))*cos(pi/(2*p)+pi/p+pi)+imag(b_tot(first_side(1,2)+1:mid(1,2),:))*sin(pi/(2*p)+pi/p+pi);
        b_second_side_2_mag = real(b_tot(mid(1,2)+1:second_side(1,2),:))*cos(pi/(2*p)+delta_m-sign(R_w-R_m)*(pi/2-theta_m_side)+pi/p+pi)+imag(b_tot(mid(1,2)+1:second_side(1,2),:))*sin(pi/(2*p)+delta_m-sign(R_w-R_m)*(pi/2-theta_m_side)+pi/p+pi);
        
        % elements within a single group have different size. The flux
        % density is weighted with the respective element size to evaluate
        % the average flux density level within each group
        weighted_B = sum([b_first_side_1_mag ; b_mid_1_mag ; b_second_side_1_mag ; b_first_side_2_mag ; b_mid_2_mag ; b_second_side_2_mag].*Areas,2)./(sum(Areas,2));
        
    end 
    
end

% Four-segments Halbach array    
if Halbach_1 == 4 % Even-segments Halbach (4-segments)
    if np == 1
        b_tot(c,nn) = 0;
        
        for m = 1:nn
        
            if g(m) > 0 && g(m) <= c
                % Element is on the first or second side magnet (first pole)
                % Store flux density at the element centroid for these elements			
                b_tot(g(m),m) = (mo_getb(real(z(m)),imag(z(m)))*[1;1j]);
                
            end
                  
        end
        
        % with b_tot_ones we put ones on every non-zero entry of b_tot
        b_tot_ones = (b_tot~= 0);
        Areas = a'.*b_tot_ones;
        % Projection along the magnetization direction for each magnet
        b_first_end_1_mag = real(b_tot(1:first_end(1,1),:))*cos(pi/2)+imag(b_tot(1:first_end(1,1),:))*sin(pi/2);
        b_first_side_1_mag = real(b_tot(first_end(1,1)+1:first_side(1,1),:))*cos(pi/(2*p)-delta_m+sign(R_w-R_m)*(pi/2-theta_m_side))+imag(b_tot(first_end(1,1)+1:first_side(1,1),:))*sin(pi/(2*p)-delta_m+sign(R_w-R_m)*(pi/2-theta_m_side));
        b_mid_1_mag = real(b_tot(first_side(1,1)+1:mid(1,1),:))*cos(pi/(2*p))+imag(b_tot(first_side(1,1)+1:mid(1,1),:))*sin(pi/(2*p));
        b_second_side_1_mag = real(b_tot(mid(1,1)+1:second_side(1,1),:))*cos(pi/(2*p)+delta_m-sign(R_w-R_m)*(pi/2-theta_m_side))+imag(b_tot(mid(1,1)+1:second_side(1,1),:))*sin(pi/(2*p)+delta_m-sign(R_w-R_m)*(pi/2-theta_m_side));
        b_second_end_1_mag = real(b_tot(second_side(1,1)+1:second_end(1,1),:))*cos(pi/p-pi/2)+imag(b_tot(second_side(1,1)+1:second_end(1,1),:))*sin(pi/p-pi/2);
        b_mag = [b_first_end_1_mag ; b_first_side_1_mag ; b_mid_1_mag ; b_second_side_1_mag ; b_second_end_1_mag];
        
        % elements within a single group have different size. The flux
        % density is weighted with the respective element size to evaluate
        % the average flux density level within each group
        weighted_B = sum(b_mag.*Areas,2)./(sum(Areas,2));
       
        
        if iter==1           
            Enabler = [b_first_end_1_mag ; b_first_side_1_mag ; b_mid_1_mag ; b_second_side_1_mag ; b_second_end_1_mag]<B_knee_iter(1,1) *([b_first_end_1_mag ; b_first_side_1_mag ; b_mid_1_mag ; b_second_side_1_mag ; b_second_end_1_mag]~= 0);
            Demag_portion = sum(sum(Areas*Enabler',1),2);
            Perc_demag = Demag_portion/sum(sum(Areas,1),2)*100;
            % Normalized distance to the knee-point flux density
            Weights = (B_knee_iter(1,1)-[b_first_end_1_mag ; b_first_side_1_mag ; b_mid_1_mag ; b_second_side_1_mag ; b_second_end_1_mag]).*Enabler./B_knee_iter(1,1);
            Weighted_Perc_demag = sum(sum(Areas.*Weights,2),1)/sum(sum(Areas,1),2)*100;
            B_dem_M = (B_knee_iter(1,1) - b_mag).*(b_mag~= 0);
            B_dem_M(B_dem_M <= 0) = 0;
            B_dem = sum(B_dem_M,1);
        end
        
    elseif np == 2
        % b_tot has size (total number os magnets divisions,total number of
        % mesh elements) and it's a highly sparse matrix
        b_tot(c,nn) = 0;
        
        
        for m = 1:nn
        
            if g(m) > 0 && g(m) <= c
                % Element is on the first or second side magnet (first pole)
                % Store flux density at the element centroid for these elements			
                b_tot(g(m),m) = (mo_getb(real(z(m)),imag(z(m)))*[1;1j]);
                
            end
                  
        end
        
        % with b_tot_ones we put ones on every non-zero entry of b_tot
        b_tot_ones = (b_tot~= 0);
        Areas = a'.*b_tot_ones;
        % Projection along the magnetization direction for each magnet
        b_first_end_1_mag = real(b_tot(1:first_end(1,1),:))*cos(pi/2)+imag(b_tot(1:first_end(1,1),:))*sin(pi/2);
        b_first_side_1_mag = real(b_tot(first_end(1,1)+1:first_side(1,1),:))*cos(pi/(2*p)-delta_m+sign(R_w-R_m)*(pi/2-theta_m_side))+imag(b_tot(first_end(1,1)+1:first_side(1,1),:))*sin(pi/(2*p)-delta_m+sign(R_w-R_m)*(pi/2-theta_m_side));
        b_mid_1_mag = real(b_tot(first_side(1,1)+1:mid(1,1),:))*cos(pi/(2*p))+imag(b_tot(first_side(1,1)+1:mid(1,1),:))*sin(pi/(2*p));
        b_second_side_1_mag = real(b_tot(mid(1,1)+1:second_side(1,1),:))*cos(pi/(2*p)+delta_m-sign(R_w-R_m)*(pi/2-theta_m_side))+imag(b_tot(mid(1,1)+1:second_side(1,1),:))*sin(pi/(2*p)+delta_m-sign(R_w-R_m)*(pi/2-theta_m_side));
        b_second_end_1_mag = real(b_tot(second_side(1,1)+1:second_end(1,1),:))*cos(pi/p-pi/2)+imag(b_tot(second_side(1,1)+1:second_end(1,1),:))*sin(pi/p-pi/2);
        
        b_first_end_2_mag = real(b_tot(second_end(1,1)+1:first_end(1,2),:))*cos(pi/p-pi/2)+imag(b_tot(second_end(1,1)+1:first_end(1,2),:))*sin(pi/p-pi/2);
        b_first_side_2_mag = real(b_tot(first_end(1,2)+1:first_side(1,2),:))*cos(pi/(2*p)-delta_m+sign(R_w-R_m)*(pi/2-theta_m_side+pi/p+pi))+imag(b_tot(first_end(1,2)+1:first_side(1,2),:))*sin(pi/(2*p)-delta_m+sign(R_w-R_m)*(pi/2-theta_m_side)+pi/p+pi);
        b_mid_2_mag = real(b_tot(first_side(1,2)+1:mid(1,2),:))*cos(pi/(2*p)+pi/p+pi)+imag(b_tot(first_side(1,2)+1:mid(1,2),:))*sin(pi/(2*p)+pi/p+pi);
        b_second_side_2_mag = real(b_tot(mid(1,2)+1:second_side(1,2),:))*cos(pi/(2*p)+delta_m-sign(R_w-R_m)*(pi/2-theta_m_side)+pi/p+pi)+imag(b_tot(mid(1,2)+1:second_side(1,2),:))*sin(pi/(2*p)+delta_m-sign(R_w-R_m)*(pi/2-theta_m_side)+pi/p+pi);
        b_second_end_2_mag = real(b_tot(second_side(1,2)+1:second_end(1,2),:))*cos(pi/p-pi/2+pi/p+pi)+imag(b_tot(second_side(1,2)+1:second_end(1,2),:))*sin(pi/p-pi/2+pi/p+pi);
        
        % elements within a single group have different size. The flux
        % density is weighted with the respective element size to evaluate
        % the average flux density level within each group
        weighted_B = sum([b_first_end_1_mag ; b_first_side_1_mag ; b_mid_1_mag ; b_second_side_1_mag ; b_second_end_1_mag ; b_first_end_2_mag ; b_first_side_2_mag ; b_mid_2_mag ; b_second_side_2_mag ; b_second_end_2_mag].*Areas,2)./(sum(Areas,2));
        
        
        
    end 
    
end

% Five-segments Halbach array    
if Halbach_1 == 5 % Odd-segments Halbach (5-segments)
    if np == 1
        b_tot(c,nn) = 0;
        
        for m = 1:nn
        
            if g(m) > 0 && g(m) <= c
                % Element is on the first or second side magnet (first pole)
                % Store flux density at the element centroid for these elements			
                b_tot(g(m),m) = (mo_getb(real(z(m)),imag(z(m)))*[1;1j]);
                
            end
                  
        end
        
        % with b_tot_ones we put ones on every non-zero entry of b_tot
        b_tot_ones = (b_tot~= 0);
        Areas = a'.*b_tot_ones;
        % Projection along the magnetization direction for each magnet
        b_first_end_1_mag = real(b_tot(1:first_end(1,1),:))*cos(pi/(2*p)-delta_m1+sign(R_w-R_m)*(pi/2-theta_m_end))+imag(b_tot(1:first_end(1,1),:))*sin(pi/(2*p)-delta_m1+sign(R_w-R_m)*(pi/2-theta_m_end));
        b_first_side_1_mag = real(b_tot(first_end(1,1)+1:first_side(1,1),:))*cos(pi/(2*p)-delta_m+sign(R_w-R_m)*(pi/2-theta_m_side))+imag(b_tot(first_end(1,1)+1:first_side(1,1),:))*sin(pi/(2*p)-delta_m+sign(R_w-R_m)*(pi/2-theta_m_side));
        b_mid_1_mag = real(b_tot(first_side(1,1)+1:mid(1,1),:))*cos(pi/(2*p))+imag(b_tot(first_side(1,1)+1:mid(1,1),:))*sin(pi/(2*p));
        b_second_side_1_mag = real(b_tot(mid(1,1)+1:second_side(1,1),:))*cos(pi/(2*p)+delta_m-sign(R_w-R_m)*(pi/2-theta_m_side))+imag(b_tot(mid(1,1)+1:second_side(1,1),:))*sin(pi/(2*p)+delta_m-sign(R_w-R_m)*(pi/2-theta_m_side));
        b_second_end_1_mag = real(b_tot(second_side(1,1)+1:second_end(1,1),:))*cos(pi/(2*p)+delta_m1-sign(R_w-R_m)*(pi/2-theta_m_end))+imag(b_tot(second_side(1,1)+1:second_end(1,1),:))*sin(pi/(2*p)+delta_m1-sign(R_w-R_m)*(pi/2-theta_m_end));
        b_mag = [b_first_end_1_mag ; b_first_side_1_mag ; b_mid_1_mag ; b_second_side_1_mag ; b_second_end_1_mag];
        
        % elements within a single group have different size. The flux
        % density is weighted with the respective element size to evaluate
        % the average flux density level within each group
        weighted_B = sum(b_mag.*Areas,2)./(sum(Areas,2));
        
        
        if iter==1           
            Enabler = [b_first_end_1_mag ; b_first_side_1_mag ; b_mid_1_mag ; b_second_side_1_mag ; b_second_end_1_mag]<B_knee_iter(1,1) *([b_first_end_1_mag ; b_first_side_1_mag ; b_mid_1_mag ; b_second_side_1_mag ; b_second_end_1_mag]~= 0);
            Demag_portion = sum(sum(Areas*Enabler',1),2);
            Perc_demag = Demag_portion/sum(sum(Areas,1),2)*100;
            % Normalized distance to the knee-point flux density
            Weights = (B_knee_iter(1,1)-[b_first_end_1_mag ; b_first_side_1_mag ; b_mid_1_mag ; b_second_side_1_mag ; b_second_end_1_mag]).*Enabler./B_knee_iter(1,1);
            Weighted_Perc_demag = sum(sum(Areas.*Weights,2),1)/sum(sum(Areas,1),2)*100;
            B_dem_M = (B_knee_iter(1,1) - b_mag).*(b_mag~= 0);
            B_dem_M(B_dem_M <= 0) = 0;
            B_dem = sum(B_dem_M,1);
        end
        
        
    elseif np == 2
        % b_tot has size (total number os magnets divisions,total number of
        % mesh elements) and it's a highly sparse matrix
        b_tot(c,nn) = 0;
        
        
        for m = 1:nn
        
            if g(m) > 0 && g(m) <= c
                % Element is on the first or second side magnet (first pole)
                % Store flux density at the element centroid for these elements			
                b_tot(g(m),m) = (mo_getb(real(z(m)),imag(z(m)))*[1;1j]);
                
            end
                  
        end
        
        % with b_tot_ones we put ones on every non-zero entry of b_tot
        b_tot_ones = (b_tot~= 0);
        Areas = a'.*b_tot_ones;
        % Projection along the magnetization direction for each magnet
        b_first_end_1_mag = real(b_tot(1:first_end(1,1),:))*cos(pi/(2*p)-delta_m1+sign(R_w-R_m)*(pi/2-theta_m_end))+imag(b_tot(1:first_end(1,1),:))*sin(pi/(2*p)-delta_m1+sign(R_w-R_m)*(pi/2-theta_m_end));
        b_first_side_1_mag = real(b_tot(first_end(1,1)+1:first_side(1,1),:))*cos(pi/(2*p)-delta_m+sign(R_w-R_m)*(pi/2-theta_m_side))+imag(b_tot(first_end(1,1)+1:first_side(1,1),:))*sin(pi/(2*p)-delta_m+sign(R_w-R_m)*(pi/2-theta_m_side));
        b_mid_1_mag = real(b_tot(first_side(1,1)+1:mid(1,1),:))*cos(pi/(2*p))+imag(b_tot(first_side(1,1)+1:mid(1,1),:))*sin(pi/(2*p));
        b_second_side_1_mag = real(b_tot(mid(1,1)+1:second_side(1,1),:))*cos(pi/(2*p)+delta_m-sign(R_w-R_m)*(pi/2-theta_m_side))+imag(b_tot(mid(1,1)+1:second_side(1,1),:))*sin(pi/(2*p)+delta_m-sign(R_w-R_m)*(pi/2-theta_m_side));
        b_second_end_1_mag = real(b_tot(second_side(1,1)+1:second_end(1,1),:))*cos(pi/(2*p)+delta_m1-sign(R_w-R_m)*(pi/2-theta_m_end))+imag(b_tot(second_side(1,1)+1:second_end(1,1),:))*sin(pi/(2*p)+delta_m1-sign(R_w-R_m)*(pi/2-theta_m_end));
        
        b_first_end_2_mag = real(b_tot(second_end(1,1)+1:first_end(1,2),:))*cos(pi/(2*p)-delta_m1+sign(R_w-R_m)*(pi/2-theta_m_end)+pi/p+pi)+imag(b_tot(second_end(1,1)+1:first_end(1,2),:))*sin(pi/(2*p)-delta_m1+sign(R_w-R_m)*(pi/2-theta_m_end)+pi/p+pi);
        b_first_side_2_mag = real(b_tot(first_end(1,2)+1:first_side(1,2),:))*cos(pi/(2*p)-delta_m+sign(R_w-R_m)*(pi/2-theta_m_side+pi/p+pi))+imag(b_tot(first_end(1,2)+1:first_side(1,2),:))*sin(pi/(2*p)-delta_m+sign(R_w-R_m)*(pi/2-theta_m_side)+pi/p+pi);
        b_mid_2_mag = real(b_tot(first_side(1,2)+1:mid(1,2),:))*cos(pi/(2*p)+pi/p+pi)+imag(b_tot(first_side(1,2)+1:mid(1,2),:))*sin(pi/(2*p)+pi/p+pi);
        b_second_side_2_mag = real(b_tot(mid(1,2)+1:second_side(1,2),:))*cos(pi/(2*p)+delta_m-sign(R_w-R_m)*(pi/2-theta_m_side)+pi/p+pi)+imag(b_tot(mid(1,2)+1:second_side(1,2),:))*sin(pi/(2*p)+delta_m-sign(R_w-R_m)*(pi/2-theta_m_side)+pi/p+pi);
        b_second_end_2_mag = real(b_tot(second_side(1,2)+1:second_end(1,2),:))*cos(pi/(2*p)+delta_m1-sign(R_w-R_m)*(pi/2-theta_m_end)+pi/p+pi)+imag(b_tot(second_side(1,1)+1:second_end(1,2),:))*sin(pi/(2*p)+delta_m1-sign(R_w-R_m)*(pi/2-theta_m_end)+pi/p+pi);
        
        % elements within a single group have different size. The flux
        % density is weighted with the respective element size to evaluate
        % the average flux density level within each group
        weighted_B = sum([b_first_end_1_mag ; b_first_side_1_mag ; b_mid_1_mag ; b_second_side_1_mag ; b_second_end_1_mag ; b_first_end_2_mag ; b_first_side_2_mag ; b_mid_2_mag ; b_second_side_2_mag ; b_second_end_2_mag].*Areas,2)./(sum(Areas,2));
        
        
    end 
    
end

    B_avg(:,iter) = weighted_B;

    % the two following variables are useful to update the knee point flux
    % density for only those groups where the flux density drops below the
    % knee point flux density
    comp = B_avg(:,iter)<B_knee_iter(:,iter);
    comp_opp = ~comp; % groups entries where the detected flux density is below the knee point flux density are marked with zero



if iter == 1 
    if np==2
        k = k+1;
        mo_addcontour(R_w+sign(R_w-R_m)*s_t/2,0);
        mo_addcontour((R_w+sign(R_w-R_m)*s_t/2)*cos((alpha)*2*pi/360),(R_w+sign(R_w-R_m)*s_t/2)*sin((alpha)*2*pi/360));
        mo_bendcontour(alpha,1);

        mo_makeplot(2,2000,'myfile.txt',1);
        xf=-importdata('myfile.txt');
        x2=xf(:,2);
        N = length(x2);
        % plot(x(:,1),x2)


        [ah,bh]=fft_femm(xf);
        x1=linspace(0,2*pi,N/2); 
        figure
        hold on
        h = ((1:N/2) - 1)';
        Bgf = ah'*cos(h*x1) + bh'*sin(h*x1);
        plot(x1,Bgf)
        plot(x1,ah(np)*cos((np-1)*x1) + bh(np)*sin((np-1)*x1))

        figure
        bar(h,sqrt(ah.^2+bh.^2))
        axis([0 10 -0.5 1.5])

        Bg1(k) = sqrt(ah(2).^2+bh(2).^2);
    elseif np==1
        k = k+1;
        mo_addcontour(R_w+sign(R_w-R_m)*s_t/2,0);
        mo_addcontour((R_w+sign(R_w-R_m)*s_t/2)*cos((alpha)*2*pi/360),(R_w+sign(R_w-R_m)*s_t/2)*sin((alpha)*2*pi/360));
        mo_bendcontour(alpha,1);

        mo_makeplot(2,2000,'myfile.txt',1);
        xf=-importdata('myfile.txt');
%         xf = [xf ; -xf(2:size(xf,1),:)];
        xf = [xf ; -flip(xf(1:size(xf,1)-1,:))];
        x2=xf(:,2);
        N = length(x2);
        % plot(x(:,1),x2)


        [ah,bh]=fft_femm(xf);
        x1=linspace(0,2*pi,N/2); 
        figure
        hold on
        h = ((1:N/2) - 1)';
        Bgf = ah'*cos(h*x1) + bh'*sin(h*x1);
        plot(x1,Bgf)
        plot(x1,ah(np+1)*cos((np)*x1) + bh(np)*sin((np)*x1))

        figure
        bar(h,sqrt(ah.^2+bh.^2))
        axis([0 10 -0.5 1.5])

        Bg1(k) = sqrt(ah(2).^2+bh(2).^2);
        
    end
    
    %% heating map
    dt = delaunayTriangulation(real(z),imag(z));
%     limLow = min(B_dem);
%     limHi  = max(B_dem);
    limLow =0;
    limHi  = B_knee_iter(1,1);
    caxis([limLow, limHi])
    figure
    hold on
    h = trisurf(dt.ConnectivityList,real(z),imag(z),zeros(numel(z),1),B_dem);
    h.FaceColor = 'Interp';
    h.EdgeColor = 'none';
    colormap(flipud(gray));
    theta = linspace(0,2*pi/(2*p),1000); % circumferential discretization (leave the -pi/(2*p))
    plot(R_r*cos(linspace(theta(1),theta(end),100)),R_r*sin(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
    plot(R_m*cos(linspace(theta(1),theta(end),100)),R_m*sin(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
    plot([R_r*cos(theta(1)+(1-alpha_p)*pi/(2*p):pi/(p):theta(end));R_m*cos(theta(1)+(1-alpha_p)*pi/(2*p):pi/(p):theta(end))],[R_r*sin(theta(1)+(1-alpha_p)*pi/(2*p):pi/(p):theta(end));R_m*sin(theta(1)+(1-alpha_p)*pi/(2*p):pi/(p):theta(end))],'linewidth',0.8,'color','k');
    plot([R_r*cos(theta(1)+(1+alpha_p)*pi/(2*p):pi/(p):theta(end));R_m*cos(theta(1)+(1+alpha_p)*pi/(2*p):pi/(p):theta(end))],[R_r*sin(theta(1)+(1+alpha_p)*pi/(2*p):pi/(p):theta(end));R_m*sin(theta(1)+(1+alpha_p)*pi/(2*p):pi/(p):theta(end))],'linewidth',0.8,'color','k');
    plot([R_r*cos(theta(1)+(1-alpha_p1)*pi/(2*p):pi/(p):theta(end));R_m*cos(theta(1)+(1-alpha_p1)*pi/(2*p):pi/(p):theta(end))],[R_r*sin(theta(1)+(1-alpha_p1)*pi/(2*p):pi/(p):theta(end));R_m*sin(theta(1)+(1-alpha_p1)*pi/(2*p):pi/(p):theta(end))],'linewidth',0.8,'color','k');
    plot([R_r*cos(theta(1)+(1+alpha_p1)*pi/(2*p):pi/(p):theta(end));R_m*cos(theta(1)+(1+alpha_p1)*pi/(2*p):pi/(p):theta(end))],[R_r*sin(theta(1)+(1+alpha_p1)*pi/(2*p):pi/(p):theta(end));R_m*sin(theta(1)+(1+alpha_p1)*pi/(2*p):pi/(p):theta(end))],'linewidth',0.8,'color','k');
    plot(R_s*cos(linspace(theta(1),theta(end),100)),R_s*sin(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
    plot(R_w*cos(linspace(theta(1),theta(end),100)),R_w*sin(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
    plot(R_wi*cos(linspace(theta(1),theta(end),100)),R_wi*sin(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
    plot([R_wi*cos(theta(1)+pi/(6*p):pi/(3*p):theta(end)-pi/(6*p));R_w*cos(theta(1)+pi/(6*p):pi/(3*p):theta(end)-pi/(6*p))],[R_wi*sin(theta(1)+pi/(6*p):pi/(3*p):theta(end)-pi/(6*p));R_w*sin(theta(1)+pi/(6*p):pi/(3*p):theta(end)-pi/(6*p))],'linewidth',0.8,'color','k');
    caxis([limLow, limHi])
    axis image
    axis off
    title(['Demagnetization map (',num2str(Perc_demag),'% of the PM volume)'])
    camroll(90-90/p)
    
end



end

if length(Bg1) == 1
    error('Number of iterations wasn''t high eneough for convergence. Increase iter_num or Var (or both)')
end
Bg1_red_perc = (Bg1(1)-Bg1(2))/Bg1(1)*100;

%% magnets' flux density distribution
% integration line along the air gap which spans alpha (degrees)
% mo_addcontour(R_w+sign(R_w-R_m)*s_t/2,0);
% mo_addcontour((R_w+sign(R_w-R_m)*s_t/2)*cos((alpha)*2*pi/360),(R_w+sign(R_w-R_m)*s_t/2)*sin((alpha)*2*pi/360));
% mo_bendcontour(alpha,1);

% mo_addcontour(R_m,0);
% mo_addcontour((R_m)*cos((alpha)*2*pi/360),(R_m)*sin((alpha)*2*pi/360));
% mo_bendcontour(alpha,1);
% 
% mo_addcontour(R_m,0);
% mo_addcontour((R_m)*cos((alpha)*2*pi/360),(R_m)*sin((alpha)*2*pi/360));
% mo_bendcontour(alpha,1);

% mo_makeplot(2,2000,'myfile.txt',1);
% xf=-importdata('myfile.txt');
% x2=xf(:,2);
% N = length(x2);
% % plot(x(:,1),x2)
% 
% 
% [ah,bh]=fft_femm(xf);
% x1=linspace(0,2*pi,N/2); 
% figure
% hold on
% h = ((1:N/2) - 1)';
% Bgf = ah'*cos(h*x1) + bh'*sin(h*x1);
% plot(x1,Bgf)
% plot(x1,ah(np)*cos((np-1)*x1) + bh(np)*sin((np-1)*x1))
% % plot(x1,ah(2*np)*cos((2*np-1)*x1) + bh(2*np)*sin((2*np-1)*x1))
% % plot(x1,ah(3*np)*cos((3*np-1)*x1) + bh(3*np)*sin((3*np-1)*x1))
% figure
% bar(h,[sqrt(ah.^2+bh.^2)])
% axis([0 10 -0.5 1.5])
% 
% %% current density verification
% mo_groupselectblock(1001);
% Sslot = mo_blockintegral(5);
% mo_clearblock();
% Scu = Sslot*K_fill;
% nlitz = round(Scu/Slitz);
% % nc = round(nlitz/lpc); % conductors per slot
% Rs_femm =rho*(lstk+lew)*1e-3*Ns/Scu*1e-6;


closefemm