% clearvars
% clc
% close all




% back-emf and torque constant are obtained from the analytical PM field
% solution
Bemf_torque_constants
% clearvars -except k_t k_v p top runs
% inductance is obtained from the armature field solution
Inductance
% clearvars -except k_t k_v p L_ph top runs

% a target torque is set, the related target current is obtained from the
% torque constant. Given the current, and the operation in MTPA, the
% voltage equation is solved and the target voltage (reference to the
% inverter control) is obtained for  given operating peed.
k_t = k_t(1); % torque constant [Nm/Arms]
Torque = 3; % target torque [Nm]
I = Torque/k_t; % peak fundamental current [A]
L = L_ph; % phase inductance [H]
R = 10e-3; % phase resistance [ohm]
VDC = 45; % available DC voltage [V]
rpm = 2500; % rotational speed [rpm]
radsec = rpm*2*pi/60; % angular frequency [rad/sec]
radsec_e = p*radsec; % electrical angular frequency [rad/sec]

fsw = 10000;
% radsec_e = fsw*2*pi/60;

Efund = k_v(1)*radsec; % peak of fundamental back-emf [V]


% Voltage equation resulting form Id=0 & Iq=I
Ud = -radsec_e*L*I;
Uq = R*I + Efund;

U = sqrt(Ud^2 + Uq^2);

m = U/VDC;

% m = 1/sqrt(3);
% U = m*VDC;

if m>1/sqrt(3)
   error('Error. \n modulation index is higher than the limit 1/sqrt(3)')
end
    
% PARAMETERS


mult = 3;
T_fund = 2*pi/radsec_e;
Tsw = 1/fsw;
shift = -Tsw/2*(1); % IT MUST BE NEGATIVE (set it to -Tsw/2 to have the reference in the first Tsw at zero deg/rad)
Ns = ceil((mult*T_fund+abs(shift))*fsw)+1; % switching periods over the specified fundamental periods

bounds = linspace(0,(Ns)*Tsw,Ns+1)+shift; % switching boundary instants [s]
mid = bounds(1:end-1)+1/fsw/2; % mid point over switching periods [s]
mid_rad = mid*radsec_e; % mid point over switching periods [rad]
mid_rad = wrapTo2Pi(mid_rad);


Volt_avg_a = U*cos(mid_rad); % average inverter output
Volt_avg_b = U*cos(mid_rad-2/3*pi); % average inverter output

t1 = zeros(1,Ns);
t2 = zeros(1,Ns);
t0 = zeros(1,Ns);
Va = zeros(8,Ns);
Vb = zeros(8,Ns);
%% range 0<theta<60  FIRST SECTOR
t1(mid_rad<pi/3) = m*sqrt(3)*Tsw/2*sin(pi/3-mid_rad(mid_rad<pi/3)); % time to keep vector (1,0,0) over each switching period
t2(mid_rad<pi/3) = m*sqrt(3)*Tsw/2*sin(mid_rad(mid_rad<pi/3)); % time to keep vector (1,1,0) over each switching period
t0(mid_rad<pi/3) = Tsw/2-(t1(mid_rad<pi/3)+t2(mid_rad<pi/3)); % time to keep vector (0,0,0) over each switching period
Vsa = [0 ; 2/3 ; 1/3 ; 0]*VDC; 
Vsa = [Vsa ; flip(Vsa)]; % voltage applied over each switching instant
Vsa = repmat(Vsa,1,Ns);
Va(:,mid_rad<pi/3) = Vsa(:,mid_rad<pi/3);
Vsb = [0 ; -1/3 ; 1/3 ; 0]*VDC; 
Vsb = [Vsb ; flip(Vsb)]; % voltage applied over each switching instant
Vsb = repmat(Vsb,1,Ns);
Vb(:,mid_rad<pi/3) = Vsb(:,mid_rad<pi/3);


%% range 60<theta<120 SECOND SECTOR
t2(mid_rad>=pi/3&mid_rad<2*pi/3) = m*sqrt(3)*Tsw/2*sin(2*pi/3-mid_rad(mid_rad>=pi/3&mid_rad<2*pi/3)); % time to keep vector (1,1,0) over each switching period
t1(mid_rad>=pi/3&mid_rad<2*pi/3) = m*sqrt(3)*Tsw/2*sin(mid_rad(mid_rad>=pi/3&mid_rad<2*pi/3)-pi/3); % time to keep vector (0,1,0) over each switching period
t0(mid_rad>=pi/3&mid_rad<2*pi/3) = Tsw/2-(t1(mid_rad>=pi/3&mid_rad<2*pi/3)+t2(mid_rad>=pi/3&mid_rad<2*pi/3)); % time to keep vector (0,0,0) over each switching period
Vsa = [0 ; -1/3 ; 1/3 ; 0]*VDC; 
Vsa = [Vsa ; flip(Vsa)]; % voltage applied over each switching instant
Vsa = repmat(Vsa,1,Ns);
Va(:,mid_rad>=pi/3&mid_rad<2*pi/3) = Vsa(:,mid_rad>=pi/3&mid_rad<2*pi/3);
Vsb = [0 ; 2/3 ; 1/3 ; 0]*VDC; 
Vsb = [Vsb ; flip(Vsb)]; % voltage applied over each switching instant
Vsb = repmat(Vsb,1,Ns);
Vb(:,mid_rad>=pi/3&mid_rad<2*pi/3) = Vsb(:,mid_rad>=pi/3&mid_rad<2*pi/3);


%% range 120<theta<180 THIRD SECTOR
t1(mid_rad>=2*pi/3&mid_rad<pi) = m*sqrt(3)*Tsw/2*sin(mid_rad(mid_rad>=2*pi/3&mid_rad<pi)); % time to keep vector (0,1,0) over each switching period
t2(mid_rad>=2*pi/3&mid_rad<pi) = m*sqrt(3)*Tsw/2*sin(mid_rad(mid_rad>=2*pi/3&mid_rad<pi)+4*pi/3); % time to keep vector (0,1,1) over each switching period
t0(mid_rad>=2*pi/3&mid_rad<pi) = Tsw/2-(t1(mid_rad>=2*pi/3&mid_rad<pi)+t2(mid_rad>=2*pi/3&mid_rad<pi)); % time to keep vector (0,0,0) over each switching period
Vsa = [0 ; -1/3 ; -2/3 ; 0]*VDC; 
Vsa = [Vsa ; flip(Vsa)]; % voltage applied over each switching instant
Vsa = repmat(Vsa,1,Ns);
Va(:,mid_rad>=2*pi/3&mid_rad<pi) = Vsa(:,mid_rad>=2*pi/3&mid_rad<pi);
Vsb = [0 ; 2/3 ; 1/3 ; 0]*VDC; 
Vsb = [Vsb ; flip(Vsb)]; % voltage applied over each switching instant
Vsb = repmat(Vsb,1,Ns);
Vb(:,mid_rad>=2*pi/3&mid_rad<pi) = Vsb(:,mid_rad>=2*pi/3&mid_rad<pi);


%% range 180<theta<240 FOURTH SECTOR
t1(mid_rad>=pi&mid_rad<4*pi/3) = -m*sqrt(3)*Tsw/2*sin(mid_rad(mid_rad>=pi&mid_rad<4*pi/3)); % time to keep vector (0,0,1) over each switching period
t2(mid_rad>=pi&mid_rad<4*pi/3) = m*sqrt(3)*Tsw/2*sin(mid_rad(mid_rad>=pi&mid_rad<4*pi/3)-pi/3); % time to keep vector (0,1,1) over each switching period
t0(mid_rad>=pi&mid_rad<4*pi/3) = Tsw/2-(t1(mid_rad>=pi&mid_rad<4*pi/3)+t2(mid_rad>=pi&mid_rad<4*pi/3)); % time to keep vector (0,0,0) over each switching period
Vsa = [0 ; -1/3 ; -2/3 ; 0]*VDC; 
Vsa = [Vsa ; flip(Vsa)]; % voltage applied over each switching instant
Vsa = repmat(Vsa,1,Ns);
Va(:,mid_rad>=pi&mid_rad<4*pi/3) = Vsa(:,mid_rad>=pi&mid_rad<4*pi/3);
Vsb = [0 ; -1/3 ; 1/3 ; 0]*VDC; 
Vsb = [Vsb ; flip(Vsb)]; % voltage applied over each switching instant
Vsb = repmat(Vsb,1,Ns);
Vb(:,mid_rad>=pi&mid_rad<4*pi/3) = Vsb(:,mid_rad>=pi&mid_rad<4*pi/3);


%% range 240<theta<300 FIFTH SECTOR
t1(mid_rad>=4*pi/3&mid_rad<5*pi/3) = m*sqrt(3)*Tsw/2*sin(4*pi/3+mid_rad(mid_rad>=4*pi/3&mid_rad<5*pi/3)); % time to keep vector (0,0,1) over each switching period
t2(mid_rad>=4*pi/3&mid_rad<5*pi/3) = m*sqrt(3)*Tsw/2*sin(pi/3-mid_rad(mid_rad>=4*pi/3&mid_rad<5*pi/3)); % time to keep vector (1,0,1) over each switching period
t0(mid_rad>=4*pi/3&mid_rad<5*pi/3) = Tsw/2-(t1(mid_rad>=4*pi/3&mid_rad<5*pi/3)+t2(mid_rad>=4*pi/3&mid_rad<5*pi/3)); % time to keep vector (0,0,0) over each switching period
Vsa = [0 ; -1/3 ; 1/3 ; 0]*VDC; 
Vsa = [Vsa ; flip(Vsa)]; % voltage applied over each switching instant
Vsa = repmat(Vsa,1,Ns);
Va(:,mid_rad>=4*pi/3&mid_rad<5*pi/3) = Vsa(:,mid_rad>=4*pi/3&mid_rad<5*pi/3);
Vsb = [0 ; -1/3 ; -2/3 ; 0]*VDC; 
Vsb = [Vsb ; flip(Vsb)]; % voltage applied over each switching instant
Vsb = repmat(Vsb,1,Ns);
Vb(:,mid_rad>=4*pi/3&mid_rad<5*pi/3) = Vsb(:,mid_rad>=4*pi/3&mid_rad<5*pi/3);

%% range 300<theta<360 SIXTH SECTOR
t2(mid_rad>=5*pi/3&mid_rad<2*pi) = -m*sqrt(3)*Tsw/2*sin(mid_rad(mid_rad>=5*pi/3&mid_rad<2*pi)); % time to keep vector (1,0,1) over each switching period
t1(mid_rad>=5*pi/3&mid_rad<2*pi) = m*sqrt(3)*Tsw/2*sin(pi/3+mid_rad(mid_rad>=5*pi/3&mid_rad<2*pi)); % time to keep vector (1,0,0) over each switching period
t0(mid_rad>=5*pi/3&mid_rad<2*pi) = Tsw/2-(t1(mid_rad>=5*pi/3&mid_rad<2*pi)+t2(mid_rad>=5*pi/3&mid_rad<2*pi)); % time to keep vector (0,0,0) over each switching period
Vsa = [0 ; 2/3 ; 1/3 ; 0]*VDC; 
Vsa = [Vsa ; flip(Vsa)]; % voltage applied over each switching instant
Vsa = repmat(Vsa,1,Ns);
Va(:,mid_rad>=5*pi/3&mid_rad<2*pi) = Vsa(:,mid_rad>=5*pi/3&mid_rad<2*pi);
Vsb = [0 ; -1/3 ; -2/3 ; 0]*VDC; 
Vsb = [Vsb ; flip(Vsb)]; % voltage applied over each switching instant
Vsb = repmat(Vsb,1,Ns);
Vb(:,mid_rad>=5*pi/3&mid_rad<2*pi) = Vsb(:,mid_rad>=5*pi/3&mid_rad<2*pi);

A = [t0/2 ; t1 ; t2 ; t0/2]; % timings for the first half of each switching period
PERM = [A ; flip(A)]; % timings for each switching period (second half is the first half flipped -algorithm peculiarity-) FIRST SECTOR


times = reshape(PERM,[],1); % vector holding the time intervals where the different vectors are applied

Full_time = zeros(8*Ns+1,1)-Tsw/2;
Full_time(2:end) = Full_time(2:end) + cumsum(times);
mid_time = Full_time(2:end)-times/2; % mid-point of each switching interval

Volt_a = U*cos(radsec_e*mid_time);
Va = reshape(Va,[],1);
Vb = reshape(Vb,[],1);
% figure; stairs(Full_time(1:end-1),Va-Vb);
Va = Va-Volt_a;

Volt_b = U*cos(radsec_e*mid_time-2*pi/3);
Vb = Vb-Volt_b;

Full_time = Full_time(5:end);
Ripple_a = 1/(L)*cumsum(times.*Va);
% RIPPLE_a = [0;Ripple_a];
RIPPLE_a = Ripple_a(4:end);

Ripple_b = 1/(L)*cumsum(times.*Vb);
RIPPLE_b = Ripple_b(4:end);
x = [abs(diff(Full_time))>0;true];
Full_time = Full_time(x);
RIPPLE_a = RIPPLE_a(x);
RIPPLE_b = RIPPLE_b(x);

FULL_a = RIPPLE_a+I*cos(radsec_e*Full_time);
FULL_b = RIPPLE_b+I*cos(radsec_e*Full_time-2*pi/3);
FULL_c = -(FULL_a+FULL_b);

%% FFT
samp_t = 1/(50*fsw); % sampling time for the resampling
N = round(T_fund/samp_t); % number of samples
N = round2even(N); % an even number of samples is needed for the following code
samp = linspace(Full_time(1),T_fund,N); % new time-sampled vector
n = ((1:N/2-1))'; % harmonic orders
t = linspace(0,mult*T_fund,floor(mult*N/2)); % time discretization


%Phase a post-processing
series = timeseries(FULL_a,Full_time);
Current_per_a = resample(series,samp,'linear'); % one period current  resampled
[ah_a,bh_a] = fft_femm([Current_per_a.Time Current_per_a.Data]);
ch_a = sqrt(ah_a(2:end).^2+bh_a(2:end).^2);
phi_h_a = -atan2(bh_a(2:end),ah_a(2:end));

% Phase b post-processing
series = timeseries(FULL_b,Full_time);
Current_per_b = resample(series,samp,'linear'); % one period current resampled
[ah_b,bh_b] = fft_femm([Current_per_b.Time Current_per_b.Data]);
ch_b = sqrt(ah_b(2:end).^2+bh_b(2:end).^2);
phi_h_b = -atan2(bh_b(2:end),ah_b(2:end));

% Phase c post-processing
series = timeseries(FULL_c,Full_time);
Current_per_c = resample(series,samp,'linear'); % one period current resampled
[ah_c,bh_c] = fft_femm([Current_per_c.Time Current_per_c.Data]);
ch_c = sqrt(ah_c(2:end).^2+bh_c(2:end).^2);
phi_h_c = -atan2(bh_c(2:end),ah_c(2:end));


% Symmetric components analysis
coeff = cos(2/3*pi)+1i*sin(2/3*pi);
coeff_sq = cos(4/3*pi)+1i*sin(4/3*pi);

dir = 1/3*((ah_a(2:end)-1i*bh_a(2:end))+(ah_b(2:end)-1i*bh_b(2:end)).*coeff+(ah_c(2:end)-1i*bh_c(2:end)).*coeff_sq);
inv = 1/3*((ah_a(2:end)-1i*bh_a(2:end))+(ah_b(2:end)-1i*bh_b(2:end)).*coeff_sq+(ah_c(2:end)-1i*bh_c(2:end)).*coeff);
zero = 1/3*((ah_a(2:end)-1i*bh_a(2:end))+(ah_b(2:end)-1i*bh_b(2:end))+(ah_c(2:end)-1i*bh_c(2:end)));
ch_a_pos = dir;
ch_b_pos = dir.*coeff_sq;
ch_c_pos = dir.*coeff;

ch_a_neg = inv;
ch_b_neg = inv.*coeff;
ch_c_neg = inv.*coeff_sq;

pos = abs(ch_a_pos);
neg = abs(ch_a_neg);
phi_p = angle(ch_a_pos);
phi_n = angle(ch_a_neg);
