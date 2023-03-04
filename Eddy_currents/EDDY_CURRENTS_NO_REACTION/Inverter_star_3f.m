% clearvars
% clc
% close all



% top = 'Inrunner'; % Choose either 'Inrunner' or 'Outrunner'

% back-emf and torque constant are obtained from the analytical PM field
% solution
Bemf_torque_constants
clearvars -except k_t k_v p top runs
% inductance is obtained from the armature field solution
Inductance
clearvars -except k_t k_v p L_ph top runs

Inrunner
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


T_fund = 2*pi/radsec_e;
Tsw = 1/fsw;
Ns = ceil(1*T_fund*fsw)+1; % switching periods over the specified fundamental periods

bounds = linspace(0,Ns*Tsw,Ns+1)-Tsw/2; % switching boundary instants [s]
mid = bounds(1:end-1)+1/fsw/2; % mid point over switching periods [s]
mid_rad = mid*radsec_e; % mid point over switching periods [rad]
mid_rad(mid_rad>=2*pi) = mid_rad(mid_rad>=2*pi)-2*pi;


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
% A = [t0/2 ; t1 ; t2 ; t0/2]; % timings for the first half of each switching period
% PERM1 = [A ; flip(A)]; % timings for each switching period (second half is the first half flipped -algorithm peculiarity-) FIRST SECTOR
Vsa = [0 ; 2/3 ; 1/3 ; 0]*VDC; 
Vsa = [Vsa ; flip(Vsa)]; % voltage applied over each switching instant
Vsa = repmat(Vsa,1,Ns);
Va(:,mid_rad<pi/3) = Vsa(:,mid_rad<pi/3);
Vsb = [0 ; -1/3 ; 1/3 ; 0]*VDC; 
Vsb = [Vsb ; flip(Vsb)]; % voltage applied over each switching instant
Vsb = repmat(Vsb,1,Ns);
Vb(:,mid_rad<pi/3) = Vsb(:,mid_rad<pi/3);

% Vl1 = [0 ; 1 ; 1 ; 1]*VDC; 
% Vl1 = [Vl1 ; flip(Vl1)]; % voltage applied over each switching instant
% Vl1 = repmat(Vl1,1,length(mid_rad(mid_rad<pi/3)));

%% range 60<theta<120 SECOND SECTOR
t2(mid_rad>=pi/3&mid_rad<2*pi/3) = m*sqrt(3)*Tsw/2*sin(2*pi/3-mid_rad(mid_rad>=pi/3&mid_rad<2*pi/3)); % time to keep vector (1,1,0) over each switching period
t1(mid_rad>=pi/3&mid_rad<2*pi/3) = m*sqrt(3)*Tsw/2*sin(mid_rad(mid_rad>=pi/3&mid_rad<2*pi/3)-pi/3); % time to keep vector (0,1,0) over each switching period
t0(mid_rad>=pi/3&mid_rad<2*pi/3) = Tsw/2-(t1(mid_rad>=pi/3&mid_rad<2*pi/3)+t2(mid_rad>=pi/3&mid_rad<2*pi/3)); % time to keep vector (0,0,0) over each switching period
% A = [t0/2 ; t2 ; t1 ; t0/2]; % timings for the first half of each switching period
% PERM2 = [A ; flip(A)]; % timings for each switching period (second half is the first half flipped -algorithm peculiarity-) 
Vsa = [0 ; -1/3 ; 1/3 ; 0]*VDC; 
Vsa = [Vsa ; flip(Vsa)]; % voltage applied over each switching instant
Vsa = repmat(Vsa,1,Ns);
Va(:,mid_rad>=pi/3&mid_rad<2*pi/3) = Vsa(:,mid_rad>=pi/3&mid_rad<2*pi/3);
Vsb = [0 ; 2/3 ; 1/3 ; 0]*VDC; 
Vsb = [Vsb ; flip(Vsb)]; % voltage applied over each switching instant
Vsb = repmat(Vsb,1,Ns);
Vb(:,mid_rad>=pi/3&mid_rad<2*pi/3) = Vsb(:,mid_rad>=pi/3&mid_rad<2*pi/3);

% Vl2 = [0 ; 0 ; 1 ; 1]*VDC; 
% Vl2 = [Vl2 ; flip(Vl2)]; % voltage applied over each switching instant
% Vl2 = repmat(Vl2,1,length(mid_rad(mid_rad>=pi/3&mid_rad<2*pi/3)));

%% range 120<theta<180 THIRD SECTOR
t1(mid_rad>=2*pi/3&mid_rad<pi) = m*sqrt(3)*Tsw/2*sin(mid_rad(mid_rad>=2*pi/3&mid_rad<pi)); % time to keep vector (0,1,0) over each switching period
t2(mid_rad>=2*pi/3&mid_rad<pi) = m*sqrt(3)*Tsw/2*sin(mid_rad(mid_rad>=2*pi/3&mid_rad<pi)+4*pi/3); % time to keep vector (0,1,1) over each switching period
t0(mid_rad>=2*pi/3&mid_rad<pi) = Tsw/2-(t1(mid_rad>=2*pi/3&mid_rad<pi)+t2(mid_rad>=2*pi/3&mid_rad<pi)); % time to keep vector (0,0,0) over each switching period
% A = [t0/2 ; t1 ; t2 ; t0/2]; % timings for the first half of each switching period
% PERM3 = [A ; flip(A)]; % timings for each switching period (second half is the first half flipped -algorithm peculiarity-) SECOND SECTOR
Vsa = [0 ; -1/3 ; -2/3 ; 0]*VDC; 
Vsa = [Vsa ; flip(Vsa)]; % voltage applied over each switching instant
Vsa = repmat(Vsa,1,Ns);
Va(:,mid_rad>=2*pi/3&mid_rad<pi) = Vsa(:,mid_rad>=2*pi/3&mid_rad<pi);
Vsb = [0 ; 2/3 ; 1/3 ; 0]*VDC; 
Vsb = [Vsb ; flip(Vsb)]; % voltage applied over each switching instant
Vsb = repmat(Vsb,1,Ns);
Vb(:,mid_rad>=2*pi/3&mid_rad<pi) = Vsb(:,mid_rad>=2*pi/3&mid_rad<pi);

% Vl3 = [0 ; 0 ; 0 ; 1]*VDC; 
% Vl3 = [Vl3 ; flip(Vl3)]; % voltage applied over each switching instant
% Vl3 = repmat(Vl3,1,length(mid_rad(mid_rad>=2*pi/3&mid_rad<pi)));

%% range 180<theta<240 FOURTH SECTOR
t1(mid_rad>=pi&mid_rad<4*pi/3) = -m*sqrt(3)*Tsw/2*sin(mid_rad(mid_rad>=pi&mid_rad<4*pi/3)); % time to keep vector (0,0,1) over each switching period
t2(mid_rad>=pi&mid_rad<4*pi/3) = m*sqrt(3)*Tsw/2*sin(mid_rad(mid_rad>=pi&mid_rad<4*pi/3)-pi/3); % time to keep vector (0,1,1) over each switching period
t0(mid_rad>=pi&mid_rad<4*pi/3) = Tsw/2-(t1(mid_rad>=pi&mid_rad<4*pi/3)+t2(mid_rad>=pi&mid_rad<4*pi/3)); % time to keep vector (0,0,0) over each switching period
% A = [t0/2 ; t1 ; t2 ; t0/2]; % timings for the first half of each switching period
% PERM4 = [A ; flip(A)]; % timings for each switching period (second half is the first half flipped -algorithm peculiarity-) SECOND SECTOR
Vsa = [0 ; -1/3 ; -2/3 ; 0]*VDC; 
Vsa = [Vsa ; flip(Vsa)]; % voltage applied over each switching instant
Vsa = repmat(Vsa,1,Ns);
Va(:,mid_rad>=pi&mid_rad<4*pi/3) = Vsa(:,mid_rad>=pi&mid_rad<4*pi/3);
Vsb = [0 ; -1/3 ; 1/3 ; 0]*VDC; 
Vsb = [Vsb ; flip(Vsb)]; % voltage applied over each switching instant
Vsb = repmat(Vsb,1,Ns);
Vb(:,mid_rad>=pi&mid_rad<4*pi/3) = Vsb(:,mid_rad>=pi&mid_rad<4*pi/3);

% Vl4 = [0 ; 0 ; 0 ; 1]*VDC; 
% Vl4 = [Vl4 ; flip(Vl4)]; % voltage applied over each switching instant
% Vl4 = repmat(Vl4,1,length(mid_rad(mid_rad>=pi&mid_rad<4*pi/3)));

%% range 240<theta<300 FIFTH SECTOR
t1(mid_rad>=4*pi/3&mid_rad<5*pi/3) = m*sqrt(3)*Tsw/2*sin(4*pi/3+mid_rad(mid_rad>=4*pi/3&mid_rad<5*pi/3)); % time to keep vector (0,0,1) over each switching period
t2(mid_rad>=4*pi/3&mid_rad<5*pi/3) = m*sqrt(3)*Tsw/2*sin(pi/3-mid_rad(mid_rad>=4*pi/3&mid_rad<5*pi/3)); % time to keep vector (1,0,1) over each switching period
t0(mid_rad>=4*pi/3&mid_rad<5*pi/3) = Tsw/2-(t1(mid_rad>=4*pi/3&mid_rad<5*pi/3)+t2(mid_rad>=4*pi/3&mid_rad<5*pi/3)); % time to keep vector (0,0,0) over each switching period
% A = [t0/2 ; t1 ; t2 ; t0/2]; % timings for the first half of each switching period
% PERM5 = [A ; flip(A)]; % timings for each switching period (second half is the first half flipped -algorithm peculiarity-) SECOND SECTOR
Vsa = [0 ; -1/3 ; 1/3 ; 0]*VDC; 
Vsa = [Vsa ; flip(Vsa)]; % voltage applied over each switching instant
Vsa = repmat(Vsa,1,Ns);
Va(:,mid_rad>=4*pi/3&mid_rad<5*pi/3) = Vsa(:,mid_rad>=4*pi/3&mid_rad<5*pi/3);
Vsb = [0 ; -1/3 ; -2/3 ; 0]*VDC; 
Vsb = [Vsb ; flip(Vsb)]; % voltage applied over each switching instant
Vsb = repmat(Vsb,1,Ns);
Vb(:,mid_rad>=4*pi/3&mid_rad<5*pi/3) = Vsb(:,mid_rad>=4*pi/3&mid_rad<5*pi/3);

% Vl5 = [0 ; 0 ; 1 ; 1]*VDC; 
% Vl5 = [Vl5 ; flip(Vl5)]; % voltage applied over each switching instant
% Vl5 = repmat(Vl5,1,length(mid_rad(mid_rad>=4*pi/3&mid_rad<5*pi/3)));

%% range 300<theta<360 SIXTH SECTOR
t2(mid_rad>=5*pi/3&mid_rad<2*pi) = -m*sqrt(3)*Tsw/2*sin(mid_rad(mid_rad>=5*pi/3&mid_rad<2*pi)); % time to keep vector (1,0,1) over each switching period
t1(mid_rad>=5*pi/3&mid_rad<2*pi) = m*sqrt(3)*Tsw/2*sin(pi/3+mid_rad(mid_rad>=5*pi/3&mid_rad<2*pi)); % time to keep vector (1,0,0) over each switching period
t0(mid_rad>=5*pi/3&mid_rad<2*pi) = Tsw/2-(t1(mid_rad>=5*pi/3&mid_rad<2*pi)+t2(mid_rad>=5*pi/3&mid_rad<2*pi)); % time to keep vector (0,0,0) over each switching period
% A = [t0/2 ; t2 ; t1 ; t0/2]; % timings for the first half of each switching period
% PERM6 = [A ; flip(A)]; % timings for each switching period (second half is the first half flipped -algorithm peculiarity-) SECOND SECTOR
Vsa = [0 ; 2/3 ; 1/3 ; 0]*VDC; 
Vsa = [Vsa ; flip(Vsa)]; % voltage applied over each switching instant
Vsa = repmat(Vsa,1,Ns);
Va(:,mid_rad>=5*pi/3&mid_rad<2*pi) = Vsa(:,mid_rad>=5*pi/3&mid_rad<2*pi);
Vsb = [0 ; -1/3 ; -2/3 ; 0]*VDC; 
Vsb = [Vsb ; flip(Vsb)]; % voltage applied over each switching instant
Vsb = repmat(Vsb,1,Ns);
Vb(:,mid_rad>=5*pi/3&mid_rad<2*pi) = Vsb(:,mid_rad>=5*pi/3&mid_rad<2*pi);

% Vl6 = [0 ; 1 ; 1 ; 1]*VDC; 
% Vl6 = [Vl6 ; flip(Vl6)]; % voltage applied over each switching instant
% Vl6 = repmat(Vl6,1,length(mid_rad(mid_rad>=5*pi/3&mid_rad<2*pi)));

% %% range 0<theta<60  FIRST SECTOR
% t1 = m*sqrt(3)*Tsw/2*sin(pi/3-mid_rad(mid_rad>=2*pi)-2*pi); % time to keep vector (1,0,0) over each switching period
% t2 = m*sqrt(3)*Tsw/2*sin(mid_rad(mid_rad>=2*pi)-2*pi); % time to keep vector (1,1,0) over each switching period
% t0 = Tsw/2-(t1+t2); % time to keep vector (0,0,0) over each switching period
% A = [t0/2 ; t1 ; t2 ; t0/2]; % timings for the first half of each switching period
% PERM7 = [A ; flip(A)]; % timings for each switching period (second half is the first half flipped -algorithm peculiarity-) FIRST SECTOR
% V7a = [0 ; 2/3 ; 1/3 ; 0]*VDC; 
% V7a = [V7a ; flip(V7a)]; % voltage applied over each switching instant
% V7a = repmat(V7a,1,length(mid_rad(mid_rad>=2*pi)));
% V7b = [0 ; -1/3 ; 1/3 ; 0]*VDC; 
% V7b = [V7b ; flip(V7b)]; % voltage applied over each switching instant
% V7b = repmat(V7b,1,length(mid_rad(mid_rad>=2*pi)));
% Vl7 = [0 ; 1 ; 1 ; 1]*VDC; 
% Vl7 = [Vl7 ; flip(Vl7)]; % voltage applied over each switching instant
% Vl7 = repmat(Vl7,1,length(mid_rad(mid_rad>=2*pi)));


A = [t0/2 ; t1 ; t2 ; t0/2]; % timings for the first half of each switching period
PERM = [A ; flip(A)]; % timings for each switching period (second half is the first half flipped -algorithm peculiarity-) FIRST SECTOR

% PERM = [PERM1 PERM2 PERM3 PERM4 PERM5 PERM6 PERM7];
times = reshape(PERM,[],1); % vector holding the time intervals where the different vectors are applied
% Va = [Vsa Vsa Vsa Vsa Vsa Vsa V7a];
% Vb = [Vsb Vsb Vsb Vsb Vsb Vsb V7b];
% Vla = [Vl1 Vl2 Vl3 Vl4 Vl5 Vl6 Vl7];

Full_time = zeros(8*Ns+1,1)-Tsw/2;
Full_time(2:end) = Full_time(2:end) + cumsum(times);
mid_time = Full_time(2:end)-times/2; % mid-point of each switching interval

Volt_a = U*cos(radsec_e*mid_time);
Va = reshape(Va,[],1);
% figure; stairs(Va);
Va = Va-Volt_a;

Volt_b = U*cos(radsec_e*mid_time-2*pi/3);
Vb = reshape(Vb,[],1);
Vb = Vb-Volt_b;

Full_time = Full_time(5:end);
Ripple_a = 1/(L)*cumsum(times.*Va);
% RIPPLE_a = [0;Ripple_a];
RIPPLE_a = Ripple_a(4:end);
% idx = find(Full_time>=T_fund);
% FUND_curr = interp1(Full_time(idx-1:idx),RIPPLE_a(idx-1:idx),T_fund); 
% RIPPLE_a = [RIPPLE_a(1:idx-1); FUND_curr];

Ripple_b = 1/(L)*cumsum(times.*Vb);
% RIPPLE_b = [0;Ripple_b];
RIPPLE_b = Ripple_b(4:end);
% FUND_curr = interp1(Full_time(idx-1:idx),RIPPLE_b(idx-1:idx),T_fund);
% RIPPLE_b = [RIPPLE_b(1:idx-1); FUND_curr];

% Full_time = [Full_time(1:idx-1); T_fund];
x = [abs(diff(Full_time))>0;true];
Full_time = Full_time(x);
RIPPLE_a = RIPPLE_a(x);
RIPPLE_b = RIPPLE_b(x);

FULL_a = RIPPLE_a+I*cos(radsec_e*Full_time);
FULL_b = RIPPLE_b+I*cos(radsec_e*Full_time-2*pi/3);
FULL_c = -(FULL_a+FULL_b);

% writematrix([Full_time, FULL_a] , 'I1.txt')
% writematrix([Full_time, FULL_b], 'I2.txt')
% writematrix([Full_time, FULL_c], 'I3.txt')
% times = [linspace(0,length(Full_time)-1,length(Full_time))' Full_time];
% writematrix(times , 'Full_time_vec')
% TEST_rms = sqrt(sum(FULL_a.^2)/length(FULL_a));

% figure
% plot(Full_time,RIPPLE_a)
% figure
% plot(Full_time,FULL_a)
% figure
% plot(Full_time,FULL_a,Full_time,FULL_b,Full_time,FULL_c)
% hold on
% plot(Full_time,FULL_b)
% plot(Full_time,FULL_c)


%% FFT

samp_t = 1/(100*fsw); % sampling time for the resampling

% Phase current resampled for FFT
N = round(T_fund/samp_t); % number of samples
N = round2even(N); % an even number of samples is needed for the following code
samp = linspace(Full_time(1),T_fund,N); % new time-sampled vector

series = timeseries(FULL_a,Full_time);
int = interp1(Full_time,FULL_a,samp);
Current_per_a = resample(series,samp,'linear'); % one period current 
                                              %   resampled

[ah_a,bh_a] = fft_femm([Current_per_a.Time Current_per_a.Data]);


n = ((1:N/2-1))'; % harmonic orders

ch_a = sqrt(ah_a(2:end).^2+bh_a(2:end).^2);
phi_h_a = -atan2(bh_a(2:end),ah_a(2:end));


% ch_a = [0 0 ch_a(2)/2]' ;
% 
% I_fft_a_ph = ch_a'*cos(h(2:4).*radsec_e*t+phi_h_a(2:4));
% I_fft_b_ph = ch_a'*cos(h(2:4).*radsec_e*t+phi_h_a(2:4)-h(2:4)*2*pi/3);
% I_fft_c_ph = ch_a'*cos(h(2:4).*radsec_e*t+phi_h_a(2:4)+h(2:4)*2*pi/3);




series = timeseries(FULL_b,Full_time);
int = interp1(Full_time,FULL_b,samp);
Current_per_b = resample(series,samp,'linear'); % one period current 
                                              %   resampled
% 
[ah_b,bh_b] = fft_femm([Current_per_b.Time Current_per_b.Data]);
% % [ah,bh] = fft_femm([Full_time FULL_a]);
% 
% I_fft_b = ah_b'*cos(h.*radsec_e*t)+bh_b'*sin(h.*radsec_e*t);
% % toc
% 
% figure
% plot(t,I_fft_b)
% 
ch_b = sqrt(ah_b(2:end).^2+bh_b(2:end).^2);
phi_h_b = -atan2(bh_b(2:end),ah_b(2:end));
% I_fft_b_ph = ch_b'*cos(h.*radsec_e*t+phi_h_b);
% 
series = timeseries(FULL_c,Full_time);
int = interp1(Full_time,FULL_c,samp);
Current_per_c = resample(series,samp,'linear'); % one period current 
                                              %   resampled
% 
[ah_c,bh_c] = fft_femm([Current_per_c.Time Current_per_c.Data]);
% % [ah,bh] = fft_femm([Full_time FULL_a]);
% 
% I_fft_c = ah_c'*cos(h.*radsec_e*t)+bh_c'*sin(h.*radsec_e*t);
% % toc
% 
% figure
% plot(t,I_fft_c)
% 
ch_c = sqrt(ah_c(2:end).^2+bh_c(2:end).^2);
phi_h_c = -atan2(bh_c(2:end),ah_c(2:end));

h = ((1:N/2-1))'; % harmonic orders
t = linspace(0,T_fund,floor(N/2));

% I_fft_a_ph = ch_a'*cos(h.*radsec_e*t+phi_h_a);
% I_fft_b_ph = ch_b'*cos(h.*radsec_e*t+phi_h_b);
% I_fft_c_ph = ch_c'*cos(h.*radsec_e*t+phi_h_c);
% SAVE THE DATA FOR COMSOL
% writematrix([t', I_fft_a_ph'] , 'I1_fft.txt')
% writematrix([t', I_fft_b_ph'], 'I2_fft.txt')
% writematrix([t', I_fft_c_ph'], 'I3_fft.txt')
% times = [linspace(0,length(t)-1,length(t))' t'];
% writematrix(times , 'Full_time_fft')

% SEQUNCE SELECTION/DELETION

% coeff = cos(-2/3*pi*h)+1i*sin(-2/3*pi*h);
% coeff_sq = cos(-4/3*pi*h)+1i*sin(-4/3*pi*h);

coeff = cos(2/3*pi)+1i*sin(2/3*pi);
coeff_sq = cos(4/3*pi)+1i*sin(4/3*pi);

dir = 1/3*((ah_a(2:end)-1i*bh_a(2:end))+(ah_b(2:end)-1i*bh_b(2:end)).*coeff+(ah_c(2:end)-1i*bh_c(2:end)).*coeff_sq);
inv = 1/3*((ah_a(2:end)-1i*bh_a(2:end))+(ah_b(2:end)-1i*bh_b(2:end)).*coeff_sq+(ah_c(2:end)-1i*bh_c(2:end)).*coeff);
zero = 1/3*((ah_a(2:end)-1i*bh_a(2:end))+(ah_b(2:end)-1i*bh_b(2:end))+(ah_c(2:end)-1i*bh_c(2:end)));
% mask = abs(dir)>abs(inv); % ones where the positive sequence has major contribution
% mask = 1; % ones where the positive sequence has major contribution
% filt_d = 
% filt = (abs(dir)>1e-11 & abs(inv)>1e-11);
% filt = 1;

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

% figure;hold on
% plot(t,(pos)'*cos(h.*radsec_e*t+phi_p)+(neg)'*cos(h.*radsec_e*t+phi_n))
% plot(t,(pos)'*cos(h.*radsec_e*t+phi_p-2/3*pi)+(neg)'*cos(h.*radsec_e*t+phi_n+2/3*pi))
% plot(t,(pos)'*cos(h.*radsec_e*t+phi_p-4/3*pi)+(neg)'*cos(h.*radsec_e*t+phi_n+4/3*pi))

% I_fft_c_ph = ch_c'*cos(h.*radsec_e*t+phi_h_c);
% 
% figure
% plot(t,ah_c(4)'*cos(h(4).*radsec_e*t)+bh_c(4)'*sin(h(4).*radsec_e*t)+ah_b(4)'*cos(h(4).*radsec_e*t)+bh_b(4)'*sin(h(4).*radsec_e*t)+ah_a(4)'*cos(h(4).*radsec_e*t)+bh_a(4)'*sin(h(4).*radsec_e*t))
