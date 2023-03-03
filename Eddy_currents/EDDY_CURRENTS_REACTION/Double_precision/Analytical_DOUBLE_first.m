%% Machine parameters
% clearvars
% clc
% close all



%% inrunner example
% mp.Digits(500);
top = 'Inrunner'; % Choose either 'Inrunner' or 'Outrunner'

% inductance is obtained from the armature field solution
Inverter_star_3f
clearvars -except n top pos neg phi_p phi_n N t runs

plotting = "no";

run(top)
% pi = mp('pi');
mu_0 = 4*pi*1e-7; % air permeability
cond = 0.667*1e6;
speed_rpm = 2500;
speed_radsec = speed_rpm*(2*pi)/60;
%% useful indices for series harmonics definition

m_J = 3; % total number of harmonics better to be limited to 7 at most. 
          % m_J=3 gives good accuracy, fast computation and low memory usage 
m = (1:2:2*m_J);

pos = repmat(pos,1,length(m));
neg = repmat(neg,1,length(m));
phi_p = repmat(phi_p,1,length(m));
phi_n = repmat(phi_n,1,length(m));
n = repmat(n,1,length(m));
m = repmat(m,size(n,1),1);

TL_m = m-1;
TL_p = m+1;

k_minus_n = m-n;
k_plus_n = n+m;
pos_m = k_minus_n;
neg_m = pos_m;
pos_p = k_plus_n;
neg_p = pos_p;
% pos_m = pos_m.*(mod(TL_m,3)==0);
% neg_m(mod(TL_p,3)~=0)=0;
% neg_p = neg_p.*(mod(TL_m,3)==0);
% pos_p(mod(TL_p,3)~=0)=0;

filt = [(mod(TL_m,3)==0), mod(TL_p,3)==0, mod(TL_p,3)==0, (mod(TL_m,3)==0)];

H = [pos_m, pos_p, neg_m, neg_p].*filt;
ch = [pos, pos, neg, neg].*filt;
phi = [-phi_p, phi_p, -phi_n, phi_n].*filt;

m = [m, m, m, m]*p;
n = [n, n, n, n];

% consider only those harmonics generating eddy currents loss
m(filt==0) = [];
m = m';

n(filt==0) = [];
n = n';

H(filt==0) = [];
H = H';

ch(filt==0) = [];
ch = ch';

phi(filt==0) = [];
phi = phi';

m_dim = length(m);


tau = sqrt(-1i*cond*mu_r*mu_0*(H)*speed_radsec*p);
% tau = (1-1i)./sqrt(2./(cond*mu_r*mu_0*H*speed_radsec*p));

%% Harmonic fliters for Gibbs phenomenon reduction
order_n = 0;
order_m = 0;
sigma_n = (sin(pi*n./n(end))./(pi*n./n(end))).^order_n; % Lanczos sigma for Gibbs phenomenon reduction
sigma_m = (sin(pi*m./m(end))./(pi*m./m(end))).^order_m; % Lanczos sigma for Gibbs phenomenon reduction

ch = sigma_n.*ch;

%% Current density distribution
%% INDEPENDENT HARMONICS CONTRIBUTION

J_m = 4*3*p*N_tc*q*ch/b/S_Ph/(2*pi);                                        % constant coefficient for current density distribution
J_m = (J_m.*sin(pi*m/(6*p))./m);
A_mJ = mu_0*sigma_m.*J_m./(m.^2-4);
if p == 2
   A_mJ(:,1) = -mu_0*J_m(:,1)/4; 
end

%% Analytical coefficient expressions in the magnets region
% J_r = besselj(m,R_r*tau);
% Y_r = bessely(m,R_r*tau);
% J1_r = -tau.*(besselj(m+1,R_r*tau)-besselj(m-1,R_r*tau))./2;
% Y1_r = -tau.*(bessely(m+1,R_r*tau)-bessely(m-1,R_r*tau))./2;
% J_r(isnan(J_r)) = 0; 
% J1_r(isnan(J1_r)) = 0; 
% Y_r(isnan(Y_r)) = 0; 
% Y1_r(isnan(Y1_r)) = 0;
% 
% J_m = besselj(m,R_m*tau);
% Y_m = bessely(m,R_m*tau);
% J1_m = -tau.*(besselj(m+1,R_m*tau)-besselj(m-1,R_m*tau))./2;
% Y1_m = -tau.*(bessely(m+1,R_m*tau)-bessely(m-1,R_m*tau))./2;
% J_m(isnan(J_m)) = 0; 
% J1_m(isnan(J1_m)) = 0; 
% Y_m(isnan(Y_m)) = 0; 
% Y1_m(isnan(Y1_m)) = 0;
% 
% 
% COMM = (2.*R_wi.^2.*(R_w/R_wi).^m - R_w.^2.*m - 2.*R_w.^2 - 2.*R_w.^2.*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) + R_wi.^2.*m.*(R_w/R_wi).^m + 2.*R_wi.^2.*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^m - R_wi.^2.*m.*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^m + R_w.^2.*m.*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m));
% COMM = A_mJ.*R_m.*(R_m/R_w).^(m - 1).*COMM;
% COMM = COMM./(R_w.*(J_m.*Y_r.*m.^2 - J_r.*Y_m.*m.^2 - J_m.*Y_r.*m.^2.*(R_i/R_r).^(2.*m) + J_r.*Y_m.*m.^2.*(R_i/R_r).^(2.*m) - J1_m.*R_m.*R_r.*Y1_r + J1_r.*R_m.*R_r.*Y1_m + J1_m.*R_m.*Y_r.*m - J_r.*R_m.*Y1_m.*m + J1_r.*R_r.*Y_m.*m - J_m.*R_r.*Y1_r.*m - J1_m.*R_m.*R_r.*Y1_r.*(R_i/R_r).^(2.*m) + J1_r.*R_m.*R_r.*Y1_m.*(R_i/R_r).^(2.*m) - J1_m.*R_m.*Y_r.*m.*(R_i/R_r).^(2.*m) + J_r.*R_m.*Y1_m.*m.*(R_i/R_r).^(2.*m) + J1_r.*R_r.*Y_m.*m.*(R_i/R_r).^(2.*m) - J_m.*R_r.*Y1_r.*m.*(R_i/R_r).^(2.*m) - J_m.*Y_r.*m.^2.*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) + J_r.*Y_m.*m.^2.*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) + J_m.*Y_r.*m.^2.*(R_i/R_r).^(2.*m).*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) - J_r.*Y_m.*m.^2.*(R_i/R_r).^(2.*m).*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) - J1_m.*R_m.*R_r.*Y1_r.*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) + J1_r.*R_m.*R_r.*Y1_m.*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) + J1_m.*R_m.*Y_r.*m.*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) - J_r.*R_m.*Y1_m.*m.*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) - J1_r.*R_r.*Y_m.*m.*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) + J_m.*R_r.*Y1_r.*m.*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) - J1_m.*R_m.*Y_r.*m.*(R_i/R_r).^(2.*m).*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) + J_r.*R_m.*Y1_m.*m.*(R_i/R_r).^(2.*m).*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) - J1_r.*R_r.*Y_m.*m.*(R_i/R_r).^(2.*m).*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) + J_m.*R_r.*Y1_r.*m.*(R_i/R_r).^(2.*m).*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) - J1_m.*R_m.*R_r.*Y1_r.*(R_i/R_r).^(2.*m).*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) + J1_r.*R_m.*R_r.*Y1_m.*(R_i/R_r).^(2.*m).*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m)));

J_r = besselj(m,(R_r*tau));
Y_r = bessely(m,(R_r*tau));
J1_r = -tau.*(besselj(m+1,(R_r*tau))-besselj(m-1,(R_r*tau)))./2;
Y1_r = -tau.*(bessely(m+1,(R_r*tau))-bessely(m-1,(R_r*tau)))./2;
J_r(isnan(J_r)) = 0; 
J1_r(isnan(J1_r)) = 0; 
Y_r(isnan(Y_r)) = 0; 
Y1_r(isnan(Y1_r)) = 0;

J_m = besselj(m,(R_m*tau));
Y_m = bessely(m,(R_m*tau));
J1_m = -tau.*(besselj(m+1,(R_m*tau))-besselj(m-1,(R_m*tau)))./2;
Y1_m = -tau.*(bessely(m+1,(R_m*tau))-bessely(m-1,(R_m*tau)))./2;
J_m(isnan(J_m)) = 0; 
J1_m(isnan(J1_m)) = 0; 
Y_m(isnan(Y_m)) = 0; 
Y1_m(isnan(Y1_m)) = 0;

COMM_double = (2.*R_wi.^2.*(R_w/R_wi).^m - R_w.^2.*m - 2.*R_w.^2 - 2.*R_w.^2.*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) + R_wi.^2.*m.*(R_w/R_wi).^m + 2.*R_wi.^2.*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^m - R_wi.^2.*m.*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^m + R_w.^2.*m.*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m));
COMM_double = A_mJ.*R_m.*(R_m/R_w).^(m - 1).*COMM_double;
COMM_double = COMM_double./(R_w.*(J_m.*Y_r.*m.^2 - J_r.*Y_m.*m.^2 - J_m.*Y_r.*m.^2.*(R_i/R_r).^(2.*m) + J_r.*Y_m.*m.^2.*(R_i/R_r).^(2.*m) - J1_m.*R_m.*R_r.*Y1_r + J1_r.*R_m.*R_r.*Y1_m + J1_m.*R_m.*Y_r.*m - J_r.*R_m.*Y1_m.*m + J1_r.*R_r.*Y_m.*m - J_m.*R_r.*Y1_r.*m - J1_m.*R_m.*R_r.*Y1_r.*(R_i/R_r).^(2.*m) + J1_r.*R_m.*R_r.*Y1_m.*(R_i/R_r).^(2.*m) - J1_m.*R_m.*Y_r.*m.*(R_i/R_r).^(2.*m) + J_r.*R_m.*Y1_m.*m.*(R_i/R_r).^(2.*m) + J1_r.*R_r.*Y_m.*m.*(R_i/R_r).^(2.*m) - J_m.*R_r.*Y1_r.*m.*(R_i/R_r).^(2.*m) - J_m.*Y_r.*m.^2.*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) + J_r.*Y_m.*m.^2.*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) + J_m.*Y_r.*m.^2.*(R_i/R_r).^(2.*m).*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) - J_r.*Y_m.*m.^2.*(R_i/R_r).^(2.*m).*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) - J1_m.*R_m.*R_r.*Y1_r.*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) + J1_r.*R_m.*R_r.*Y1_m.*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) + J1_m.*R_m.*Y_r.*m.*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) - J_r.*R_m.*Y1_m.*m.*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) - J1_r.*R_r.*Y_m.*m.*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) + J_m.*R_r.*Y1_r.*m.*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) - J1_m.*R_m.*Y_r.*m.*(R_i/R_r).^(2.*m).*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) + J_r.*R_m.*Y1_m.*m.*(R_i/R_r).^(2.*m).*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) - J1_r.*R_r.*Y_m.*m.*(R_i/R_r).^(2.*m).*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) + J_m.*R_r.*Y1_r.*m.*(R_i/R_r).^(2.*m).*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) - J1_m.*R_m.*R_r.*Y1_r.*(R_i/R_r).^(2.*m).*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) + J1_r.*R_m.*R_r.*Y1_m.*(R_i/R_r).^(2.*m).*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m)));

C = ((R_r.*Y1_r - Y_r.*m + R_r.*Y1_r.*(R_i/R_r).^(2.*m) + Y_r.*m.*(R_i/R_r).^(2.*m)).*COMM_double);
D = -((J1_r.*R_r - J_r.*m + J1_r.*R_r.*(R_i/R_r).^(2.*m) + J_r.*m.*(R_i/R_r).^(2.*m)).*COMM_double);

A_pl = -(A_mJ.*mu_r.*(R_m/R_w).^(m - 1).*(mu_r + mu_r.*(R_i/R_r).^(2.*m) - (R_i/R_r).^(2.*m) + 1).*(2.*R_wi.^2.*(R_w/R_wi).^m - R_w.^2.*m - 2.*R_w.^2 - 2.*R_w.^2.*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) + R_wi.^2.*m.*(R_w/R_wi).^m + 2.*R_wi.^2.*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^m - R_wi.^2.*m.*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^m + R_w.^2.*m.*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m)))./(R_w.*(2.*mu_r - 2.*mu_r.*(R_r/R_m).^(2.*m) - (R_i/R_r).^(2.*m).*(R_r/R_m).^(2.*m) - (R_i/R_r).^(2.*m) + (R_r/R_m).^(2.*m) + mu_r.^2.*(R_i/R_r).^(2.*m) + mu_r.^2.*(R_r/R_m).^(2.*m) + mu_r.^2 + (R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) + mu_r.^2.*(R_i/R_r).^(2.*m).*(R_r/R_m).^(2.*m) - (R_i/R_r).^(2.*m).*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) + (R_r/R_m).^(2.*m).*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) - mu_r.^2.*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) + 2.*mu_r.*(R_i/R_r).^(2.*m).*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) - (R_i/R_r).^(2.*m).*(R_r/R_m).^(2.*m).*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) - mu_r.^2.*(R_i/R_r).^(2.*m).*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) - mu_r.^2.*(R_r/R_m).^(2.*m).*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) - 2.*mu_r.*(R_i/R_r).^(2.*m).*(R_r/R_m).^(2.*m).*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) - mu_r.^2.*(R_i/R_r).^(2.*m).*(R_r/R_m).^(2.*m).*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) + 1));
A_pl(H~=0) = 0;
A_mi = -(A_mJ.*mu_r.*(R_r/R_m).^(m - 1).*(R_m/R_w).^(m - 1).*(mu_r + mu_r.*(R_i/R_r).^(2.*m) + (R_i/R_r).^(2.*m) - 1).*(2.*R_wi.^2.*(R_w/R_wi).^m - R_w.^2.*m - 2.*R_w.^2 - 2.*R_w.^2.*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) + R_wi.^2.*m.*(R_w/R_wi).^m + 2.*R_wi.^2.*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^m - R_wi.^2.*m.*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^m + R_w.^2.*m.*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m)))./(R_w.*(2.*mu_r - 2.*mu_r.*(R_r/R_m).^(2.*m) - (R_i/R_r).^(2.*m).*(R_r/R_m).^(2.*m) - (R_i/R_r).^(2.*m) + (R_r/R_m).^(2.*m) + mu_r.^2.*(R_i/R_r).^(2.*m) + mu_r.^2.*(R_r/R_m).^(2.*m) + mu_r.^2 + (R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) + mu_r.^2.*(R_i/R_r).^(2.*m).*(R_r/R_m).^(2.*m) - (R_i/R_r).^(2.*m).*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) + (R_r/R_m).^(2.*m).*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) - mu_r.^2.*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) + 2.*mu_r.*(R_i/R_r).^(2.*m).*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) - (R_i/R_r).^(2.*m).*(R_r/R_m).^(2.*m).*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) - mu_r.^2.*(R_i/R_r).^(2.*m).*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) - mu_r.^2.*(R_r/R_m).^(2.*m).*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) - 2.*mu_r.*(R_i/R_r).^(2.*m).*(R_r/R_m).^(2.*m).*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) - mu_r.^2.*(R_i/R_r).^(2.*m).*(R_r/R_m).^(2.*m).*(R_m/R_w).^(2.*m).*(R_wi/R_s).^(2.*m).*(R_w/R_wi).^(2.*m) + 1));
A_mi(H~=0) = 0;
%% ------------------------------------------------------------------------
P_loss = conj((C.*J1_m+D.*Y1_m)./(mu_r*mu_0)); % H_theta* contribution to the loss equation
P_loss = (1i*H.*(p*speed_radsec).*(C.*J_m + D.*Y_m)).*P_loss;
P_loss = pi*R_m*l_a.*real(P_loss);
P_loss(isnan(P_loss)) = 0;
P_loss(isinf(P_loss)) = 0;
P_loss = sum(sum(P_loss));


%% PLOTS & MAPS
if contains(plotting, 'yes')
    T = 2*pi/(p*speed_radsec); % electrical period [s]
    shift = - 2*(t(2)-t(1));
    r = (R_m+R_r)/2; % test radius 
    % t = linspace (T,2*T,floor(N/2)+1);
    % t = t(round(length(t)/2):end);
    theta = 0;

    J_m = besselj(m,r*tau);
    Y_m = bessely(m,r*tau);
    J_m(isnan(J_m)) = 0;  
    Y_m(isnan(Y_m)) = 0; 


    Amp_Az = C.*J_m + D.*Y_m;
    Amp_Az(isnan(Amp_Az)) = 0;
    Amp_Az(isinf(Amp_Az)) = 0;
    Amp_Az = Amp_Az + (A_pl./m*R_m.*(r/R_m).^m + A_mi./m*R_r.*(R_r/r).^(m));


    Az = real(Amp_Az.'*(cos(H*p*speed_radsec.*t+phi+m.*(theta)) + 1i*sin(H*p*speed_radsec.*t+phi+m.*(theta))));
        figure
        plot(t,Az)

    J0 = real((1i*p*speed_radsec*cond*Amp_Az.*H).'*(cos(H*p*speed_radsec.*(t+shift)+phi+m.*(theta)) + 1i*sin(H*p*speed_radsec.*(t+shift)+phi+m.*(theta)))) ;
    figure
    plot(t,-J0)

end