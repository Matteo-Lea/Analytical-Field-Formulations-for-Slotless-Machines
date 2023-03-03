function [E_n,T_avg,t_spec,Magnets_weight,Joule,Tot_weight,Weighted_Perc_demag] = Field_Solution(el,m_PM,r_dis,x,params_geo,M_r_n,M_theta_n,params_mach,Const)
% Field_Solution (Field solution derived quantities
% -torque,back-emf,back-iron flux density-)
%   Detailed explanation goes here
mu_0 = 4*pi*1e-7; % air permeability

alpha_p = x(:,1); % mid-magnet to pole ratio [-]
l_a = x(:,2); % active lenght [m]
R_s = x(:,3);  % stator radius (facing winding) [m]
In_Out = x(:,8); % 1 for inrunner -1 for outrunner
R_wi = R_s; % winding radius (facing back-iron) [m]
R_w = R_wi - In_Out.*x(:,4); % winding radius (facing air-gap)[m]
R_m = R_w - In_Out.*x(:,5); % magnets array radius (facing air-gap)[m]
R_r = R_m - In_Out.*x(:,6); % magnets array radius (facing magnets support)[m]

% This definition of R_i shouldn't vary iteratively in the optimization
% process to keep the geometries with and without backing iron
% R_i = R_r; % infinite permeability boundary behind the magnets [m] (magnetic support is default)
% R_i(In_Out==1) = R_i(In_Out==1).*randi([0,1],nnz(In_Out==1),1); % non-magnetic support set for some in-runner topologies (randomly)
% R_i(In_Out==-1) = R_i(In_Out==-1).*1./randi([0,1],nnz(In_Out==-1),1); % non-magnetic support set for some out-runner topologies (randomly)
R_i = params_geo.R_i;
% R_sleeve = (x(:,12) - x(:,10).*x(:,16)- x(:,10).*x(:,15)- x(:,10).*x(:,14)- x(:,10).*x(:,17)).*(R_i~=0&R_i~=Inf); 
% R_sleeve = params_geo.R_sleeve;

B_r = Const(:,3);
B_rs = Const(:,4);
mu_r = params_mach.mu_r; % PM recoil permeability [-]

p = x(:,9); % pole pairs
% params_geo.In_Out = params_mag.In_Out; 



S_Ph = abs(pi*(R_wi.^2-R_w.^2)./(6*p)); % phase belt cross section
N_tc = params_mach.N_tc; % number of conductors per coil
q = params_mach.q ; %coils per pole and per phase
b = params_mach.b; % parallel ways per phase
% I = params_mach.I; % phase current peak value [A]
n_cs = N_tc*q/b;
omega = params_mach.omega*p; % electrical angular frequency [rad/s]
K_fill = params_mach.K_fill;
% J_Ph = n_cs*I./(S_Ph*K_fill); % peak copper current density [A/m2]
J_Ph = x(:,7); % peak copper current density [A/m2]
I = J_Ph.*(S_Ph*K_fill); % phase current peak value [A]
L_ew = 2.5*(R_w+R_wi)./p;
R_Ph = (params_mach.res_copp*1e-6)*q*2*p./3*n_cs.*(l_a+L_ew)./(S_Ph*K_fill);
B_sat = params_mach.B_sat;

mu_rc = 7.969; % PM coercive permeability (from knee-point to Hc)
H_c = 760*1e3; % PM hysteresis curve coercivity [A/m]
H_int = (B_r-mu_0*mu_rc*H_c)/(mu_0*(mu_rc-mu_r));
B_int = mu_0*mu_rc*(H_int+H_c);
% B_knee = B_r - abs(mu_0*mu_r*H_c);
B_knee = B_int;

Joule = 3*R_Ph.*I.^2;
%% useful indices for series harmonics definition
xi = 0:1:m_PM;
% x = linspace(0,40,m_PM+1);
n = p*(2*xi+1); % existing harmonic series components (n=1,3,5,...)


% %% circumferential discretization
% % tic
% m_th = 2000; % points along the circumference/arc
% theta = linspace(0,2*pi/(p),m_th)-pi/(2*p); % circumferential discretization (leave the -pi/(2*p))
% % dt = theta(2)-theta(1);
% % theta_mid = 
% Theta = repmat(theta,m_PM+1,1);


% HARMONICS COEFFICIENTS FROM THE MAGNETIC VECTOR POTENTIAL
% (WINDINGS/AIR-GAP REGIONS)

% coefficient from the particular solution of Poisson's equation accounting
% for the magnetization distribution
A_zm_n = mu_0*(n.*M_r_n+M_theta_n)./(n.^2-1);
% singularities are adjusted
A_zm_n(p==1,1) = -mu_0/2*(M_r_n(p==1,1)+M_theta_n(p==1,1));

% Coefficients from the field formulations
% Neither K_Bn nor K_Bn_out make any sense for np=1, the field solution
% exhibits a singularity for that very case which must be treated
% separately
DEN_I = ((mu_r.^2-1).*(-(R_i./R_r).^(2*n)-(R_m./R_s).^(2*n)+(R_r./R_s).^(2*n)+(R_i./R_m).^(2*n))+(mu_r+1).^2.*(1-(R_i./R_s).^(2*n))+(mu_r-1).^2.*((R_m./R_s).^(2*n).*(R_i./R_r).^(2*n)-(R_r./R_m).^(2*n)));

DEN_O = ((mu_r.^2-1).*(-(R_s./R_r).^(2*n)-(R_m./R_i).^(2*n)+(R_r./R_i).^(2*n)+(R_s./R_m).^(2*n))+(mu_r+1).^2.*((R_s./R_i).^(2*n)-1)+(mu_r-1).^2.*((R_m./R_r).^(2*n)-(R_r./R_i).^(2*n).*(R_s./R_m).^(2*n)));

K_Bn = (((-A_zm_n+n.*A_zm_n-mu_0*M_theta_n).*((mu_r+1)-(R_i./R_r).^(2*n).*(mu_r-1))+(A_zm_n+n.*A_zm_n+mu_0*M_theta_n).*((R_r./R_m).^(2*n).*(mu_r-1)-(R_i./R_m).^(2*n).*(mu_r+1))+2*(R_r./R_m).^(n+1).*(A_zm_n+mu_0*M_theta_n-mu_r.*n.*A_zm_n)+2*(R_r./R_m).^(n+1).*(R_i./R_r).^(2*n).*(A_zm_n+mu_0*M_theta_n+mu_r.*n.*A_zm_n))./DEN_I);

K_Bn_out = (((-A_zm_n+n.*A_zm_n-mu_0.*M_theta_n).*((R_m./R_i).^(2*n).*(mu_r+1)-(R_m./R_r).^(2*n).*(mu_r-1))+(A_zm_n+n.*A_zm_n+mu_0*M_theta_n).*((R_r./R_i).^(2*n).*(mu_r-1)-(mu_r+1))+2*(R_r./R_m).*(R_r./R_i).^(n).*(R_m./R_i).^(n).*(A_zm_n+mu_0*M_theta_n-mu_r.*n.*A_zm_n)+2*(R_m./R_r).^(n-1).*(A_zm_n+mu_0*M_theta_n+mu_r.*n.*A_zm_n))./DEN_O);

  % correction of singularities
   DEN_I_1 = -((R_s.^2+R_m.^2).*(R_m.^2-R_r.^2).*(R_r.^2+R_i.^2)-2.*(R_m.^2+R_r.^2).*(R_m.^2.*R_i.^2-R_s.^2.*R_r.^2).*mu_r+(R_s.^2-R_m.^2).*(R_m.^2-R_r.^2).*(R_r.^2-R_i.^2).*mu_r.^2);
   DEN_O_1 = DEN_I_1;
   K_Bn_1 = R_s.^2.*((R_m.^2-R_r.^2).*(A_zm_n(1)+mu_0*M_theta_n(1)).*(R_r.^2.*(mu_r+1)-R_i.^2.*(mu_r-1))+2.*R_r.^2.*A_zm_n(1).*(R_r.^2.*(mu_r-1)-R_i.^2.*(mu_r+1)).*log(R_r./R_m))./DEN_I_1;
   K_Bn_out_1 = K_Bn_1.*R_m.^2./R_s.^2;  
       DEN_O_1_inf = -((R_s.^2+R_m.^2).*(R_m.^2-R_r.^2)-2*(R_m.^2+R_r.^2).*R_m.^2.*mu_r-(R_s.^2-R_m.^2).*(R_m.^2-R_r.^2).*mu_r.^2);
       K_Bn_out_1_inf =  -R_m.^2.*((R_m.^2-R_r.^2).*(A_zm_n(1)+mu_0*M_theta_n(1)).*(mu_r-1)+2*R_r.^2.*A_zm_n(1).*(mu_r+1).*log(R_r./R_m))./DEN_O_1_inf;
    
K_Bn(p==1,1) = K_Bn_1(p==1,1);
K_Bn_out(p==1,1) = K_Bn_out_1(p==1,1);
K_Bn_out(p==1 & R_i==Inf,1) = K_Bn_out_1_inf(p==1 & R_i==Inf,1);



C_2 = R_m.*log(R_wi./R_w) + R_m./(4*R_s.^4).*(R_wi.^4-R_w.^4);
C = (R_m./R_s).^(n).*(R_wi.^2.*(R_wi./R_s).^(n)-R_w.^2.*(R_w./R_s).^(n))./((n+2).*R_m)+(R_wi.^2.*(R_m./R_wi).^(n)-R_w.^2.*(R_m./R_w).^(n))./((2-n).*R_m);
C(p==2,1) = C_2(p==2);

C_out_2 = 1./R_m.^3.*(R_s.^4.*log(R_w./R_wi)+R_w.^4/4-R_wi.^4/4); 
C_out = (R_w.^2.*(R_w./R_m).^(n)-R_wi.^2.*(R_wi./R_m).^(n))./((n+2).*R_m)+(R_w.^2./R_s.*(R_s./R_w).^(n)-R_wi.^2./R_s.*(R_s./R_wi).^(n))./(2-n).*(R_s./R_m).^(n+1);
C_out(p==2,1) = C_out_2(p==2);

% B_g_r_m = Amp_r_m*cos(n'.*Theta);

% Winding back-iron flux computaion
Flux_in = dot((K_Bn./(n).*R_m).*2.*(R_m./R_s).^(n),sin((n.*pi./(2*p))),2); % inrunner

Flux_out = dot((K_Bn_out./(n).*R_m).*2.*(R_s./R_m).^(n),sin((n.*pi./(2*p))),2); % outrunner


% B_bi = abs(Flux_in./(R_se-R_s)); % back-iron flux density (inrunner formulation)
% B_bi_out = abs(Flux_out./(R_se-R_s)); % back-iron flux density (outrunner formulation)
% B_bi(In_Out==-1) = B_bi_out(In_Out==-1); % back-iron flux density (outrunner adjustion)

t_bi = abs(Flux_in./(B_sat)); % back-iron flux density (inrunner formulation)
t_bi_out = abs(Flux_out./(B_sat)); % back-iron flux density (outrunner formulation)
t_bi(In_Out==-1) = t_bi_out(In_Out==-1); % back-iron flux density (outrunner adjustion)

% Back-emf and Torque computation

E_n = -(n_cs*l_a./S_Ph.*omega*4.*sin(n./p.*pi/6).*C.*(K_Bn.*R_m.^2./(n))); % back-emf harmonics (inrunner formulation)
E_n_out = -(n_cs*l_a./S_Ph.*omega*4.*sin(n./p.*pi/6).*C_out.*(K_Bn_out.*R_m.^2./(n))); % back-emf harmonics (outrunner formulation)
E_n(In_Out==-1,:) = E_n_out(In_Out==-1,:); % back-emf harmonics (outrunner adjustion)

% emf_a = E_n'.*cos(n'.*Theta);
% emf_b = E_n'.*cos(n'.*Theta-2*pi/3*(n/p)');
% emf_c = E_n'.*cos(n'.*Theta-4*pi/3*(n/p)');

% resulting radial component along the circumference
% emf_a(isnan(emf_a))=0; % higher order harmonics are not represented in Matlab double precision floating point representation (NaN numbers are therefore zeroed out)
% emf_a = sum(emf_a,1);
% emf_b(isnan(emf_b))=0; % higher order harmonics are not represented in Matlab double precision floating point representation (NaN numbers are therefore zeroed out)
% emf_b = sum(emf_b,1);
% emf_c(isnan(emf_c))=0; % higher order harmonics are not represented in Matlab double precision floating point representation (NaN numbers are therefore zeroed out)
% emf_c = sum(emf_c,1);
% 
% T = -p/omega*(emf_a*I.*cos(p*theta+theta_i)+emf_b*I.*cos(p*theta-2*pi/3+theta_i)+emf_c*I.*cos(p*theta-4*pi/3+theta_i));

T_avg = -p./omega*3/2.*E_n(:,1).*I;


m_th = 100; % points along the circumference/arc
theta = (linspace(0,pi,m_th)-pi/(2))./p; % circumferential discretization (leave the -pi/(2*p))
theta = reshape(theta,1,[],el);
sigma = (sin(pi*n./n(:,end))./(pi*n./n(:,end))).^1; % Lanczos sigma for Gibbs phenomenon reduction
% (Amp_m_r_m)*cos(reshape(n,[],1,el).*theta);
% ADAPTED FIELD SOLUTION IN THE MAGNETS REGION WITH MAGNETIC SUPPORT TO
% PREVENT THE MAGNETIC SLEEVE FROM SATURATUING
% FLUX DENSITY IN THE MAGNETS REGION

A_m_pl = (-((mu_r+1)-(R_i./R_r).^(2*n).*(mu_r-1)).*((mu_0*M_theta_n+A_zm_n+mu_r.*n.*A_zm_n)+(R_m./R_s).^(2*n).*(mu_0*M_theta_n+A_zm_n-mu_r.*n.*A_zm_n))+(R_r./R_m).^(n+1).*((R_m./R_s).^(2*n).*(mu_r+1)-(mu_r-1)).*((mu_0*M_theta_n+A_zm_n-mu_r.*n.*A_zm_n)+(R_i./R_r).^(2*n).*(mu_0*M_theta_n+A_zm_n+mu_r.*n.*A_zm_n)))./DEN_I;

A_m_mi = (-(R_r./R_m).^(n-1).*(-(mu_r-1)+(R_i./R_r).^(2*n).*(mu_r+1)).*((mu_0*M_theta_n+A_zm_n+mu_r.*n.*A_zm_n)+(R_m./R_s).^(2*n).*(mu_0*M_theta_n+A_zm_n-mu_r.*n.*A_zm_n))-((R_m./R_s).^(2*n).*(mu_r-1)-(mu_r+1)).*((mu_0*M_theta_n+A_zm_n-mu_r.*n.*A_zm_n)+(R_i./R_r).^(2*n).*(mu_0*M_theta_n+A_zm_n+mu_r.*n.*A_zm_n)))./DEN_I;

Flux_in_sl = dot(((R_m.*A_m_pl.*(R_r./R_m).^n+R_r.*A_m_mi)./n+A_zm_n.*R_r),sin((n.*pi./(2*p))),2); % inrunner


r_m = R_r + linspace(0,1,r_dis).*(R_m-R_r);
r_m = permute(shiftdim(r_m,-1),[2 1 3]);

Amp_m_r_m = permute(sigma.*(A_m_pl.*(r_m./R_m).^(n-1)+A_m_mi.*(R_r./r_m).^(n+1)+repmat(n.*A_zm_n,1,1,r_dis)),[3 2 1]);
Amp_m_theta_m = permute(sigma.*(A_m_pl.*(r_m./R_m).^(n-1)-A_m_mi.*(R_r./r_m).^(n+1)+repmat(A_zm_n,1,1,r_dis)),[3 2 1]);

theta_1 = (linspace(0,pi,m_th)-pi/(2)).*alpha_p./p; % circumferential discretization mid-magnet (leave the -pi/(2))
theta_2 = -pi./(2*p) + linspace(0,1,m_th).*(-alpha_p*pi./(2*p)+pi./(2*p)); % circumferential discretization left-side-magnet (leave the -pi/(2))
theta_3 = alpha_p*pi./(2*p) + linspace(0,1,m_th).*(pi./(2*p)-alpha_p*pi./(2*p)); % circumferential discretization right-side-magnet (leave the -pi/(2))
Theta_1 = shiftdim(theta_1',-1);
Theta_2 = shiftdim(theta_2',-1);
Theta_3 = shiftdim(theta_3',-1);
Theta_1_mag = repmat(Theta_1,r_dis,1,1);
Theta_2_mag = repmat(Theta_2,r_dis,1,1)+shiftdim(pi./(2*p)',-1)+pi/2;
Theta_3_mag = repmat(Theta_3,r_dis,1,1)-shiftdim(pi./(2*p)',-1)-pi/2;
B_r_1 = pagemtimes(Amp_m_r_m,cos(permute(shiftdim(n',-1),[2 1 3]).*Theta_1));
B_theta_1 = pagemtimes(-Amp_m_theta_m,sin(permute(shiftdim(n',-1),[2 1 3]).*Theta_1));
B_r_2 = pagemtimes(Amp_m_r_m,cos(permute(shiftdim(n',-1),[2 1 3]).*Theta_2));
B_theta_2 = pagemtimes(-Amp_m_theta_m,sin(permute(shiftdim(n',-1),[2 1 3]).*Theta_2));
B_r_3 = pagemtimes(Amp_m_r_m,cos(permute(shiftdim(n',-1),[2 1 3]).*Theta_3));
B_theta_3 = pagemtimes(-Amp_m_theta_m,sin(permute(shiftdim(n',-1),[2 1 3]).*Theta_3));

B_mag_1 = B_r_1.*cos(Theta_1_mag)+B_theta_1.*sin(Theta_1_mag);
B_mag_2 = B_r_2.*cos(Theta_2_mag)+B_theta_2.*sin(Theta_2_mag);
B_mag_3 = B_r_3.*cos(Theta_3_mag)+B_theta_3.*sin(Theta_3_mag);

B_knee = B_knee(1);
B_dem_1 = B_knee - B_mag_1;
B_dem_2 = B_knee - B_mag_2;
B_dem_3 = B_knee - B_mag_3;
B_dem_1(B_dem_1 <= 0) = 0;
B_dem_2(B_dem_2 <= 0) = 0;
B_dem_3(B_dem_3 <= 0) = 0;

Weight = [B_dem_1./B_knee B_dem_2./B_knee B_dem_3./B_knee];

Weighted_Perc_demag = sum(sum(Weight))./(3*m_th*r_dis)*100;
Weighted_Perc_demag = squeeze(Weighted_Perc_demag);
Perc_demag = sum(squeeze(sum(Weight~=0,1)))./(3*m_th*r_dis)*100;


% tic
% Weighted_Perc_demag_loop(el) = 0;
% Perc_demag_loop(el) = 0;
% for ii = 1:el
%     r_m_loop = linspace(R_r(ii),R_m(ii),r_dis)';
%     Amp_m_r_m_loop = sigma(ii,:).*(A_m_pl(ii,:).*(r_m_loop./R_m(ii)).^(n(ii,:)-1)+A_m_mi(ii,:).*(R_r(ii)./r_m_loop).^(n(ii,:)+1)+n(ii,:).*A_zm_n(ii,:));
%     Amp_m_theta_m_loop = sigma(ii,:).*(A_m_pl(ii,:).*(r_m_loop./R_m(ii)).^(n(ii,:)-1)-A_m_mi(ii,:).*(R_r(ii)./r_m_loop).^(n(ii,:)+1)+A_zm_n(ii,:));
%     theta_1_loop = (linspace(0,pi,m_th)-pi/(2)).*alpha_p(ii)./p(ii); % circumferential discretization mid-magnet (leave the -pi/(2))
%     theta_2_loop = linspace(-pi./(2*p(ii)),-alpha_p(ii)*pi./(2*p(ii)),m_th); % circumferential discretization left-side-magnet (leave the -pi/(2))
%     theta_3_loop = linspace(alpha_p(ii)*pi./(2*p(ii)),pi./(2*p(ii)),m_th); % circumferential discretization right-side-magnet (leave the -pi/(2))
%     B_r_1_loop = Amp_m_r_m_loop*cos(n(ii,:).*theta_1_loop')';
%     B_theta_1_loop = -Amp_m_theta_m_loop*sin(n(ii,:).*theta_1_loop')';
%     B_r_2_loop = Amp_m_r_m_loop*cos(n(ii,:).*theta_2_loop')';
%     B_theta_2_loop = -Amp_m_theta_m_loop*sin(n(ii,:).*theta_2_loop')';
%     B_r_3_loop = Amp_m_r_m_loop*cos(n(ii,:).*theta_3_loop')';
%     B_theta_3_loop = -Amp_m_theta_m_loop*sin(n(ii,:).*theta_3_loop')';
%     
%     B_mag_1_loop = B_r_1_loop.*cos(theta_1_loop)+B_theta_1_loop.*sin(theta_1_loop);
%     B_mag_2_loop = B_r_2_loop.*cos(theta_2_loop+pi/(2*p(ii))+pi/2)+B_theta_2_loop.*sin(theta_2_loop+pi/(2*p(ii))+pi/2);
%     B_mag_3_loop = B_r_3_loop.*cos(theta_3_loop-pi/(2*p(ii))-pi/2)+B_theta_3_loop.*sin(theta_3_loop-pi/(2*p(ii))-pi/2);
%     
%     B_dem_1_loop = B_knee - B_mag_1_loop;
%     B_dem_1_loop(B_dem_1_loop <= 0) = 0;
%     B_dem_2_loop = B_knee - B_mag_2_loop;
%     B_dem_2_loop(B_dem_2_loop <= 0) = 0;
%     B_dem_3_loop = B_knee - B_mag_3_loop;
%     B_dem_3_loop(B_dem_3_loop <= 0) = 0;
%     
%     Weight_loop = [B_dem_1_loop./B_knee B_dem_2_loop./B_knee B_dem_3_loop./B_knee];
%     Weighted_Perc_demag_loop(ii) = sum(sum(Weight_loop))./(3*m_th*r_dis)*100;
%     Perc_demag_loop(ii) = nnz(Weight_loop)/(3*m_th*r_dis)*100;
%     
% end
% toc


% FLUX DENSITY IN THE MAGNETS REGION OUTRUNNER
A_m_pl = (-(-(mu_r+1)+(R_s./R_m).^(2*n).*(mu_r-1)).*((mu_0*M_theta_n+A_zm_n+mu_r.*n.*A_zm_n)+(R_r./R_i).^(2*n).*(mu_0*M_theta_n+A_zm_n-mu_r.*n.*A_zm_n))-(R_m./R_r).^(n+1).*((R_r./R_i).^(2*n).*(mu_r+1)-(mu_r-1)).*((mu_0*M_theta_n+A_zm_n-mu_r.*n.*A_zm_n)+(R_s./R_m).^(2*n).*(mu_0*M_theta_n+A_zm_n+mu_r.*n.*A_zm_n)))./DEN_O;

A_m_mi = ((R_m./R_r).^(n-1).*(-(mu_r-1)+(R_s./R_m).^(2*n).*(mu_r+1)).*((mu_0*M_theta_n+A_zm_n+mu_r.*n.*A_zm_n)+(R_r./R_i).^(2*n).*(mu_0*M_theta_n+A_zm_n-mu_r.*n.*A_zm_n))+((R_r./R_i).^(2*n).*(mu_r-1)-(mu_r+1)).*((mu_0*M_theta_n+A_zm_n-mu_r.*n.*A_zm_n)+(R_s./R_m).^(2*n).*(mu_0*M_theta_n+A_zm_n+mu_r.*n.*A_zm_n)))./DEN_O;

Flux_out_sl = dot(((R_r.*A_m_pl+R_m.*A_m_mi.*(R_m./R_r).^n)./n+A_zm_n.*R_r),sin((n.*pi./(2*p))),2); % outrunner

r_m = R_m + linspace(0,1,r_dis).*(R_r-R_m);

% B_sleeve = abs(Flux_in_sl./(R_r-R_sleeve)); % back-iron flux density (inrunner formulation)
% B_sleeve_out = abs(Flux_out_sl./(R_r-R_sleeve)); % back-iron flux density (outrunner formulation)
% B_sleeve(In_Out==-1) = B_sleeve_out(In_Out==-1); % back-iron flux density (outrunner adjustion)

t_sleeve = abs(Flux_in_sl./(B_sat)); % back-iron flux density (inrunner formulation)
t_sleeve_out = abs(Flux_out_sl./(B_sat)); % back-iron flux density (outrunner formulation)
t_sleeve(In_Out==-1) = t_sleeve_out(In_Out==-1); % back-iron flux density (outrunner adjustion)
t_sleeve(t_sleeve<2e-3) = 2e-3;

% the case pn=1 is adjusted in the following, as the general solution
% exhibits a singularity in that very case
% A_z_1_pl = (A_zm_n(:,1)+mu_0*M_theta_n(:,1)).*(R_m.^4+R_s.^2.*(R_m.^2+R_r.^2.*(mu_r-1))-R_r.^2.*R_m.^2.*(mu_r+1))+R_r.^2.*A_zm_n(:,1).*(R_s.^2.*(mu_r-1)-R_m.^2.*(mu_r+1)).*log(R_r)+A_zm_n(:,1).*R_m.^2.*(-R_m.^2.*(mu_r-1)+R_s.^2.*(1+mu_r)).*log(R_m)./(R_m.^4.*(mu_r-1)+R_r.^2.*R_m.^2.*(mu_r+1)-R_s.^2.*(R_r.^2.*(mu_r-1)+R_m.^2.*(mu_r+1)));
% A_z_1_mi = R_r.^2.*R_m.^2.*((R_m.^2-R_s.^2).*(A_zm_n(:,1)+mu_0*M_theta_n(:,1)).*mu_r+A_zm_n(:,1).*(R_m.^2.*(mu_r-1)-R_s.^2.*(1+mu_r)).*log(R_r)-A_zm_n(:,1).*R_m.^2.*(mu_r-1).*log(R_m)+A_zm_n(:,1).*R_s.^2.*(mu_r+1).*log(R_m))./(R_m.^4.*(mu_r-1)+R_r.^2.*R_m.^2.*(mu_r+1)-R_s.^2.*(R_r.^2.*(mu_r-1)+R_m.^2.*(mu_r+1)));
A_z_1_pl = ((A_zm_n(:,1)+mu_0.*M_theta_n(:,1)).*(-R_m.^2.*R_r.^2.*(R_r.^2+R_i.^2).*(1+mu_r)+R_m.^4.*(R_r.^2.*(mu_r+1)-R_i.^2.*(mu_r-1))+R_s.^2.*(R_r.^2.*(R_r.^2+R_i.^2).*(mu_r-1)+R_m.^2.*(R_r.^2.*(mu_r+1)-R_i.^2.*(mu_r-1))))+R_m.^2.*A_zm_n(:,1).*(-R_m.^2.*(mu_r-1)+R_s.^2.*(mu_r+1)).*(-R_i.^2.*(mu_r-1)+R_r.^2.*(mu_r+1)).*log(R_m)-R_r.^2.*A_zm_n(:,1).*(R_s.^2.*(mu_r-1)-R_m.^2.*(mu_r+1)).*(R_r.^2.*(mu_r-1)-R_i.^2.*(mu_r+1)).*log(R_r))./(R_m.^2.*(R_r.^2.*R_i.^2.*(mu_r+1).^2-R_r.^4.*(mu_r.^2-1)+R_m.^2.*(mu_r-1).*(R_r.^2.*(mu_r+1)-R_i.^2.*(mu_r-1)))+R_s.^2.*(-R_m.^2.*(mu_r+1).*(R_r.^2.*(mu_r+1)-R_i.^2.*(mu_r-1))+R_r.^2.*(mu_r-1).*(R_r.^2.*(mu_r-1)-R_i.^2.*(mu_r+1))));
A_z_1_mi = -(R_m.^2.*R_r.^2.*(2.*(R_s.^2.*R_r.^2-R_m.^2.*R_i.^2).*(A_zm_n(:,1)+mu_0.*M_theta_n(:,1)).*mu_r+A_zm_n(:,1).*(-R_m.^2.*(mu_r-1)+R_s.^2.*(mu_r+1)).*(R_r.^2.*(mu_r-1)-R_i.^2.*(mu_r+1)).*log(R_m)-A_zm_n(:,1).*(-R_m.^2.*(mu_r-1)+R_s.^2.*(mu_r+1)).*(R_r.^2.*(mu_r-1)-R_i.^2.*(mu_r+1)).*log(R_r)))/(R_m.^2.*(R_r.^2.*R_i.^2.*(mu_r+1).^2-R_r.^4.*(mu_r.^2-1)+R_m.^2.*(mu_r-1).*(R_r.^2.*(mu_r+1)-R_i.^2.*(mu_r-1)))+R_s.^2.*(-R_m.^2.*(mu_r+1).*(R_r.^2.*(mu_r+1)-R_i.^2.*(mu_r-1))+R_r.^2.*(mu_r-1).*(R_r.^2.*(mu_r-1)-R_i.^2.*(mu_r+1))));

% A_z_1_pl_inf = ((A_zm_n(:,1)+mu_0.*M_theta_n(:,1)).*(-R_m.^2.*R_r.^2.*(1+mu_r)+R_m.^4.*(mu_r-1)+R_s.^2.*(R_r.^2.*(mu_r-1)-R_m.^2.*(mu_r-1)))-R_m.^2.*A_zm_n(:,1).*(-R_m.^2.*(mu_r-1)+R_s.^2.*(mu_r+1)).*(mu_r-1).*log(R_m)+R_r.^2.*A_zm_n(:,1).*(R_s.^2.*(mu_r-1)-R_m.^2.*(mu_r+1)).*(mu_r+1).*log(R_r))./(R_m.^2.*(R_r.^2.*(mu_r+1).^2-R_m.^2.*(mu_r-1).^2)+R_s.^2.*(R_m.^2.*(mu_r+1).*(mu_r-1)-R_r.^2.*(mu_r-1).*(mu_r+1)));
% A_z_1_mi_inf = (R_m.^2.*R_r.^2.*(2.*R_m.^2.*(A_zm_n(:,1)+mu_0.*M_theta_n(:,1)).*mu_r+A_zm_n(:,1).*(-R_m.^2.*(mu_r-1)+R_s.^2.*(mu_r+1)).*(mu_r+1).*log(R_m)-A_zm_n(:,1).*(-R_m.^2.*(mu_r-1)+R_s.^2.*(mu_r+1)).*(mu_r+1).*log(R_r)))./(R_m.^2.*(R_r.^2.*(mu_r+1).^2-R_m.^2.*(mu_r-1).^2)+R_s.^2.*(R_m.^2.*(mu_r+1).*(mu_r-1)-R_r.^2.*(mu_r-1).*(mu_r+1)));
 

Amp_Az_1 = A_z_1_pl.*R_r+A_z_1_mi./R_r+A_zm_n(:,1).*R_r.*log(R_r);
% Amp_Az_1_inf = A_z_1_pl_inf.*R_r+A_z_1_mi./R_r+A_zm_n(:,1).*R_r.*log(R_r);
t_sleeve(p==1) = abs((Amp_Az_1(p==1,1)*sin(pi/2))./(B_sat));

% B_sleeve(R_sleeve==0) = B_sat; % B_sleeve is not relevant if no magnetic support is there

rho_mag = params_mach.rho_mag; % PM density [kg/m3] (NdFeB=7500;Ferrite=5000;BondedNdFeB=5100;SmCo=8400;Alnico=7300)
rho_iron = params_mach.rho_iron; % iron density [kg/m3] 
rho_copper = params_mach.rho_copper; % copper density [kg/m3] 
rho_epoxy = params_mach.rho_epoxy; % epoxy density [kg/m3] 
rho_alu = params_mach.rho_alu; % aluminum density [kg/m3]
K_fill = params_mach.K_fill; % copper filling factor [-]

Winding_weight = pi*abs(R_wi.^2-R_w.^2).*(l_a+L_ew)*(K_fill*rho_copper + (1-K_fill)*rho_epoxy);
Stat_bi_weight = pi*abs(R_s.^2-In_Out.*(R_s+t_bi).^2).*l_a*rho_iron;
Rot_sleeve_weight = (pi*abs(R_r.^2-In_Out.*(R_r-t_sleeve).^2).*l_a*rho_iron).*(R_i~=0&R_i~=Inf);
Rot_sleeve_weight = Rot_sleeve_weight + (pi*abs(R_r.^2-In_Out.*(R_r-t_sleeve).^2).*l_a*rho_alu).*(R_i==0|R_i==Inf);
Magnets_weight = (pi*abs(R_r.^2-R_m.^2).*l_a*rho_mag);

Tot_weight = Winding_weight + Stat_bi_weight + Rot_sleeve_weight + Magnets_weight;

t_spec = T_avg./Tot_weight; % specific torque [Nm/kg]

end

