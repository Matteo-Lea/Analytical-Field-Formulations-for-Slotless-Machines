function [t_bi,peak,T_avg,t_spec,Magnets_weight,Joule,Tot_weight] = Field_Solution(el,m_PM,x,params_geo,M_r_n,M_theta_n,params_mach)
% Field_Solution (Field solution derived quantities
% -torque,back-emf,back-iron flux density-)
%   Detailed explanation goes here
mu_0 = 4*pi*1e-7; % air permeability

R_s = x(:,3);  % stator radius (facing winding) [m]
l_a = x(:,2).*R_s; % active lenght [m]
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
% R_i = params_geo.R_i;
R_i = 0;
% R_sleeve = (x(:,12) - x(:,10).*x(:,16)- x(:,10).*x(:,15)- x(:,10).*x(:,14)- x(:,10).*x(:,17)).*(R_i~=0&R_i~=Inf); 
% R_sleeve = params_geo.R_sleeve;


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

m_th = 300; % points along the circumference/arc

Joule = 3*R_Ph.*I.^2;
%% useful indices for series harmonics definition 
xi = 0:1:m_PM;

n = p*(2*xi+1); % existing harmonic series components (n=1,3,5,...)

%% (INRUNNER ONLY IMPLEMENTED) CURRENT FIELD CONTRIBUTION
m_J = 4;                                                                  % the total number of harmonics will be 2*(m_J+1)
z = 0:1:m_J;
y = z+1;
% current harmonics order
h = p*(6*z+1);
k = p*(6*y-1);
m = [h,k];                                                                  % the harmponics are put together in a single vector
 

% Current density distribution
I_tot = N_tc*q*I/b;                                                         % total current through the phase belt
J_Ph = I_tot./S_Ph;                                                          % phase current density
J_coeff = 4*3*p.*J_Ph./(2*pi);                                                % constant coefficient for current density distribution

J_h = (J_coeff.*sin(pi*h./(6*p))./(h));
J_k = (J_coeff.*sin(pi*k./(6*p))./(k));

J_m = [J_h,J_k];

A_mJ = mu_0*J_m./(m.^2-4);


% if R_s>R_m %  inrunner
    DEN_JI = 2*((R_i./R_s).^(2*m)-1);


    A_pl_J_w = A_mJ.*(R_wi.*(2-(R_wi./R_s).^(2*m).*(m-2)+m)-R_w.*(R_wi./R_s).^(2*m).*(R_w./R_wi).^(m+1).*(2-m+(R_i./R_w).^(2*m).*(2+m)))./DEN_JI;

    A_mi_J_w = A_mJ.*(R_wi.*(R_w./R_wi).^(m-1).*(R_i./R_w).^(2*m).*(2-(R_wi./R_s).^(2*m).*(m-2)+m)-R_w.*(2-m+(R_i./R_w).^(2*m).*(2+m)))./DEN_JI;

%     if p == 2
%        A_pl_J_w(1) = R_wi*A_mJ(1)*(R_w^4+R_i^4-R_s^4-R_wi^4-4*R_s^4*log(R_wi)+4*R_i^4*log(R_w))/DEN_JI(1);
%        A_mi_J_w(1) = A_mJ(1)/R_w^3*(R_s^4*R_w^4-R_wi^4*R_i^4+4*R_s^4*R_i^4*log(R_w/R_wi))/DEN_JI(1);
%     end

    % one-component magnetic potential 
    Amp_Az_J_w = (A_pl_J_w./m.*R_wi+A_mi_J_w./m.*R_w.*(R_w./R_s).^(m)+A_mJ.*R_s.^2);
%     if p == 2
%         Amp_Az_J_w(:,1) = (A_pl_J_w(1)./2*R_wi.*(r_w./R_wi).^(2)+A_mi_J_w(1)./2*R_w.*(R_w./r_w).^(2)+A_mJ(1).*r_w.^2.*log(r_w));
%     end
%     A_z_J_w = sigma_m.*Amp_Az_J_w*cos(m'.*Theta);
    flux_J = dot(2*Amp_Az_J_w,ones(el,2*(m_J+1)),2); % inrunner







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

% WRONG TORQUE FOR PSO VALIDATION
% Amp_r_m = (K_Bn.*((R_w./R_s).^(n-1).*(R_m./R_s).^(n+1)+(R_m./R_w).^(n+1)));
% kp = sin(pi./p./6)./(pi./p./6);
% 
% TORQUE = 3*Amp_r_m(:,1).*p.*l_a.*((R_w+R_s)./2).^2.*kp.*n_cs.*I./2;


C_2 = R_m.*log(R_wi./R_w) + R_m./(4*R_s.^4).*(R_wi.^4-R_w.^4);
C = (R_m./R_s).^(n).*(R_wi.^2.*(R_wi./R_s).^(n)-R_w.^2.*(R_w./R_s).^(n))./((n+2).*R_m)+(R_wi.^2.*(R_m./R_wi).^(n)-R_w.^2.*(R_m./R_w).^(n))./((2-n).*R_m);
C(p==2,1) = C_2(p==2);

C_out_2 = 1./R_m.^3.*(R_s.^4.*log(R_w./R_wi)+R_w.^4/4-R_wi.^4/4); 
C_out = (R_w.^2.*(R_w./R_m).^(n)-R_wi.^2.*(R_wi./R_m).^(n))./((n+2).*R_m)+(R_w.^2./R_s.*(R_s./R_w).^(n)-R_wi.^2./R_s.*(R_s./R_wi).^(n))./(2-n).*(R_s./R_m).^(n+1);
C_out(p==2,1) = C_out_2(p==2);

% B_g_r_m = Amp_r_m*cos(n'.*Theta);

% Winding back-iron flux computaion
% ((K_Bn./(n)*R_m).*((r./R_s).^(n).*(R_m/R_s).^(n)+(R_m./r).^(n)));
% Amp_Az_J_w = (A_pl_J_w./m*R_wi.*(r_w./R_wi).^(m)+A_mi_J_w./m*R_w.*(R_w./r_w).^(m)+A_mJ.*r_w.^2);
Flux_in = dot((K_Bn./(n).*R_m).*2.*(R_m./R_s).^(n),sin((n.*pi./(2*p))),2); % inrunner

Flux_out = dot((K_Bn_out./(n).*R_m).*2.*(R_s./R_m).^(n),sin((n.*pi./(2*p))),2); % outrunner


% B_bi = abs(Flux_in./(R_se-R_s)); % back-iron flux density (inrunner formulation)
% B_bi_out = abs(Flux_out./(R_se-R_s)); % back-iron flux density (outrunner formulation)
% B_bi(In_Out==-1) = B_bi_out(In_Out==-1); % back-iron flux density (outrunner adjustion)

% t_bi = abs(Flux_in./(B_sat)); % back-iron flux density (inrunner formulation)
t_bi = abs((Flux_in+flux_J)./(B_sat)); % back-iron flux density (inrunner formulation)
t_bi_out = abs(Flux_out./(B_sat)); % back-iron flux density (outrunner formulation)
t_bi(In_Out==-1) = t_bi_out(In_Out==-1); % back-iron flux density (outrunner adjustion)

% Back-emf and Torque computation

E_n = -(n_cs*l_a./S_Ph.*omega*4.*sin(n./p.*pi/6).*C.*(K_Bn.*R_m.^2./(n))); % back-emf harmonics (inrunner formulation)
E_n_out = -(n_cs*l_a./S_Ph.*omega*4.*sin(n./p.*pi/6).*C_out.*(K_Bn_out.*R_m.^2./(n))); % back-emf harmonics (outrunner formulation)
E_n(In_Out==-1,:) = E_n_out(In_Out==-1,:); % back-emf harmonics (outrunner adjustion)

theta_1 = (linspace(-pi/(2),pi./(2),m_th))./p; % circumferential discretization of a half mid-magnet
Theta_1 = shiftdim(theta_1',-1);
peak = squeeze(pagemtimes(permute(E_n,[3 2 1]),cos(permute(shiftdim(n',-1),[2 1 3]).*Theta_1)));
peak = max(abs(peak))';

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



% ADAPTED FIELD SOLUTION IN THE MAGNETS REGION WITH MAGNETIC SUPPORT TO
% PREVENT THE MAGNETIC SLEEVE FROM SATURATUING
% FLUX DENSITY IN THE MAGNETS REGION

A_m_pl = (-((mu_r+1)-(R_i./R_r).^(2*n).*(mu_r-1)).*((mu_0*M_theta_n+A_zm_n+mu_r.*n.*A_zm_n)+(R_m./R_s).^(2*n).*(mu_0*M_theta_n+A_zm_n-mu_r.*n.*A_zm_n))+(R_r./R_m).^(n+1).*((R_m./R_s).^(2*n).*(mu_r+1)-(mu_r-1)).*((mu_0*M_theta_n+A_zm_n-mu_r.*n.*A_zm_n)+(R_i./R_r).^(2*n).*(mu_0*M_theta_n+A_zm_n+mu_r.*n.*A_zm_n)))./DEN_I;

A_m_mi = (-(R_r./R_m).^(n-1).*(-(mu_r-1)+(R_i./R_r).^(2*n).*(mu_r+1)).*((mu_0*M_theta_n+A_zm_n+mu_r.*n.*A_zm_n)+(R_m./R_s).^(2*n).*(mu_0*M_theta_n+A_zm_n-mu_r.*n.*A_zm_n))-((R_m./R_s).^(2*n).*(mu_r-1)-(mu_r+1)).*((mu_0*M_theta_n+A_zm_n-mu_r.*n.*A_zm_n)+(R_i./R_r).^(2*n).*(mu_0*M_theta_n+A_zm_n+mu_r.*n.*A_zm_n)))./DEN_I;

Flux_in_sl = dot(((R_m.*A_m_pl.*(R_r./R_m).^n+R_r.*A_m_mi)./n+A_zm_n.*R_r),sin((n.*pi./(2*p))),2); % inrunner

% Amp_Az_m = ((R_r*A_m_pl+R_m*A_m_mi.*(R_m./R_r).^n)./n+A_zm_n.*R_r);


% FLUX DENSITY IN THE MAGNETS REGION
A_m_pl = (-(-(mu_r+1)+(R_s./R_m).^(2*n).*(mu_r-1)).*((mu_0*M_theta_n+A_zm_n+mu_r.*n.*A_zm_n)+(R_r./R_i).^(2*n).*(mu_0*M_theta_n+A_zm_n-mu_r.*n.*A_zm_n))-(R_m./R_r).^(n+1).*((R_r./R_i).^(2*n).*(mu_r+1)-(mu_r-1)).*((mu_0*M_theta_n+A_zm_n-mu_r.*n.*A_zm_n)+(R_s./R_m).^(2*n).*(mu_0*M_theta_n+A_zm_n+mu_r.*n.*A_zm_n)))./DEN_O;

A_m_mi = ((R_m./R_r).^(n-1).*(-(mu_r-1)+(R_s./R_m).^(2*n).*(mu_r+1)).*((mu_0*M_theta_n+A_zm_n+mu_r.*n.*A_zm_n)+(R_r./R_i).^(2*n).*(mu_0*M_theta_n+A_zm_n-mu_r.*n.*A_zm_n))+((R_r./R_i).^(2*n).*(mu_r-1)-(mu_r+1)).*((mu_0*M_theta_n+A_zm_n-mu_r.*n.*A_zm_n)+(R_s./R_m).^(2*n).*(mu_0*M_theta_n+A_zm_n+mu_r.*n.*A_zm_n)))./DEN_O;

Flux_out_sl = dot(((R_r.*A_m_pl+R_m.*A_m_mi.*(R_m./R_r).^n)./n+A_zm_n.*R_r),sin((n.*pi./(2*p))),2); % outrunner

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

% wrong assumption
% UNO = ((K_Bn./(n).*R_m).*((R_w./R_s).^(n).*(R_m./R_s).^(n)+(R_m./R_w).^(n)))*2;
% UNO = dot(UNO,sin((n.*pi./(2*p))),2); % inrunner
% DUE = (A_pl_J_w./m.*R_wi.*(R_w./R_wi).^(m)+A_mi_J_w./m.*R_w+A_mJ.*R_w.^2);
% DUE = dot(2*DUE,ones(el,2*(m_J+1)),2); % inrunner
% t_bi_new = abs((UNO+DUE)./(B_sat)); % back-iron flux density (inrunner formulation)
% Stat_bi_weight = pi*abs(R_s.^2-In_Out.*(R_s+t_bi_new).^2).*l_a*rho_iron;


Winding_weight = pi*abs(R_wi.^2-R_w.^2).*(l_a+L_ew)*(K_fill*rho_copper + (1-K_fill)*rho_epoxy);
Stat_bi_weight = pi*abs(R_s.^2-In_Out.*(R_s+t_bi).^2).*l_a*rho_iron;
Rot_sleeve_weight = (pi*abs(R_r.^2-In_Out.*(R_r-t_sleeve).^2).*l_a*rho_iron).*(R_i~=0&R_i~=Inf);
Rot_sleeve_weight = Rot_sleeve_weight + (pi*abs(R_r.^2-In_Out.*(R_r-t_sleeve).^2).*l_a*rho_alu).*(R_i==0|R_i==Inf);
Magnets_weight = (pi*abs(R_r.^2-R_m.^2).*l_a*rho_mag);

Tot_weight = Winding_weight + Stat_bi_weight + Rot_sleeve_weight + Magnets_weight;

t_spec = T_avg./Tot_weight; % specific torque [Nm/kg]

end

