%% ANALYTICAL FIELD SOLUTION FROM STATOR CURRENT (SCRIPT)

% Matteo Leandro
% matteo.leandro@ntnu.no
% VERSION: 14Jan2021

% This code finds the field solution in any region of a 2-D section of a
% slotless machine, in terms of flux density distribution and one-compoment
% vector potential. The field solution in the airgap and magnets region was
% tested for both inrunner and outrunner topology. The solution treats
% solely the case with unitary magnets permeability. Air-cored topologies
% can be studied. Singolarities are treated as well. The iron permeability
% is assumed infinity


%% Machine parameters
% clearvars
% clc
% close all

mu_0 = 4*pi*1e-7; % air permeability

% run(top)
%% useful indices for series harmonics definition
m_J = 5;                                                                  % the total number of harmonics will be 2*(m_J+1)
x = 0:1:m_J;
y = x+1;

% current harmonics order
m = p*(2*x+1); % existing single-phase harmonic series components (n=1,3,5,...)

% sigma_n = (sin(pi*m./m(end))./(pi*m./m(end))).^3;


% sigma_m = [sigma_h,sigma_k];
% sigma_m = sigma_n;
%% circumferential discretization

% sec = 2;                                                                    % number of poles to be modeled
% m_th = 1000*sec;                                                            % points along the modeled sector
% mechanical angle
% theta = linspace(0,sec*pi/(p),m_th);
% Theta = repmat(theta,2*(m_J+1),1);
% Theta = repmat(theta,(m_J+1),1);

%% Current density distribution
I_tot = N_tc*q*I/b;                                                         % total current through the phase belt
J_Ph = I_tot/S_Ph;                                                          % phase current density
J_1_Ph = 4*p./(m*pi).*J_Ph.*sin(pi*m/(6*p));


% J_theta = J_m*cos(m'.*Theta); % three-phase contribution
% J_theta = sigma_n.*J_1_Ph*cos(m'.*Theta); % single-phase contribution

% A_mJ = mu_0*J_m./(m.^2-4);
A_nJ = mu_0*J_1_Ph./(m.^2-4);

if p == 2
   A_nJ(1) = -mu_0*J_1_Ph(1)/4; 
end


% HARMONICS COEFFICIENTS FROM THE MAGNETIC VECTOR POTENTIAL
% coefficients which keep constant in the different regions

DEN_JI = 2*((R_i/R_s).^(2*m)-1);

DEN_JO = 2*((R_s/R_i).^(2*m)-1);

if p == 2
    DEN_JI(1) = 2*(R_s^4-R_i^4);
    DEN_JO(1) = DEN_JI(1);
    if R_i == Inf
       DEN_JO(1) = -2;
   end 
end

%% WINDING REGION 
% HARMONICS COEFFICIENTS FROM THE MAGNETIC VECTOR POTENTIAL
% (WINDING REGION)
if R_s>R_m %  inrunner

    A_pl_J_w = A_nJ.*(R_ws.*(2-(R_ws/R_s).^(2*m).*(m-2)+m)-R_w*(R_ws/R_s).^(2*m).*(R_w/R_ws).^(m+1).*(2-m+(R_i/R_w).^(2*m).*(2+m)))./DEN_JI;

    A_mi_J_w = A_nJ.*(R_ws*(R_w/R_ws).^(m-1).*(R_i/R_w).^(2*m).*(2-(R_ws/R_s).^(2*m).*(m-2)+m)-R_w*(2-m+(R_i/R_w).^(2*m).*(2+m)))./DEN_JI;

    if p == 2
       A_pl_J_w(1) = R_ws*A_nJ(1)*(R_w^4+R_i^4-R_s^4-R_ws^4-4*R_s^4*log(R_ws)+4*R_i^4*log(R_w))/DEN_JI(1);
       A_mi_J_w(1) = A_nJ(1)/R_w^3*(R_s^4*R_w^4-R_ws^4*R_i^4+4*R_s^4*R_i^4*log(R_w/R_ws))/DEN_JI(1);
    end

    i = 0;
    FLUX_const = (n_cs*l_a*2*p/S_Ph*(1./(m+2).*A_pl_J_w./m.*R_ws^3.*(1-(R_w/R_ws).^(m+2))+1./(2-m).*A_mi_J_w./m.*R_w^3.*((R_w/R_ws).^(m-2)-1)+A_nJ./4.*(R_ws^4-R_w^4)));
    if p == 2
        FLUX_const(1) = n_cs*l_a*2*p/S_Ph*(A_pl_J_w(1)/8.*R_ws^3.*(1-(R_w/R_ws)^(4))+A_mi_J_w(1)/2*R_w^3*ln(R_ws/R_w)+A_nJ/4*(R_ws^4-R_w^4));
    end
    
    FLUX_self = FLUX_const.*1./m.*(sin(m*(2*i+1)*pi./(6*p))-sin(m*(2*i-1)*pi./(6*p)));
    i = 1;
    FLUX_mut = FLUX_const.*1./m.*(sin(m*(2*i+1)*pi./(6*p))-sin(m*(2*i-1)*pi./(6*p)));
    L_self = sum(FLUX_self)/I;
    L_mut = sum(FLUX_mut)/I;
    L_sync = L_self + L_mut;
    pitch = pi* (2*R_w+2*R_ws)/2/p;
    Ls_ew = 4*pi*1e-7*pitch/2*1.428*n_cs^2*p/2;
    L_ph = (Ls_ew+L_sync);
    L_star = 2*(Ls_ew+L_sync);
    
else % outrunner
    
    A_pl_J_w = A_nJ.*(R_w.*(2-(R_w/R_i).^(2*m).*(m-2)+m)-R_ws*(R_w/R_i).^(2*m).*(R_ws/R_w).^(m+1).*(2-m+(R_s/R_ws).^(2*m).*(2+m)))./DEN_JO;

    A_mi_J_w = A_nJ.*(R_w*(R_ws/R_w).^(m-1).*(R_s/R_ws).^(2*m).*(2-(R_w/R_i).^(2*m).*(m-2)+m)+R_ws*(m-2-(R_s/R_ws).^(2*m).*(2+m)))./DEN_JO;

    if p == 2
       A_pl_J_w(1) = R_w*A_nJ(1)*(R_w^4+R_i^4-R_s^4-R_ws^4-4*R_s^4*log(R_ws)+4*R_i^4*log(R_w))/DEN_JO(1);
       A_mi_J_w(1) = A_nJ(1)/R_ws^3*(R_s^4*R_w^4-R_ws^4*R_i^4+4*R_s^4*R_i^4*log(R_w/R_ws))/DEN_JO(1);
       if R_i == Inf 
          A_pl_J_w(1) = R_w*A_nJ(1)*(1+4*log(R_w))/DEN_JO(1); 
          A_mi_J_w(1) = A_nJ(1)/R_ws^3*(4*R_s^4*log(R_w/R_ws)-R_ws^4)/DEN_JO(1);
       end
    end
    
    i = 0;
    FLUX_const = n_cs*l_a*2*p/S_Ph*(1./(m+2).*A_pl_J_w./m.*R_w^3.*(1-(R_ws/R_w).^(m+2))+1./(2-m).*A_mi_J_w./m.*R_ws^3.*((R_ws/R_w).^(m-2)-1)+A_nJ./4.*(R_w^4-R_ws^4));
    if p == 2
        FLUX_const(1) = n_cs*l_a*2*p/S_Ph*(A_pl_J_w(1)/8.*R_w^3.*(1-(R_ws/R_w)^(4))+A_mi_J_w(1)/2*R_ws^3*ln(R_w/R_ws)+A_nJ/4*(R_w^4-R_ws^4));
    end
    
    FLUX_self = FLUX_const.*1./m.*(sin(m*(2*i+1)*pi./(6*p))-sin(m*(2*i-1)*pi./(6*p)));
    i = 1;
    FLUX_mut = FLUX_const.*1./m.*(sin(m*(2*i+1)*pi./(6*p))-sin(m*(2*i-1)*pi./(6*p)));
    L_self = sum(FLUX_self)/I;
    L_mut = sum(FLUX_mut)/I;
    L_sync = L_self + L_mut;
    pitch = pi* (2*R_w+2*R_ws)/2/p;
    Ls_ew = 4*pi*1e-7*pitch/2*1.428*n_cs^2*p/2;
    L_ph = (Ls_ew+L_sync);
    L_star = 2*(Ls_ew+L_sync);
 
end
    

