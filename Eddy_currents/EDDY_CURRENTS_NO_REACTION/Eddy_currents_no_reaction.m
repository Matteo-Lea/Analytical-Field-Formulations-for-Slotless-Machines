%% ANALYTICAL EDDY CURRENTS SOLUTION FROM STATOR CURRENT (SCRIPT)

% Matteo Leandro
% matteo.leandro@ntnu.no
% VERSION: 28Sep2021

% This code uses the armature field solution accounting for inverter
% harmonics obtained from a dedicated code to calculate induced eddy
% current losses in the magnets region. The eddy current reaction field is
% neglected, but circumferential segmentation is considered. In the
% analytical expressions, the interaction between different harmonic orders
% is neglected. Inrunner and Outrunner topologies can be considered.
% However, singularities in the eddy-current density and power
% integrals need to be fixed.


%% Machine parameters
% clearvars
% clc
% close all



%% inrunner example
% top = 'Inrunner'; % Choose either 'Inrunner' or 'Outrunner'


% inductance is obtained from the armature field solution
Inverter_star_3f
clearvars -except n top pos neg phi_p phi_n T_fund runs

map = 'no'; % 'yes' if current density map is wanted 'no' otherwise!
plotting = 'yes'; % 'yes' if debugging is needed (plots and other stuff)

Inrunner
omega = rpm*2*pi/60; % mechanical angular frequency [rad/s]

if contains(plotting, 'yes')
    t = linspace(0,T_fund,length(n));
    figure;hold on
    plot(t*1000,(pos)'*cos(n.*p*omega*t+phi_p)+(neg)'*cos(n.*p*omega*t+phi_n))
    plot(t*1000,(pos)'*cos(n.*p*omega*t+phi_p-2/3*pi)+(neg)'*cos(n.*p*omega*t+phi_n+2/3*pi))
    plot(t*1000,(pos)'*cos(n.*p*omega*t+phi_p-4/3*pi)+(neg)'*cos(n.*p*omega*t+phi_n+4/3*pi))
    xlabel('Time [ms]')
    ylabel('Current [A]')
end

cond = 0.667*1e6;                                                           % PM conductivity [S/m]
mu_0 = 4*pi*1e-7; % air permeability


alpha_m = 0.5;
theta_mid = alpha_m*pi/p;
theta_side = (1-alpha_m)*pi/p;
%% outrunner example
% Outrunner
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

filt = [(mod(TL_m,3)==0), mod(TL_p,3)==0, mod(TL_p,3)==0, (mod(TL_m,3)==0)];

H = [pos_m, pos_p, neg_m, neg_p].*filt;
ch = [pos, pos, neg, neg].*filt;
phi = [-phi_p, phi_p, -phi_n, phi_n].*filt;

TL_m(mod(TL_m,3)~=0)=0;
TL_p(mod(TL_p,3)~=0)=0;

m = [m, m, m, m]*p;

n = [n, n, n, n];

%% Harmonic fliters for Gibbs phenomenon reduction
order_n = 0;
order_m = 0;
sigma_n = (sin(pi*n./n(end))./(pi*n./n(end))).^order_n; % Lanczos sigma for Gibbs phenomenon reduction
sigma_m = (sin(pi*m./m(end))./(pi*m./m(end))).^order_m; % Lanczos sigma for Gibbs phenomenon reduction

ch = sigma_n.*ch;

%% INDEPENDENT HARMONICS CONTRIBUTION

J_m = 4*3*p*N_tc*q*ch/b/S_Ph/(2*pi);                                        % constant coefficient for current density distribution
J_m = (J_m.*sin(pi*m/(6*p))./m);
A_mJ = mu_0*sigma_m.*J_m./(m.^2-4);
if p == 2
   A_mJ(:,1) = -mu_0*J_m(:,1)/4; 
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

% Phase-a field coefficient
    A_J = ((R_w/R_ws).^m.*(R_ws/R_s).^(2*m).*(R_w^2*(R_w/R_ws).^m-R_ws^2).*A_mJ.*(m-2)-((R_w^2-R_ws^2*(R_w/R_ws).^m).*A_mJ.*(m+2)))./(R_w*DEN_JI);   
    if p == 2
       A_J(1) = A_mJ(1)*R_w*(R_w^4-R_ws^4+4*R_s^4*log(R_w/R_ws))/DEN_JI(1);
    end

P1 = l_a*cond*(R_w^2*A_J./m.*H*p*omega).^2./2;
P1 = P1.*(1./(2*m+2).*((R_m/R_w).^(2*m+2)-(R_r/R_w).^(2*m+2))+1./(2-2*m).*(R_i/R_w).^(2*m+2).*((R_i/R_m).^(2*m-2)-(R_i/R_r).^(2*m-2))+(R_i/R_w).^(2*m).*((R_m/R_w)^2-(R_r/R_w)^2));
P1(H==0) = 0;
P1_mid = P1*theta_mid;
P1_side = P1*theta_side;

P2 = l_a*cond.*(A_J./m.^2*R_w^3).^2.*4./(R_m^2-R_r^2).*(H*p*omega).^2;
P2 = P2.*(1./(m+2).*((R_m/R_w).^(m+2)-(R_r/R_w).^(m+2))+1./(2-m).*(R_i/R_w).^(m+2).*((R_i/R_m).^(m-2)-(R_i/R_r).^(m-2))).^2;
P2(H==0) = 0;
P2_mid = -P2.*sin(m*theta_mid./2).^2./theta_mid;
P2_side = -P2.*sin(m*theta_side./2).^2./theta_side;

P_mid = 2*p*(P1_mid+P2_mid);
P_mid_main = sum(sum(P_mid)); % Contribution from each harmonic

P_side = 2*p.*(P1_side+P2_side);
P_side_main = sum(sum(P_side)); % Contribution from each harmonic


%% The 2 different cases
%% 1st case: (H)i=(H)j


% THE FOLLOWING PART IS COMPUTATIONALLY MORE EFFICIENT BUT THE CASES OF
% Hi=Hj WHERE SPACE AND TIME HARMONICS ARE THE SAME SHOULD BE TREATED
% SEPARATELY AS THEY ARE NOT INCLUDED. HOWEVWER, THE RELATED EXPRESSIONS ARE
% MUCH SIMPLER.
%% CASE 3Li=3Lj
[S, idx] = sort(H(:).'); % the harmonic matrix TL is sorted over a row vector
ind = [false, diff(S) == 0]; % logic vector placing a zero whenever two     
                             % consecutive elements give zero difference, 
                             % i.e., they are equal (a zero is added at the
                             % beginning for the first element)
% zero harmonics are disregarded for they do not contribute to the losses
idx(S==0) = [];
ind(S==0) = [];

start = strfind(ind, [false, true]);
ending = [start(2)-1 start(3:end)-1 length(idx)];

seq = cell(length(start),1);
for k = 1:length(start)
    seq{k} = nchoosek(idx(start(k):ending(k)),2);
end
%--------------------------------------------------------------------------

seq = cell2mat(seq);

m_i = m(seq(:,1));
m_j = m(seq(:,2));
phi_ni = phi(seq(:,1));
phi_nj = phi(seq(:,2));
tl_mut = H(seq(:,1));
ch_i = ch(seq(:,1));
ch_j = ch(seq(:,2));
sigma_mi = (sin(pi*m_i./m(end))./(pi*m_i./m(end))).^order_m; % Lanczos sigma for Gibbs phenomenon reduction
sigma_mj = (sin(pi*m_j./m(end))./(pi*m_j./m(end))).^order_m; % Lanczos sigma for Gibbs phenomenon reduction

%% Current density distribution
% PHASE A
J_mi = 4*3*p*N_tc*q*ch_i/b/S_Ph/(2*pi);                                    % constant coefficient for current density distribution
J_mj = 4*3*p*N_tc*q*ch_j/b/S_Ph/(2*pi);
J_mi = (J_mi.*sin(pi*m_i/(6*p))./m_i);
J_mj = (J_mj.*sin(pi*m_j/(6*p))./m_j);
A_mJi = mu_0*sigma_mi.*J_mi./(m_i.^2-4);
A_mJj = mu_0*sigma_mj.*J_mj./(m_j.^2-4);
if p == 2
   A_mJ(:,1) = -mu_0*J_m(:,1)/4; 
end

alpha_m = 0.5;
theta_mid = alpha_m*pi/p;
theta_side = (1-alpha_m)*pi/p;
mid_pos = 0;
side_pos = pi/(2*p);

    
%% AIRGAP/MAGNETS/AIR-BACKING REGION    
% HARMONICS COEFFICIENTS FROM THE MAGNETIC VECTOR POTENTIAL
% (AIR-GAP/MAGNETS/(NON-MAGNETIC-BACKING) REGIONS)
% ONLY INRUNNER IMPLEMENTED 
    
    % HARMONICS COEFFICIENTS FROM THE MAGNETIC VECTOR POTENTIAL
    % coefficients which keep constant in the different regions
    
    DEN_JIi = 2*((R_i/R_s).^(2*m_i)-1);
    DEN_JIj = 2*((R_i/R_s).^(2*m_j)-1);

    if p == 2
        DEN_JIi(1) = 2*(R_s^4-R_i^4);
    end
    
    % Phase-a field coefficient
    A_Ji = ((R_w/R_ws).^m_i.*(R_ws/R_s).^(2*m_i).*(R_w^2*(R_w/R_ws).^m_i-R_ws^2).*A_mJi.*(m_i-2)-((R_w^2-R_ws^2*(R_w/R_ws).^m_i).*A_mJi.*(m_i+2)))./(R_w*DEN_JIi); 
    A_Jj = ((R_w/R_ws).^m_j.*(R_ws/R_s).^(2*m_j).*(R_w^2*(R_w/R_ws).^m_j-R_ws^2).*A_mJj.*(m_j-2)-((R_w^2-R_ws^2*(R_w/R_ws).^m_j).*A_mJj.*(m_j+2)))./(R_w*DEN_JIj);
    if p == 2
       A_J(1) = A_mJ(1)*R_w*(R_w^4-R_ws^4+4*R_s^4*log(R_w/R_ws))/DEN_JIi(1);
    end
 
    R0 = 1./(m_j+m_i+2).*((R_m/R_w).^(m_j+m_i)*R_m^2-(R_r/R_w).^(m_j+m_i)*R_r^2);
    R0 = R0 + 1./(2-m_j+m_i).*(R_i/R_w).^(m_j).*((R_i/R_m).^(m_j).*(R_m/R_w).^(m_i)*R_m^2-(R_i/R_r).^(m_j).*(R_r/R_w).^(m_i)*R_r^2);
    R0 = R0 + 1./(2-m_i+m_j).*(R_i/R_w).^(m_i).*((R_i/R_m).^(m_i).*(R_m/R_w).^(m_j)*R_m^2-(R_i/R_r).^(m_i).*(R_r/R_w).^(m_j)*R_r^2);
    R0 = R0 + 1./(2-m_j-m_i).*(R_i/R_w).^(m_i+m_j).*((R_i/R_m).^(m_j+m_i)*R_m^2-(R_i/R_r).^(m_j+m_i)*R_r^2);
    
    P1 = R_w^2*l_a*cond*(p*omega)^2.*R0./(m_i.*m_j).*A_Ji.*A_Jj.*tl_mut.^2;
    P1_mid = P1./(m_i-m_j).*sin((m_i-m_j)*theta_mid./2).*cos((m_i-m_j)*mid_pos+phi_ni-phi_nj);
    P1_mid((m_i-m_j)==0,:) = P1((m_i-m_j)==0,:)./2*theta_mid.*cos(phi_ni((m_i-m_j)==0,:)-phi_nj((m_i-m_j)==0,:));
    P1_side = P1./(m_i-m_j).*sin((m_i-m_j)*theta_side./2).*cos((m_i-m_j)*side_pos+phi_ni-phi_nj);
    P1_side((m_i-m_j)==0,:) = P1((m_i-m_j)==0,:)./2*theta_side.*cos(phi_ni((m_i-m_j)==0,:)-phi_nj((m_i-m_j)==0,:));
    
    
    ej = 1./(m_j+2).*((R_m/R_w).^(m_j+2)-(R_r/R_w).^(m_j+2));
    ej = ej + 1./(2-m_j).*(R_i/R_w).^(m_j+2).*((R_i/R_m).^(m_j-2)-(R_i/R_r).^(m_j-2));
    ei = 1./(m_i+2).*((R_m/R_w).^(m_i+2)-(R_r/R_w).^(m_i+2));
    ei = ei + 1./(2-m_i).*(R_i/R_w).^(m_i+2).*((R_i/R_m).^(m_i-2)-(R_i/R_r).^(m_i-2));
    P2 = -l_a*cond*(R_w^3*tl_mut*p*omega).^2.*A_Ji.*A_Jj./(m_i.*m_j)*4./((R_m^2-R_r^2).*(m_i).*(m_j)).*ej.*ei;
    P2_mid = P2.*sin(m_j*theta_mid./2).*sin(m_i*theta_mid./2).*cos(phi_ni-phi_nj+(m_i-m_j)*mid_pos)./theta_mid;
    P2_side = P2.*sin(m_j*theta_side./2).*sin(m_i*theta_side./2).*cos(phi_ni-phi_nj+(m_i-m_j)*p*side_pos)./theta_side;
    
    P_mid = 2*p.*(P1_mid+P2_mid);
    
    P_side = 2*p.*(P1_side+P2_side);
    %% 2nd case: (H)i=-(H)j
    
    idx = find(H<0); % isolate indices of negative harmonics 

    % find in TL all the combinations of harmonics with opposite sign
    seq = cell(length(idx),1);
    for k = 1:length(idx)
        seq{k} = find(H==-H(idx(k)));
        seq{k} = [seq{k} repmat(idx(k),length(seq{k}),1)];
    end

    seq = cell2mat(seq);

    m_i = m(seq(:,1));
    m_j = m(seq(:,2));
    phi_ni = phi(seq(:,1));
    phi_nj = phi(seq(:,2));
    tl_mut = H(seq(:,1));
    ch_i = ch(seq(:,1));
    ch_j = ch(seq(:,2));
    sigma_mi = (sin(pi*m_i./m(end))./(pi*m_i./m(end))).^order_m; % Lanczos sigma for Gibbs phenomenon reduction
    sigma_mj = (sin(pi*m_j./m(end))./(pi*m_j./m(end))).^order_m; % Lanczos sigma for Gibbs phenomenon reduction

    %% Current density distribution
    % PHASE A
    J_mi = 4*3*p*N_tc*q*ch_i/b/S_Ph/(2*pi);                                        % constant coefficient for current density distribution
    J_mj = 4*3*p*N_tc*q*ch_j/b/S_Ph/(2*pi);
    J_mi = (J_mi.*sin(pi*m_i/(6*p))./m_i);
    J_mj = (J_mj.*sin(pi*m_j/(6*p))./m_j);
    A_mJi = mu_0*sigma_mi.*J_mi./(m_i.^2-4);
    A_mJj = mu_0*sigma_mj.*J_mj./(m_j.^2-4);
    if p == 2
       A_mJ(:,1) = -mu_0*J_m(:,1)/4; 
    end
    
    % HARMONICS COEFFICIENTS FROM THE MAGNETIC VECTOR POTENTIAL
    % coefficients which keep constant in the different regions
    
    DEN_JIi = 2*((R_i/R_s).^(2*m_i)-1);
    DEN_JIj = 2*((R_i/R_s).^(2*m_j)-1);

    if p == 2
        DEN_JIi(1) = 2*(R_s^4-R_i^4);
    end
    
    % Phase-a field coefficient
    A_Ji = ((R_w/R_ws).^m_i.*(R_ws/R_s).^(2*m_i).*(R_w^2*(R_w/R_ws).^m_i-R_ws^2).*A_mJi.*(m_i-2)-((R_w^2-R_ws^2*(R_w/R_ws).^m_i).*A_mJi.*(m_i+2)))./(R_w*DEN_JIi); 
    A_Jj = ((R_w/R_ws).^m_j.*(R_ws/R_s).^(2*m_j).*(R_w^2*(R_w/R_ws).^m_j-R_ws^2).*A_mJj.*(m_j-2)-((R_w^2-R_ws^2*(R_w/R_ws).^m_j).*A_mJj.*(m_j+2)))./(R_w*DEN_JIj);
    if p == 2
       A_J(1) = A_mJ(1)*R_w*(R_w^4-R_ws^4+4*R_s^4*log(R_w/R_ws))/DEN_JIi(1);
    end
    
    R0 = 1./(m_j+m_i+2).*((R_m/R_w).^(m_j+m_i)*R_m^2-(R_r/R_w).^(m_j+m_i)*R_r^2);
    R0 = R0 + 1./(2-m_j+m_i).*(R_i/R_w).^(m_j).*((R_i/R_m).^(m_j).*(R_m/R_w).^(m_i)*R_m^2-(R_i/R_r).^(m_j).*(R_r/R_w).^(m_i)*R_r^2);
    R0 = R0 + 1./(2-m_i+m_j).*(R_i/R_w).^(m_i).*((R_i/R_m).^(m_i).*(R_m/R_w).^(m_j)*R_m^2-(R_i/R_r).^(m_i).*(R_r/R_w).^(m_j)*R_r^2);
    R0 = R0 + 1./(2-m_j-m_i).*(R_i/R_w).^(m_i+m_j).*((R_i/R_m).^(m_j+m_i)*R_m^2-(R_i/R_r).^(m_j+m_i)*R_r^2);
    
    P1 = R_w^2*l_a*cond*(p*omega)^2.*R0./(m_i.*m_j).*A_Ji.*A_Jj.*tl_mut.^2;
    P1_mid2 = P1./(m_i+m_j).*sin((m_i+m_j)*theta_mid./2).*cos((m_i+m_j)*mid_pos+phi_ni+phi_nj);
    P1_side2 = P1./(m_i+m_j).*sin((m_i+m_j)*theta_side./2).*cos((m_i+m_j)*side_pos+phi_ni+phi_nj);
    
    
    ej = 1./(m_j+2).*((R_m/R_w).^(m_j+2)-(R_r/R_w).^(m_j+2));
    ej = ej + 1./(2-m_j).*(R_i/R_m).^(m_j+2).*((R_i/R_m).^(m_j-2)-(R_i/R_r).^(m_j-2));
    ei = 1./(m_i+2).*((R_m/R_w).^(m_i+2)-(R_r/R_w).^(m_i+2));
    ei = ei + 1./(2-m_i).*(R_i/R_m).^(m_i+2).*((R_i/R_m).^(m_i-2)-(R_i/R_r).^(m_i-2));
    P2 = -l_a*cond*(R_w^3*tl_mut*p*omega).^2.*A_Ji.*A_Jj./(m_i.*m_j)*4./((R_m^2-R_r^2).*(m_i).*(m_j)).*ej.*ei;
    P2_mid2 = P2.*sin(m_j*theta_mid./2).*sin(m_i*theta_mid./2).*cos(phi_ni+phi_nj+(m_i+m_j)*mid_pos)./theta_mid;
    P2_side2 = P2.*sin(m_j*theta_side./2).*sin(m_i*theta_side./2).*cos(phi_ni+phi_nj+(m_i+m_j)*p*side_pos)./theta_side;
    
%% Indipendent contribution of existing harmonics
P1 = sum(sum(2*p*[P1_mid; P1_mid2]));
P2 = sum(sum(2*p*[P2_mid; P2_mid2]));
P_tot_mid = P1 + P2;

P1 = sum(sum(2*p*[P1_side; P1_side2]));
P2 = sum(sum(2*p*[P2_side; P2_side2]));
P_tot_side = P1 + P2;

TOT_mid = P_mid_main+P_tot_mid;

TOT_side = P_side_main+P_tot_side;

TOT = P_mid_main+P_tot_mid + P_side_main+P_tot_side;


%% SOME DEBUGGING
if contains(plotting, 'yes')
r = (R_m+R_r)/2; % test radius 
Amp = R_w.*A_J./m.*((r/R_w).^m+(R_i/R_w).^m.*(R_i/r).^m);
    
t = repmat(shiftdim(linspace (0,T_fund,size(n,1)+1),-1),size(n,1),size(m,2));
theta = 0;
J0 = p*omega*cond*Amp.*H.*sin(H*p*omega.*t+phi+m.*(theta)) ;
J0 = squeeze(sum(sum(J0,2),1));
figure
plot(linspace (0,T_fund,size(n,1)+1), J0)


Az = Amp.*(cos(H*p*omega.*t+phi+m.*(theta)));
    test = squeeze(sum(sum(Az,2),1));
    figure
    plot(linspace (0,T_fund,size(n,1)+1),test)

end

%% MAPPPING
if contains(map, 'yes')

    r_mid = shiftdim(linspace(R_r,R_m,10),-1); % test radius 
    theta_m = repmat(shiftdim(linspace(-alpha_p*pi/(2*p),alpha_p*pi/(2*p),100),-2),1,1,size(r_mid,3),1);
    r_mid = repmat(r_mid,1,1,1,size(theta_m,4));
    Amp = R_w.*A_J./m.*((r_mid/R_w).^m+(R_i/R_w).^m.*(R_i/r_mid).^m);
    
    t = 1.822277847309136e-04;
    
    J0 = p*omega*cond*Amp.*H.*sin(H*p*omega.*t+phi+m.*(theta_m)) ;
    J0 = squeeze(sum(sum(J0,2),1));


    Ja = -cond*(4*R_w^3*A_J./(theta_mid*(R_m^2-R_r^2)*m)).*(1./(m+2).*((R_m/R_w).^(m+2)-(R_r/R_w).^(m+2))+1./(2-m).*(R_i/R_w).^(m+2).*((R_i/R_m).^(m-2)-(R_i/R_r).^(m-2))).*(H*p*omega./(m)).*sin(m*theta_mid./2).*sin(H*p*omega*t+phi+m*mid_pos);
    J_mid = J0 + sum(sum(Ja));
    
    r_side = shiftdim(linspace(R_r,R_m,10),-1); % test radius 
    theta_s = repmat(shiftdim(linspace(alpha_p*pi/(2*p),(1-alpha_p/2)*pi/(p),100),-2),1,1,size(r_side,3),1);
    r_side = repmat(r_side,1,1,1,size(theta_s,4));
    Amp = R_w.*A_J./m.*((r_side/R_w).^m+(R_i/R_w).^m.*(R_i/r_side).^m);
    
    J0 = p*omega*cond*Amp.*H.*sin(H*p*omega.*t+phi+m.*(theta_s)) ;
    J0 = squeeze(sum(sum(J0,2),1));


    Ja = -cond*(4*R_w^3*A_J./(theta_side*(R_m^2-R_r^2)*m)).*(1./(m+2).*((R_m/R_w).^(m+2)-(R_r/R_w).^(m+2))+1./(2-m).*(R_i/R_w).^(m+2).*((R_i/R_m).^(m-2)-(R_i/R_r).^(m-2))).*(H*p*omega./(m)).*sin(m*theta_side./2).*sin(H*p*omega*t+phi+m*side_pos);
    J_side = J0 + sum(sum(Ja));
    
    [r_m,t_m] = meshgrid(linspace(R_r,R_m,10),-linspace(-alpha_p*pi/(p),0,100));
    [r_s,t_s] = meshgrid(linspace(R_r,R_m,10),-linspace(0,alpha_p*pi/(p),100));
    x_m = r_m.*sin(t_m);
    y_m = r_m.*cos(t_m);
    x_s = r_s.*sin(t_s);
    y_s = r_s.*cos(t_s);
    
    theta = -linspace(-alpha_p*pi/(p),(alpha_p)*pi/(p),10); % circumferential discretization (leave the -pi/(2*p))
    
    J_max = max(abs([min(min(J_mid)) min(min(J_side)) max(max(J_mid)) max(max(J_side))]));
    J_levels = 100;
    th1 = -(alpha_m*pi/p);
    th2 = +(alpha_m*pi/p);
    %---------------------------------------------------------------------
    
    figure;
    hold on;
    contourf(x_m,y_m,J_mid',J_levels,'edgecolor','none');
    contourf(x_s,y_s,J_side',J_levels,'edgecolor','none');
    plot(R_r*sin(linspace(th1,th2,100)),R_r*cos(linspace(th1,th2,100)),'linewidth',0.8,'color','k');
    plot(R_m*sin(linspace(th1,th2,100)),R_m*cos(linspace(th1,th2,100)),'linewidth',0.8,'color','k');
    plot(R_s*sin(linspace(th1,th2,100)),R_s*cos(linspace(th1,th2,100)),'linewidth',0.8,'color','k');
    plot(R_se*sin(linspace(th1,th2,100)),R_se*cos(linspace(th1,th2,100)),'linewidth',0.8,'color','k');
    plot([R_r*sin(th1);R_m*sin(th1)],[R_r*cos(th1);R_m*cos(th1)],'linewidth',0.8,'color','k');
    plot([R_r*sin(th2);R_m*sin(th2)],[R_r*cos(th2);R_m*cos(th2)],'linewidth',0.8,'color','k');
    plot([R_r*sin(0);R_m*sin(0)],[R_r*cos(0);R_m*cos(0)],'linewidth',0.8,'color','k');
    % Winding slots
    plot(R_w*sin(linspace(th1,th2,100)),R_w*cos(linspace(th1,th2,100)),'linewidth',0.8,'color','k');
    plot(R_ws*sin(linspace(th1,th2,100)),R_ws*cos(linspace(th1,th2,100)),'linewidth',0.8,'color','k');
    plot([R_ws*sin(th1+pi/(6*p):pi/(3*p):th2-pi/(6*p));R_w*sin(th1+pi/(6*p):pi/(3*p):th2-pi/(6*p))],[R_ws*cos(th1+pi/(6*p):pi/(3*p):th2-pi/(6*p));R_w*cos(th1+pi/(6*p):pi/(3*p):th2-pi/(6*p))],'linewidth',0.8,'color','k');
    set(gca,'visible','off');
    colormap(jet)
    c = colorbar;
    c.Label.String = 'Flux density norm [T]';
    caxis([-2e6 2e6])
    axis off
    axis image
    % camroll(-90)
    title({'Induced current density map',' in the magnets region'})
    
end

