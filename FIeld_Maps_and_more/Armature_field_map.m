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
clearvars
clc
close all

mu_0 = 4*pi*1e-7; % air permeability

%% inrunner example
Inrunner


%% outrunner example
% Outrunner

mapping = "yes";
%% useful indices for series harmonics definition
m_J = 4;                                                                  % the total number of harmonics will be 2*(m_J+1)
x = 0:1:m_J;
y = x+1;
% current harmonics order
h = p*(6*x+1);
k = p*(6*y-1);
m = [h,k];                                                                  % the harmponics are put together in a single vector
 
sigma_h = (sin(pi*h./h(end))./(pi*h./h(end))).^3;                           % Lanczos sigma for Gibbs phenomenon reduction
sigma_k = (sin(pi*k./k(end))./(pi*k./k(end))).^3;                           % Lanczos sigma for Gibbs phenomenon reduction

sigma_m = [sigma_h,sigma_k];
%% circumferential discretization

sec = 1;                                                                    % number of poles to be modeled
m_th = 1000*sec;                                                            % points along the modeled sector
% mechanical angle
theta = linspace(0,sec*pi/(p),m_th);
Theta = repmat(theta,2*(m_J+1),1);

%% Current density distribution
I_tot = N_tc*q*I/b;                                                         % total current through the phase belt
J_Ph = I_tot/S_Ph;                                                          % phase current density
J_coeff = 4*3*p*J_Ph/(2*pi);                                                % constant coefficient for current density distribution

J_h = (J_coeff*sin(pi*h/(6*p))./(h));
J_k = (J_coeff*sin(pi*k/(6*p))./(k));

J_m = [J_h,J_k];

J_theta = J_m*cos(m'.*Theta);

A_mJ = mu_0*J_m./(m.^2-4);

if p == 2
   A_mJ(1) = -mu_0*J_m(1)/4; 
end


%% RADIAL DISCRETIZATION 
%iron backing
if R_i == R_r 
    % inrunner
    if R_s >R_m 
        R_ie = 0.9*R_r;
        r_ext = linspace(R_ie,R_r,50)';
        r_g = linspace(R_r,R_w,50)';                                        % airgap/magnets/air-backing radial discretization
        % outrunner
    else 
        R_ie = 1.1*R_r; 
        r_ext = linspace(R_r,R_ie,50)';
        r_g = linspace(R_w,R_r,50)';                                        % airgap/magnets/air-backing radial discretization
    end
% no iron backing

else % no iron backing
    if R_s >R_m % inrunner
        R_ie = R_r-0.25*pi*R_r/p;
        r_ext = linspace(R_ie,R_r,50)';
        r_g = linspace(R_ie,R_w,50)';                                       % airgap/magnets/air-backing radial discretization
    else % outrunner
        R_ie = R_r+0.5*pi*R_r/p;
        r_ext = linspace(R_r,R_ie,50)';
        r_g = linspace(R_w,R_ie,50)';                                       % airgap/magnets/air-backing radial discretization
    end
end

% inrunner
if R_s>R_m 
r_w = linspace(R_w,R_ws,50)';                                               % winding radial discretization
r_bw = linspace(R_ws,R_s,50)';                                              % air back-winding region discretization
r_s = linspace(R_s,R_se,50)';                                               % stator iron discretization
% outrunner
elseif R_s<R_m 
r_w = linspace(R_ws,R_w,50)';                                               % winding radial discretization
r_bw = linspace(R_s,R_ws,50)';                                              % air back-winding region discretization
r_s = linspace(R_se,R_s,50)';                                               % stator iron discretization
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

%% WINDING-TO-BACKIRON REGION REGION 
% HARMONICS COEFFICIENTS FROM THE MAGNETIC VECTOR POTENTIAL
% (WINDING-TO-BACKIRON REGION)
if R_s>R_m %  inrunner

    A_J_m_bw = ((R_w/R_ws).^m.*(R_ws^2*(R_w/R_ws).^m-R_w^2).*(R_i/R_w).^(2*m).*A_mJ.*(2+m)-((R_ws^2-R_w^2*(R_w/R_ws).^m).*A_mJ.*(m-2)))./(DEN_JI*R_ws); 

    if p == 2
       A_J_m_bw(1) = (R_s/R_ws)^3*A_mJ(1)*R_s*(R_w^4-R_ws^4+4*R_i^4*log(R_w/R_ws))/DEN_JI(1);
    end

    % radial component
    Amp_r_J_bw = (A_J_m_bw.*((r_bw./R_s).^(m-1).*(R_ws/R_s).^(m+1)+(R_ws./r_bw).^(m+1)));
    B_r_J_bw = -sigma_m.*Amp_r_J_bw*sin(m'.*Theta);

    % one-component magnetic potential 
    Amp_Az_J_bw = (A_J_m_bw./(m)*R_ws.*((r_bw./R_s).^(m).*(R_ws/R_s).^(m)+(R_ws./r_bw).^(m)));
    A_z_J_bw = sigma_m.*Amp_Az_J_bw*cos(m'.*Theta);

    % circumferential component
    Amp_theta_J_bw = (A_J_m_bw.*(-(r_bw./R_s).^(m-1).*(R_ws/R_s).^(m+1)+(R_ws./r_bw).^(m+1)));
    B_theta_J_bw = sigma_m.*Amp_theta_J_bw*cos(m'.*Theta);


else % outrunner
    
    A_J_m_bw = A_mJ.*((R_ws/R_w).^m.*(R_ws^2*(R_ws/R_w).^m-R_w^2).*(R_w/R_i).^(2*m).*(m-2)-((R_ws^2-R_w^2*(R_ws/R_w).^m).*(m+2)))./(DEN_JO*R_ws); 

    if p == 2
       A_J_m_bw(1) = A_mJ(1)*R_ws*(R_w^4-R_ws^4+4*R_i^4*log(R_w/R_ws))/DEN_JO(1);
       if R_i == Inf
           A_J_m_bw(1) = A_mJ(1)*R_ws*4*log(R_w/R_ws)/DEN_JO(1);
       end
    end

    % radial component
    Amp_r_J_bw = (A_J_m_bw.*((r_bw./R_ws).^(m-1)+(R_s./r_bw).^(m+1).*(R_s/R_ws).^(m-1)));
    B_r_J_bw = -sigma_m.*Amp_r_J_bw*sin(m'.*Theta);

    % one-component magnetic potential 
    Amp_Az_J_bw = (A_J_m_bw./(m)*R_ws.*((r_bw./R_ws).^(m)+(R_s./r_bw).^(m).*(R_s/R_ws).^(m)));
    A_z_J_bw = sigma_m.*Amp_Az_J_bw*cos(m'.*Theta);

    % circumferential component
    Amp_theta_J_bw = (A_J_m_bw.*(-(r_bw./R_ws).^(m-1)+(R_s./r_bw).^(m+1).*(R_s/R_ws).^(m-1)));
    B_theta_J_bw = sigma_m.*Amp_theta_J_bw*cos(m'.*Theta);
       
end

%% WINDING REGION 
% HARMONICS COEFFICIENTS FROM THE MAGNETIC VECTOR POTENTIAL
% (WINDING REGION)
if R_s>R_m %  inrunner

    A_pl_J_w = A_mJ.*(R_ws.*(2-(R_ws/R_s).^(2*m).*(m-2)+m)-R_w*(R_ws/R_s).^(2*m).*(R_w/R_ws).^(m+1).*(2-m+(R_i/R_w).^(2*m).*(2+m)))./DEN_JI;

    A_mi_J_w = A_mJ.*(R_ws*(R_w/R_ws).^(m-1).*(R_i/R_w).^(2*m).*(2-(R_ws/R_s).^(2*m).*(m-2)+m)-R_w*(2-m+(R_i/R_w).^(2*m).*(2+m)))./DEN_JI;

    if p == 2
       A_pl_J_w(1) = R_ws*A_mJ(1)*(R_w^4+R_i^4-R_s^4-R_ws^4-4*R_s^4*log(R_ws)+4*R_i^4*log(R_w))/DEN_JI(1);
       A_mi_J_w(1) = A_mJ(1)/R_w^3*(R_s^4*R_w^4-R_ws^4*R_i^4+4*R_s^4*R_i^4*log(R_w/R_ws))/DEN_JI(1);
    end

    % radial component
    Amp_r_J_w = (A_pl_J_w.*(r_w./R_ws).^(m-1)+A_mi_J_w.*(R_w./r_w).^(m+1)+m.*A_mJ.*r_w);
    if p == 2
        Amp_r_J_w(:,1) = (A_pl_J_w(1).*(r_w./R_ws)+A_mi_J_w(1).*(R_w./r_w).^(3)+2.*A_mJ(1).*r_w.*log(r_w));
    end
    B_r_J_w = -sigma_m.*Amp_r_J_w*sin(m'.*Theta);

    % one-component magnetic potential 
    Amp_Az_J_w = (A_pl_J_w./m*R_ws.*(r_w./R_ws).^(m)+A_mi_J_w./m*R_w.*(R_w./r_w).^(m)+A_mJ.*r_w.^2);
    if p == 2
        Amp_Az_J_w(:,1) = (A_pl_J_w(1)./2*R_ws.*(r_w./R_ws).^(2)+A_mi_J_w(1)./2*R_w.*(R_w./r_w).^(2)+A_mJ(1).*r_w.^2.*log(r_w));
    end
    A_z_J_w = sigma_m.*Amp_Az_J_w*cos(m'.*Theta);

    % circumferential component
    Amp_theta_J_w = (-A_pl_J_w.*(r_w./R_ws).^(m-1)+A_mi_J_w.*(R_w./r_w).^(m+1)-2*A_mJ.*r_w);
    if p == 2
        Amp_theta_J_w(:,1) = (-A_pl_J_w(1).*(r_w./R_ws)+A_mi_J_w(1).*(R_w./r_w).^(3)-2.*A_mJ(1).*r_w.*log(r_w)-A_mJ(1).*r_w);
    end
    B_theta_J_w = sigma_m.*Amp_theta_J_w*cos(m'.*Theta);

else % outrunner
    
    A_pl_J_w = A_mJ.*(R_w.*(2-(R_w/R_i).^(2*m).*(m-2)+m)-R_ws*(R_w/R_i).^(2*m).*(R_ws/R_w).^(m+1).*(2-m+(R_s/R_ws).^(2*m).*(2+m)))./DEN_JO;

    A_mi_J_w = A_mJ.*(R_w*(R_ws/R_w).^(m-1).*(R_s/R_ws).^(2*m).*(2-(R_w/R_i).^(2*m).*(m-2)+m)+R_ws*(m-2-(R_s/R_ws).^(2*m).*(2+m)))./DEN_JO;

    if p == 2
       A_pl_J_w(1) = R_w*A_mJ(1)*(R_w^4+R_i^4-R_s^4-R_ws^4-4*R_s^4*log(R_ws)+4*R_i^4*log(R_w))/DEN_JO(1);
       A_mi_J_w(1) = A_mJ(1)/R_ws^3*(R_s^4*R_w^4-R_ws^4*R_i^4+4*R_s^4*R_i^4*log(R_w/R_ws))/DEN_JO(1);
       if R_i == Inf 
          A_pl_J_w(1) = R_w*A_mJ(1)*(1+4*log(R_w))/DEN_JO(1); 
          A_mi_J_w(1) = A_mJ(1)/R_ws^3*(4*R_s^4*log(R_w/R_ws)-R_ws^4)/DEN_JO(1);
       end
    end

    % radial component
    Amp_r_J_w = (A_pl_J_w.*(r_w./R_w).^(m-1)+A_mi_J_w.*(R_ws./r_w).^(m+1)+m.*A_mJ.*r_w);
    if p == 2
        Amp_r_J_w(:,1) = (A_pl_J_w(1).*(r_w./R_w)+A_mi_J_w(1).*(R_ws./r_w).^(3)+2.*A_mJ(1).*r_w.*log(r_w));
    end
    B_r_J_w = -sigma_m.*Amp_r_J_w*sin(m'.*Theta);

    % one-component magnetic potential 
    Amp_Az_J_w = (A_pl_J_w./m*R_w.*(r_w./R_w).^(m)+A_mi_J_w./m*R_ws.*(R_ws./r_w).^(m)+A_mJ.*r_w.^2);
    if p == 2
        Amp_Az_J_w(:,1) = (A_pl_J_w(1)./2*R_w.*(r_w./R_w).^(2)+A_mi_J_w(1)./2*R_ws.*(R_ws./r_w).^(2)+A_mJ(1).*r_w.^2.*log(r_w));
    end
    A_z_J_w = sigma_m.*Amp_Az_J_w*cos(m'.*Theta);

    % circumferential component
    Amp_theta_J_w = (-A_pl_J_w.*(r_w./R_w).^(m-1)+A_mi_J_w.*(R_ws./r_w).^(m+1)-2*A_mJ.*r_w);
    if p == 2
        Amp_theta_J_w(:,1) = (-A_pl_J_w(1).*(r_w./R_w)+A_mi_J_w(1).*(R_ws./r_w).^(3)-2.*A_mJ(1).*r_w.*log(r_w)-A_mJ(1).*r_w);
    end
    B_theta_J_w = sigma_m.*Amp_theta_J_w*cos(m'.*Theta);
 
end
    
%% AIRGAP/MAGNETS/AIR-BACKING REGION    
% HARMONICS COEFFICIENTS FROM THE MAGNETIC VECTOR POTENTIAL
% (AIR-GAP/MAGNETS/(NON-MAGNETIC-BACKING) REGIONS)
if R_s>R_m %  inrunner
    
    A_J_m_g = ((R_w/R_ws).^m.*(R_ws/R_s).^(2*m).*(R_w^2*(R_w/R_ws).^m-R_ws^2).*A_mJ.*(m-2)-((R_w^2-R_ws^2*(R_w/R_ws).^m).*A_mJ.*(m+2)))./(R_w*DEN_JI);   
    
    if p == 2
       A_J_m_g(1) = A_mJ(1)*R_w*(R_w^4-R_ws^4+4*R_s^4*log(R_w/R_ws))/DEN_JI(1);
    end

    % radial component
    Amp_r_J_g = (A_J_m_g.*((r_g./R_w).^(m-1)+(R_i/R_w).^(m-1).*(R_i./r_g).^(m+1)));
    B_r_J_g = -sigma_m.*Amp_r_J_g*sin(m'.*Theta);

    % one-component magnetic potential 
    Amp_Az_J_g = (A_J_m_g./m*R_w.*((r_g./R_w).^(m)+(R_i/R_w).^(m).*(R_i./r_g).^(m)));
    A_z_J_g = sigma_m.*Amp_Az_J_g*cos(m'.*Theta);

    % circumferential component
    Amp_theta_J_g = (A_J_m_g.*(-(r_g./R_w).^(m-1)+(R_i/R_w).^(m-1).*(R_i./r_g).^(m+1)));
    B_theta_J_g = sigma_m.*Amp_theta_J_g*cos(m'.*Theta);

else % outrunner
    
    A_J_m_g = A_mJ.*((R_ws/R_w).^m.*(R_s/R_ws).^(2*m).*(R_w^2*(R_ws/R_w).^m-R_ws^2).*(m+2)+((-R_w^2+R_ws^2*(R_ws/R_w).^m).*(m-2)))./(R_w*DEN_JO);   

    if p == 2
       A_J_m_g(1) = A_mJ(1)*R_i*(R_i/R_w)^3*(R_w^4-R_ws^4+4*R_s^4*log(R_w/R_ws))/DEN_JO(1);
       if R_i == Inf
          A_J_m_g(1) = A_mJ(1)/R_w^3*(R_w^4-R_ws^4+4*R_s^4*log(R_w/R_ws))/DEN_JO(1);
       end
    end

    % radial component
    Amp_r_J_g = (A_J_m_g.*((r_g./R_i).^(m-1).*(R_w/R_i).^(m+1)+(R_w./r_g).^(m+1)));
    B_r_J_g = -sigma_m.*Amp_r_J_g*sin(m'.*Theta);

    % one-component magnetic potential 
    Amp_Az_J_g = (A_J_m_g./m*R_w.*((r_g./R_i).^(m).*(R_w/R_i).^(m)+(R_w./r_g).^(m)));
    A_z_J_g = sigma_m.*Amp_Az_J_g*cos(m'.*Theta);

    % circumferential component
    Amp_theta_J_g = (A_J_m_g.*(-(r_g./R_i).^(m-1).*(R_w/R_i).^(m+1)+(R_w./r_g).^(m+1)));
    B_theta_J_g = sigma_m.*Amp_theta_J_g*cos(m'.*Theta);

end

%% MAGNETS' BACKING
if R_i == R_r
if R_s>R_m % INRUNNER stator core region field computation
    
    G_Rm_J = (2*A_J_m_g.*(R_r./R_w).^(m-1))./(1-(R_ie/R_r).^(2*m));
    
    Amp_AzR_m = (R_r*G_Rm_J./m.*((r_ext./R_r).^(m)-(R_ie/R_r).^(m).*(R_ie./r_ext).^(m)));
    A_z_R_J = sigma_m.*Amp_AzR_m*cos(m'.*Theta);
    
    Amp_Rr_m = (G_Rm_J.*((r_ext./R_r).^(m-1)-(R_ie/R_r).^(m-1).*(R_ie./r_ext).^(m+1)));
    B_Rr_J = sigma_m.*Amp_Rr_m*sin(m'.*Theta);

    Amp_Rtheta_m = -(G_Rm_J.*((r_ext./R_r).^(m-1)+(R_ie/R_r).^(m-1).*(R_ie./r_ext).^(m+1)));
    B_Rtheta_J = sigma_m.*Amp_Rtheta_m*cos(m'.*Theta);
    
else % OUTRUNNER stator core region field computation
    
    G_Rm_out_J = (2*A_J_m_g.*(R_w./R_r).^(m+1))./((R_r/R_ie).^(2*m)-1);
    
    Amp_AzR_m = (R_r*G_Rm_out_J./m.*((r_ext./R_ie).^(m).*(R_r/R_ie).^(m)-(R_r./r_ext).^(m)));
    A_z_R_J = sigma_m.*Amp_AzR_m*cos(m'.*Theta);
    
    Amp_Rr_m = (G_Rm_out_J.*((r_ext./R_ie).^(m-1).*(R_r/R_ie).^(m+1)-(R_r./r_ext).^(m-1)));
    B_Rr_J = sigma_m.*Amp_Rr_m*sin(m'.*Theta);

    Amp_Rtheta_m = -(G_Rm_out_J.*((r_ext./R_ie).^(m-1).*(R_r/R_ie).^(m+1)+(R_r./r_ext).^(m-1)));
    B_Rtheta_J = sigma_m.*Amp_Rtheta_m*cos(m'.*Theta);
end
end

%% WINDING BACKIG
if R_ws == R_s
    if R_s>R_m % INRUNNER stator core region field computation
    
        G_wb_m_J = (A_pl_J_w + A_mi_J_w.*(R_w./R_s).^(m+1) + m.*A_mJ*R_s)./((R_s/R_se).^(2*m)-1);
        
        if p == 2
            G_wb_m_J(1) = (A_pl_J_w(1) + A_mi_J_w(1).*(R_w./R_s).^(3) + +2.*A_mJ(1).*R_s.*log(R_s))./((R_s/R_se).^(4)-1);
        end

        Amp_Az_wb_m = (R_s*G_wb_m_J./m.*((r_s./R_se).^(m).*(R_s/R_se).^(m)-(R_s./r_s).^(m)));
        A_z_wb_J = sigma_m.*Amp_Az_wb_m*cos(m'.*Theta);

        Amp_wb_r_m = (G_wb_m_J.*((r_s./R_se).^(m-1).*(R_s/R_se).^(m+1)-(R_s./r_s).^(m+1)));
        B_wb_r_J = sigma_m.*Amp_wb_r_m*sin(m'.*Theta);

        Amp_wb_theta_m = -(G_wb_m_J.*((r_s./R_se).^(m-1).*(R_s/R_se).^(m+1)+(R_s./r_s).^(m+1)));
        B_wb_theta_J = sigma_m.*Amp_wb_theta_m*cos(m'.*Theta);
        
    else % OUTRUNNER stator core region field computation

        G_wb_m_out_J = (A_pl_J_w.*(R_s/R_w).^(m-1) + A_mi_J_w + m.*A_mJ*R_s)./(1-(R_se/R_s).^(2*m));
        
        if p == 2
            G_wb_m_out_J(1) = (A_pl_J_w(1).*(R_s/R_w) + A_mi_J_w(1) +2.*A_mJ(1).*R_s.*log(R_s))./(1-(R_se/R_s).^(4));
        end

        Amp_Az_wb_m = (R_s*G_wb_m_out_J./m.*((r_s./R_s).^(m)-(R_se./r_s).^(m).*(R_se/R_s).^(m)));
        A_z_wb_J = sigma_m.*Amp_Az_wb_m*cos(m'.*Theta);

        Amp_wb_r_m = (G_wb_m_out_J.*((r_s./R_s).^(m-1)-(R_se./r_s).^(m+1).*(R_se/R_s).^(m-1)));
        B_wb_r_J = sigma_m.*Amp_wb_r_m*sin(m'.*Theta);

        Amp_wb_theta_m = -(G_wb_m_out_J.*((r_s./R_s).^(m-1)+(R_se./r_s).^(m+1).*(R_se/R_s).^(m-1)));
        B_wb_theta_J = sigma_m.*Amp_wb_theta_m*cos(m'.*Theta);
    end

    
else
    if R_s>R_m % INRUNNER stator core region field computation
    
        G_wb_m_J = (2*A_J_m_bw.*(R_ws./R_s).^(m+1))./((R_s/R_se).^(2*m)-1);

        Amp_Az_wb_m = (R_s*G_wb_m_J./m.*((r_s./R_se).^(m).*(R_s/R_se).^(m)-(R_s./r_s).^(m)));
        A_z_wb_J = sigma_m.*Amp_Az_wb_m*cos(m'.*Theta);

        Amp_wb_r_m = (G_wb_m_J.*((r_s./R_se).^(m-1).*(R_s/R_se).^(m+1)-(R_s./r_s).^(m+1)));
        B_wb_r_J = sigma_m.*Amp_wb_r_m*sin(m'.*Theta);

        Amp_wb_theta_m = -(G_wb_m_J.*((r_s./R_se).^(m-1).*(R_s/R_se).^(m+1)+(R_s./r_s).^(m+1)));
        B_wb_theta_J = sigma_m.*Amp_wb_theta_m*cos(m'.*Theta);

    else % OUTRUNNER stator core region field computation

        G_wb_m_out_J = (2*A_J_m_bw.*(R_s./R_ws).^(m-1))./(1-(R_se/R_s).^(2*m));

        Amp_Az_wb_m = (R_s*G_wb_m_out_J./m.*((r_s./R_s).^(m)-(R_se./r_s).^(m).*(R_se/R_s).^(m)));
        A_z_wb_J = sigma_m.*Amp_Az_wb_m*cos(m'.*Theta);

        Amp_wb_r_m = (G_wb_m_out_J.*((r_s./R_s).^(m-1)-(R_se./r_s).^(m+1).*(R_se/R_s).^(m-1)));
        B_wb_r_J = sigma_m.*Amp_wb_r_m*sin(m'.*Theta);

        Amp_wb_theta_m = -(G_wb_m_out_J.*((r_s./R_s).^(m-1)+(R_se./r_s).^(m+1).*(R_se/R_s).^(m-1)));
        B_wb_theta_J = sigma_m.*Amp_wb_theta_m*cos(m'.*Theta);
    end
    
end

%% FIGURES


if contains(mapping, 'yes')
% Create polar data
[r_g,t] = meshgrid(r_g,theta);
[r_w,t_w] = meshgrid(r_w,theta);
[r_bw,t_bw] = meshgrid(r_bw,theta);
[r_s,t_s] = meshgrid(r_s,theta);
[r_ext,t_ext] = meshgrid(r_ext,theta);
x = r_g.*sin(t);
y = r_g.*cos(t);
x_w = r_w.*sin(t_w);
y_w = r_w.*cos(t_w);
x_bw = r_bw.*sin(t_bw);
y_bw = r_bw.*cos(t_bw);
x_s = r_s.*sin(t_s);
y_s = r_s.*cos(t_s);
x_ext = r_ext.*sin(t_ext);
y_ext = r_ext.*cos(t_ext);

% Flux density Norm within the different domains
Norm_BgJ = sqrt(B_r_J_g'.^2+B_theta_J_g'.^2);
Norm_BwJ = sqrt(B_r_J_w'.^2+B_theta_J_w'.^2);
Norm_BbwJ = sqrt(B_r_J_bw'.^2+B_theta_J_bw'.^2);
Norm_BwbJ = sqrt(B_wb_r_J'.^2+B_wb_theta_J'.^2);
if R_i == R_r
    Norm_BRJ = sqrt(B_Rr_J'.^2+B_Rtheta_J'.^2);
    Levels = linspace(min([min(min(A_z_J_bw)) min(min(A_z_J_g)) min(min(A_z_J_w)) min(min(A_z_R_J)) min(min(A_z_wb_J))]),max([max(max(A_z_J_bw)) max(max(A_z_J_g)) max(max(A_z_J_w)) max(max(A_z_R_J)) max(max(A_z_wb_J))]),20);
else
    Levels = linspace(min([min(min(A_z_J_bw)) min(min(A_z_J_g)) min(min(A_z_J_w)) min(min(A_z_wb_J))]),max([max(max(A_z_J_bw)) max(max(A_z_J_g)) max(max(A_z_J_w)) max(max(A_z_wb_J))]),20);
end


B_min =  min( [min(min(Norm_BbwJ))  min(min(Norm_BgJ))  min(min(Norm_BwbJ))  min(min(Norm_BwJ))]);
B_max =  max( [max(max(Norm_BbwJ))  max(max(Norm_BgJ))  max(max(Norm_BwbJ))  max(max(Norm_BwJ))]);
B_levels = linspace(B_min,B_max,100);
figure('Renderer','Painters');
hold on;
contourf(x,y,Norm_BgJ,B_levels,'LineStyle','none');
contour(x, y,A_z_J_g',Levels,'LineColor','k','linewidth',1.2)  
contourf(x_w,y_w,Norm_BwJ,B_levels,'LineStyle','none');
contour(x_w, y_w,A_z_J_w',Levels,'LineColor','k','linewidth',1.2)
contourf(x_bw,y_bw,Norm_BbwJ,B_levels,'LineStyle','none');
contour(x_bw, y_bw,A_z_J_bw',Levels,'LineColor','k','linewidth',1.2)
contourf(x_s,y_s,Norm_BwbJ,B_levels,'LineStyle','none');
contour(x_s, y_s,A_z_wb_J',Levels,'LineColor','k','linewidth',1.2)
if R_i == R_r
    contourf(x_ext,y_ext,Norm_BRJ,B_levels,'LineStyle','none');
    contour(x_ext, y_ext,A_z_R_J',Levels,'LineColor','k','linewidth',1.2) 
end
plot(R_ie*sin(linspace(theta(1),theta(end),100)),R_ie*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot(R_r*sin(linspace(theta(1),theta(end),100)),R_r*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot(R_m*sin(linspace(theta(1),theta(end),100)),R_m*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot([R_r*sin(theta(1)+alpha_p*pi/(2*p):pi/(p):theta(end));R_m*sin(theta(1)+alpha_p*pi/(2*p):pi/(p):theta(end))],[R_r*cos(theta(1)+alpha_p*pi/(2*p):pi/(p):theta(end));R_m*cos(theta(1)+alpha_p*pi/(2*p):pi/(p):theta(end))],'linewidth',0.8,'color','k');
plot([R_r*sin(theta(1)+pi/p-alpha_p*pi/(2*p):pi/(p):theta(end));R_m*sin(theta(1)+pi/p-alpha_p*pi/(2*p):pi/(p):theta(end))],[R_r*cos(theta(1)+pi/p-alpha_p*pi/(2*p):pi/(p):theta(end));R_m*cos(theta(1)+pi/p-alpha_p*pi/(2*p):pi/(p):theta(end))],'linewidth',0.8,'color','k');
plot(R_s*sin(linspace(theta(1),theta(end),100)),R_s*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot(R_se*sin(linspace(theta(1),theta(end),100)),R_se*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot(R_w*sin(linspace(theta(1),theta(end),100)),R_w*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot(R_ws*sin(linspace(theta(1),theta(end),100)),R_ws*cos(linspace(theta(1),theta(end),100)),'linewidth',0.8,'color','k');
plot([R_ws*sin(theta(1)+pi/(6*p):pi/(3*p):theta(end)-pi/(6*p));R_w*sin(theta(1)+pi/(6*p):pi/(3*p):theta(end)-pi/(6*p))],[R_ws*cos(theta(1)+pi/(6*p):pi/(3*p):theta(end)-pi/(6*p));R_w*cos(theta(1)+pi/(6*p):pi/(3*p):theta(end)-pi/(6*p))],'linewidth',0.8,'color','k');
set(gca,'visible','off');
colormap(jet)
c = colorbar;
c.Label.String = 'Flux density norm [T]';
caxis([B_min B_max])
% Turn off axes and set square aspect ratio
axis off
axis image
title('Flux density norm map and isopotential lines in the whole 2D domain (optimal formulation)')
end
