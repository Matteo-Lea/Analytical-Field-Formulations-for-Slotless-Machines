function [M_r_n,M_theta_n] = Magnetization(el,m_PM,x,Const)
% Definition of the magnetization distributionfor parallel magnetized
% magnets
mu_0 = 4*pi*1e-7; % air permeability

alpha_p1 = x(:,1);
theta_m_side = Const(:,6)*pi/180;
theta_m_end = Const(:,5)*pi/180;
delta_m = (x(:,1) +alpha_p1)*pi./(4*x(:,9)); % mid-point angle of the side-magnet
delta_m1 = pi./(2*x(:,9)); % adjustment of mid-point angle of the 
                                      % end-side-magnet for even-segmetns Halbach 
B_r = Const(:,3);
B_rs = Const(:,4);
                                     
                                        

%% useful indices for series harmonics definition
xi = 0:1:m_PM;
% x = linspace(0,40,m_PM+1);
n = x(:,9)*(2*xi+1); % existing harmonic series components (n=1,3,5,...)

%% Magnetization distribution

% PARALLEL MAGNETIZED PM
% Coefficients for defining the magnetization distribution
M_1_0(el,m_PM+1) = 0;
M_2_0(el,m_PM+1) = 0;
M_1_2(el,m_PM+1) = 0;
M_2_2(el,m_PM+1) = 0;
M_3_2(el,m_PM+1) = 0;
M_4_2(el,m_PM+1) = 0;
M_1_3(el,m_PM+1) = 0;
M_2_3(el,m_PM+1) = 0;
M_3_3(el,m_PM+1) = 0;
M_4_3(el,m_PM+1) = 0;



% Mid-magnets coefficients
M_1_0 = sin((n-1).*x(:,1).*pi./(2*x(:,9)))./((n-1).*x(:,1).*pi./(2*x(:,9)));
M_2_0 = sin((n+1).*x(:,1).*pi./(2*x(:,9)))./((n+1).*x(:,1).*pi./(2*x(:,9)));
% Side-magnets coefficients
M_1_2 = cos((n-1).*x(:,1).*pi./(2.*x(:,9))+theta_m_side+delta_m)./((n-1).*x(:,1).*pi./(2*x(:,9)));
M_2_2 = cos((n+1).*x(:,1).*pi./(2*x(:,9))-theta_m_side-delta_m)./((n+1).*x(:,1).*pi./(2*x(:,9)));
M_3_2 = cos((n-1).*alpha_p1.*pi./(2*x(:,9))+theta_m_side+delta_m)./((n-1).*alpha_p1.*pi./(2*x(:,9)));
M_4_2 = cos((n+1).*alpha_p1.*pi./(2*x(:,9))-theta_m_side-delta_m)./((n+1).*alpha_p1.*pi./(2*x(:,9)));
% End-Side-magnets coefficients
M_1_3 = cos((n-1).*alpha_p1.*pi./(2*x(:,9))+theta_m_end+delta_m1)./((n-1).*alpha_p1.*pi./(2*x(:,9)));
M_2_3 = cos((n+1).*alpha_p1.*pi./(2*x(:,9))-theta_m_end-delta_m1)./((n+1).*alpha_p1.*pi./(2*x(:,9)));
M_3_3 = cos((n-1).*pi./(2*x(:,9))+theta_m_end+delta_m1)./((n-1).*pi./(2*x(:,9)));
M_4_3 = cos((n+1).*pi./(2*x(:,9))-theta_m_end-delta_m1)./((n+1).*pi./(2*x(:,9)));

% correction of singularities
M_1_0(x(:,9)==1,1) = 1; 
M_1_2(x(:,9)==1,1) = -sin(x(x(:,9)==1,8)+delta_m(x(:,9)==1));
M_3_2(x(:,9)==1,1) = -sin(x(x(:,9)==1,8)+delta_m(x(:,9)==1));
M_1_3(x(:,9)==1,1) = -sin(x(x(:,9)==1,7)+delta_m1(x(:,9)==1));
M_3_3(x(:,9)==1,1) = -sin(x(x(:,9)==1,7)+delta_m1(x(:,9)==1));

% Magnetization harmonics amplitude

M_r_n_par_mid = B_r./mu_0.*x(:,1).*(M_1_0+M_2_0);
M_theta_n_par_mid = -B_r/mu_0.*x(:,1).*(M_1_0-M_2_0);
M_r_n_par_end_side = B_rs/mu_0.*(alpha_p1.*(M_1_3-M_2_3)+(M_4_3-M_3_3));
M_theta_n_par_end_side = -B_rs/mu_0.*(alpha_p1.*(M_2_3+M_1_3)-(M_3_3+M_4_3));
M_r_n_par_side = B_rs/mu_0.*(x(:,1).*(M_1_2-M_2_2)+alpha_p1.*(M_4_2-M_3_2));
M_theta_n_par_side = -B_rs/mu_0.*(x(:,1).*(M_2_2+M_1_2)-alpha_p1.*(M_3_2+M_4_2));


    
M_r_n = M_r_n_par_mid+M_r_n_par_end_side+M_r_n_par_side;
M_theta_n= x(:,8).*(M_theta_n_par_mid+M_theta_n_par_end_side+M_theta_n_par_side);

M_r_n(x(:,9)==1,1) = M_r_n_par_mid(x(:,9)==1,1)-M_r_n_par_end_side(x(:,9)==1,1)-M_r_n_par_side(x(:,9)==1,1);
M_theta_n(x(:,9)==1,1)= (M_theta_n_par_mid(x(:,9)==1,1)-M_theta_n_par_end_side(x(:,9)==1,1)-M_theta_n_par_side(x(:,9)==1,1));
    
end

