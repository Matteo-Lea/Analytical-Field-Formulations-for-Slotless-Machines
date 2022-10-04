clc;clearvars;

el = 100;
% x structure [ Halbach , Br, Brs, mu_r, alpha_p, teta_m_end, p ]
x = [round(rand(el,1)), 0.9+0.5*rand(el,1), 0.9+0.5*rand(el,1), ones(el,1), rand(el,1)*(0.9-0.1)+0.1, rand(el,1).*90, round(rand(el,1).*(25-1)+1)];
Halbach = x(:,1); % 0 for even segments and 1 for odd segmens arrays
Halbach_1 = Halbach+2*round((1+rand(length(Halbach),1))); % number of segments per pole
B_r = x(:,2); % remanent flux density [T] (mid-magnet)
B_rs = x(:,3); % remanent flux density [T] (side-magnets)
mu_r = x(:,4); % PM recoil permeability [-]
alpha_p =  x(:,5); % mid-magnet to pole ratio [-]
alpha_p(Halbach_1==2) = alpha_p(Halbach_1==2).*round(rand(nnz(Halbach_1==2),1)); % some special Halbach are included as well
alpha_p1 = rand(el,1).*(1-x(:,5))+x(:,5); % side-magnet + mid-magnet to pole ratio [-]
theta_m_end = x(:,6); % orientation angle end-side-magnet [deg]
theta_m_side = rand(el,1).*(90-theta_m_end)+theta_m_end; % orientation angle side-magnet [deg]

if nnz(~(mod(Halbach_1,2)~=0 + Halbach==0))
    commandwindow
    error([' -)You have set an odd number of segments per pole', ...
          ' but \n selected an even-segment Halbach array.  ', ... 
          ' \n  change either of the following values:', ...
          ' \n Halbach or Halbach_1'])
elseif nnz(~(mod(Halbach_1,2)==0 + Halbach==1))
    commandwindow
    error([' -)"You have set an even number of segments per pole \n',...
           ' but selected an odd-segment Halbach array.  \n',...
           ' change either of the following values:\n',...
           ' Halbach or Halbach_1"']);
end


alpha_p1(alpha_p==0) =1;
alpha_p1(Halbach_1 == 2) = alpha_p(Halbach_1 == 2); % side-magnets + mid-magnet to pole ratio [-]
alpha_p1(Halbach_1 == 3) = 1; 


theta_m_side = theta_m_side*pi/180;
theta_m_end = theta_m_end*pi/180;
p =x(:,7); % pole pairs
delta_m = (alpha_p+alpha_p1)*pi./(4*p); % mid-point angle of the side-magnet
delta_m1 = (alpha_p1+1)*pi./(4*p); % mid-point angle of the end-side-magnet

%% useful indices for series harmonics definition
m_PM = 10; % number of harmonics or series components tionfor the magnetization functions
x = 0:1:m_PM;
% x = linspace(0,40,m_PM+1);
n = p*(2*x+1); % existing harmonic series components (n=1,3,5,...)

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
M_1_0 = sin((n-1).*alpha_p.*pi./(2*p))./((n-1).*alpha_p.*pi./(2*p));
M_2_0 = sin((n+1).*alpha_p.*pi./(2*p))./((n+1).*alpha_p.*pi./(2*p));
% Side-magnets coefficients
M_1_2 = cos((n-1).*alpha_p.*pi./(2.*p)+theta_m_side+delta_m)./((n-1).*alpha_p.*pi./(2*p));
M_2_2 = cos((n+1).*alpha_p.*pi./(2*p)-theta_m_side-delta_m)./((n+1).*alpha_p.*pi./(2*p));
M_3_2 = cos((n-1).*alpha_p1.*pi./(2*p)+theta_m_side+delta_m)./((n-1).*alpha_p1.*pi./(2*p));
M_4_2 = cos((n+1).*alpha_p1.*pi./(2*p)-theta_m_side-delta_m)./((n+1).*alpha_p1.*pi./(2*p));
% End-Side-magnets coefficients
M_1_3 = cos((n-1).*alpha_p1.*pi./(2*p)+theta_m_end+delta_m1)./((n-1).*alpha_p1.*pi./(2*p));
M_2_3 = cos((n+1).*alpha_p1.*pi./(2*p)-theta_m_end-delta_m1)./((n+1).*alpha_p1.*pi./(2*p));
M_3_3 = cos((n-1).*pi./(2*p)+theta_m_end+delta_m1)./((n-1).*pi./(2*p));
M_4_3 = cos((n+1).*pi./(2*p)-theta_m_end-delta_m1)./((n+1).*pi./(2*p));

% correction of singularities
M_1_0(p==1,1) = 1; 
M_1_2(p==1,1) = sin(theta_m_side(p==1)+delta_m(p==1));
M_3_2(p==1,1) = sin(theta_m_side(p==1)+delta_m(p==1));
M_1_3(p==1,1) = sin(theta_m_end(p==1)+delta_m1(p==1));
M_3_3(p==1,1) = sin(theta_m_end(p==1)+delta_m1(p==1));

    
% Mid-magnets coefficients
M_1_0(alpha_p ==0,:) = 0;
M_2_0(alpha_p ==0,:) = 0;
% Side-magnets coefficients
M_1_2(alpha_p ==0,:) = cos(theta_m_side(alpha_p ==0)+delta_m(alpha_p ==0))./((n(alpha_p ==0,:)-1)*pi./(2*p(alpha_p ==0)));
M_2_2(alpha_p ==0,:) = cos(theta_m_side(alpha_p ==0)+delta_m(alpha_p ==0))./((n(alpha_p ==0,:)+1)*pi./(2*p(alpha_p ==0)));
M_3_2(alpha_p ==0,:) = cos((n(alpha_p ==0,:)-1)*pi./(2*p(alpha_p ==0))+theta_m_side(alpha_p ==0)+delta_m(alpha_p ==0))./((n(alpha_p ==0,:)-1)*pi./(2*p(alpha_p ==0)));
M_4_2(alpha_p ==0,:) = cos((n(alpha_p ==0,:)+1)*pi./(2*p(alpha_p ==0))-theta_m_side(alpha_p ==0)-delta_m(alpha_p ==0))./((n(alpha_p ==0,:)+1)*pi./(2*p(alpha_p ==0)));
% End-side-magnets coefficients
M_1_3(alpha_p ==0,:) = 0;
M_2_3(alpha_p ==0,:) = 0;
M_3_3(alpha_p ==0,:) = 0;
M_4_3(alpha_p ==0,:) = 0;

% correction of singularities
M_1_2( ~(p==1+alpha_p==0),1) = 0;
M_3_2( ~(p==1+alpha_p==0),1) = sin(theta_m_side(p(alpha_p==0)==1)+delta_m(p(alpha_p==0)==1));
