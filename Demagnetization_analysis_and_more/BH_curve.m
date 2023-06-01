clearvars
close all
clc

%% IMPORTING AVAILBLE DATA
% Check file names in PM-data folder
filename = 'N30AH@220C';

CELL = readcell('PM_data.csv');

if any(any(strcmp(CELL,filename)))
    msg = 'Material data already stored.';
    error(msg)
else % Upload PM data




B=importdata(['PM_data\' filename '.tab'],'\t');



%% PROCESSING DATA (TWO SEGMENTS CURVE)
% Operations on the normal curve
BH = B.data;
B_r = BH(end,2);

v = interp1(BH(:,2),BH(:,1),linspace(BH(1,2),BH(end,2),1000));             % improve data discretization
NEW = [linspace(BH(1,2),BH(end,2),1000)' v'];                              % dense Bh curve
mu_NEW = (diff(NEW(:,1))./diff(NEW(:,2)))./(4*pi*1e-7);                    % relative permeability curve
Hc_NEW = interp1(NEW(:,1),NEW(:,2),0);                                     % coercive field (H-field at zero flux density)
mu_r_test = mu_NEW(find(mu_NEW> 1,1,'last'));
variation = abs(mu_NEW-mu_r_test)./mu_r_test*100;
H_int = NEW(find(variation> 50,1,'last'),2);                               % found when the recoil permeability changes more than 50%

x = [abs(diff(NEW(:,2)))>0;true];
x(find(x==1,1,'first')) = 0;
NEW = [[NEW(1,1);NEW(x,1)] , [NEW(1,2);NEW(x,2)]];


B20 = interp1(NEW(:,2),NEW(:,1),0.2*Hc_NEW);                               % B-value at 20% coercive field
B70 = interp1(NEW(:,2),NEW(:,1),0.7*Hc_NEW);                               % B-value at 70% coercive field
mu_r_ARNOLD = (B20-B70)/(0.2*Hc_NEW-0.7*Hc_NEW)/(4*pi*1e-7);               % recoil permeability according to Arnold definition

B095_line = 0.95*NEW(end,1)+NEW(:,2)*mu_r_ARNOLD*4*pi*1e-7;                % recoil line starting at 0.95 the Br
Brec_line = NEW(end,1)+NEW(:,2)*mu_r_ARNOLD*4*pi*1e-7;                     % recoil line
ind = sign(NEW(:,1)-B095_line);                                            % intersection between 0.95 recoil line and real BH curve
H_knee = NEW(find(ind> -1,1,'first'),2);                                   % find Hknee from the previous intersection
B_knee = interp1(NEW(:,2),Brec_line,H_knee);                               % find Bknee from the recoil line at Hknee
B_knee2 = NEW(find(ind> -1,1,'first'),1);                                  % find Bknee from the BH curve at Hknee
mu_knee = (B_knee2-NEW(1,1))/(H_knee-NEW(1,2))/(4*pi*1e-7);
q = B_knee2-mu_knee*(4*pi*1e-7)*H_knee;                                                % Intersection of knee-line with H=0 axis

start_knee = 1.01*NEW(1,2);                                                 % start point of knee-line
end_knee = 0.99*H_knee;                                                     % end point of knee-line
BH_knee_line = [mu_knee*(4*pi*1e-7)*(linspace(start_knee,end_knee,1000)-start_knee)+ mu_knee*(4*pi*1e-7)*start_knee+q ; linspace(start_knee,end_knee,1000)];
B0_line = linspace(start_knee,end_knee,1000)*mu_r_ARNOLD*4*pi*1e-7;
% B_demag = B0_line(find(B0_line-BH_knee_line(1,:)<0,1,'first'));
% B_demag = B_demag - B_knee2;
B_demag = H_knee*mu_r_ARNOLD*4*pi*1e-7-B_knee2; 
% BH_knee_line = [mu_knee*(4*pi*1e-7)*(linspace(NEW(1,2),H_knee,1000)-NEW(1,2))+NEW(1,1) ; linspace(NEW(1,2),H_knee,1000)];
Hc_app = interp1(BH_knee_line(1,:),BH_knee_line(2,:),0);

% figure
hold on
plot(NEW(:,2),NEW(:,1))
plot(NEW(:,2),Brec_line)
plot(BH_knee_line(2,:),BH_knee_line(1,:))
plot(linspace(start_knee,0,1000),linspace(start_knee,0,1000)*mu_r_ARNOLD*4*pi*1e-7)

% data structured as {'PM_material','B_r','mu_rec','B_dem'}
data = {filename,B_r,round(mu_r_ARNOLD,5),round(B_knee2,5),round(B_demag,5)};
fileID = 'PM_data.csv';
writecell(data,fileID,'WriteMode','append')
end
% writecell({'PM_material','B_r','mu_rec','B_knee','B_dem'},fileID)



% Hc = interp1(BH(:,2),BH(:,1),0);