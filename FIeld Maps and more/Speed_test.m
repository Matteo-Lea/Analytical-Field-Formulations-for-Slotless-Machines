clearvars
clc
runs = 1000;
% Bemf_torque_constants
% Inductance
tic
for ii = 1:runs
    Inductance % change this to any script you would like to check the runtime of
end

delta_t = toc;


% fprintf('The script ran %s in %d taking approx. %a per single run .\n',runs,delta_t,delta_t/runs);
X = ['The script ran ',num2str(runs),' in ',num2str(delta_t),' taking approx. ',num2str(delta_t/runs), ' per single run'];
disp(X)