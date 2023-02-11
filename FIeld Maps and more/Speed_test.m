runs = 10000;

tic
for ii = 1:runs
    Bemf_torque_constants
end

delta_t = toc;

X = ["Script 1 ran",runs,' in ',delta_t,' taking approx.',delta_t/runs, 'per single run'];
disp(X)