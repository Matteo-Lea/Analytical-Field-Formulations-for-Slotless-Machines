runs = 2000;

tic
for ii = 1:runs
    PM_field_map
end

delta_t = toc;

X = ["Script 1 ran",runs,' in ',delta_t,' taking approx.',delta_t/runs, 'per single run'];
disp(X)