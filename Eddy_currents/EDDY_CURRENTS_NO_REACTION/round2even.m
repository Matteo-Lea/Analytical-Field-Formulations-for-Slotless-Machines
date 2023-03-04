function [N] = round2even(N)
%round to nearest even integer
if mod(N,2)<1 
N = fix(N); 
else 
N =fix(N) + 1; 
end
end

