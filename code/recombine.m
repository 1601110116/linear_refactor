function h = recombine(g)
% This function calculates the decay term of g
%  the resulting term h excludes boundaries
global den_dampdt

h = -den_dampdt * g(2:end-1, 2:end-1, 2:end-1);
