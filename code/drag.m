function h = drag(g)
% This function calculates the decay term of g
%  the resulting term h excludes boundaries
global momentum_dampdt

h = -momentum_dampdt * g(2:end-1, 2:end-1, 2:end-1);
