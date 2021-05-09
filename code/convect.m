function h = convect(g)
% This function calculates the EXB convection of g
%  the resulting term h excludes boundaries
global vEx vEy convxdt

% This is equivalent to [\phi, h], [ , ] is the Poisson bracket
h = convxdt .* (vEx(3:end, 2:end-1, 2:end-1) .* g(3:end, 2:end-1, 2:end-1) ...
	- vEx(1:end-2, 2:end-1, 2:end-1) .* g(1:end-2, 2:end-1, 2:end-1)) ...
	+ convxdt .* (vEy(2:end-1, 3:end, 2:end-1) .* g(2:end-1, 3:end, 2:end-1) ...
	- vEy(2:end-1, 1:end-2, 2:end-1) .* g(2:end-1, 1:end-2, 2:end-1));
