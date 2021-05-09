function h = vdeconvect(g)
% This function calculates the electron diamagnetic convection of g
%  the resulting term h excludes boundaries
global vdex vdey convxdt

% This is equivalent to -[pe, h], [ , ] is the Poisson bracket
h = convxdt .* (vdex(3:end, 2:end-1, 2:end-1) .* g(3:end, 2:end-1, 2:end-1) ...
	- vdex(1:end-2, 2:end-1, 2:end-1) .* g(1:end-2, 2:end-1, 2:end-1)) + ...
	convxdt .* (vdey(2:end-1, 3:end, 2:end-1) .* g(2:end-1, 3:end, 2:end-1) ...
	- vdey(2:end-1, 1:end-2, 2:end-1) .* g(2:end-1, 1:end-2, 2:end-1));
