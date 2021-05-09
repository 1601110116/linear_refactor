function h = ddz(g)
% This function calculates the z-direction derivative of g, times dt
%  the resulting term h excludes boundaries
global ddzdt
h = ddzdt * (g(2:end-1, 2:end-1, 3:end) - g(2:end-1, 2:end-1, 1:end-2));
