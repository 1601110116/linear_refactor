function [vr, vtht] = vec2pol_cartgrid(vx, vy)
% This function decomposits a vector in cartician grid to r-tht directions
global rc thtc zc

vr = vx .* cos(thtc) + vy .* sin(thtc);
vtht = -vx .* sin(thtc) + vy .* cos(thtc);
