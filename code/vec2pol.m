function [vr, vtht] = vec2pol(vx, vy)
% This function converts the grid of a vector in cartician grid to polar grid
global xc yc zc xp yp zp thtp

vxp = interp3(yc, xc, zc, vx, yp, xp, zp);
vyp = interp3(yc, xc, zc, vy, yp, xp, zp);

vr = vxp .* cos(thtp) + vyp .* sin(thtp);
vtht = -vxp .* sin(thtp) + vyp .* cos(thtp);
