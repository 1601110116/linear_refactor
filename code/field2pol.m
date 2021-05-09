function h = field2pol(g)
% This function converts the grid of a field in cartecian grid to polar grid
global xc yc zc xp yp zp

h = interp3(yc, xc, zc, g, yp, xp, zp);

