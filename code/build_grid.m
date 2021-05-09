% This script generates grids used for further analysis
%  xc: x of the cartesian grid points
%  xp: x of the polar grid points
%  rp: radial position of the polar grid points
%  thtp: azimuthal position of the polar grid points
%  zp: axial position of the polar grid points

global xc yc zc xp yp zp rp thtp zp thtc rc zc

[xc, yc, zc] = ndgrid(x, x, z);
nr = nx;
ntht = 6*nr;
r = linspace(0, x_max, nr);
dr = r(2)-r(1);
tht = linspace(0, 2*pi, ntht);
dtht = tht(2) - tht(1);
[rp, thtp, zp] = ndgrid(r, tht, z);
[xp, yp, zp] = pol2cart(thtp, rp, zp);
[thtc, rc, zc] = cart2pol(xc, yc, zc);
