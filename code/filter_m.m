function f_filtered = filter_m(f, m, l_max)
% This function filters a 2D field using Fourier-Bessel series
%  should be used after build_grid_2d.m
%
% Input:
%  f: a 2D field in the original Cartesian grid
%  m: filter f to only maintain the azimuthal mode number of m
%  l_max: keep the Bessel series up to l_max

global tht2d r2d alpha dx

f_filtered = zeros(size(f));
for l = 1: l_max
	[coeff, f_ml] = fbt_ml(f, m, l);
	f_filtered = f_filtered + f_ml;
end
