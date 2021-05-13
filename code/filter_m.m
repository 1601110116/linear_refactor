function f_filtered = filter_m(f, m, l_max)
% This function filters a 2D field using Fourier-Bessel series
%  should be used after build_grid_2d.m
%
% Input:
%  f: a scalar field in the original Cartesian grid,
%     either be 2D(x-y) or 3D(x-y-z)
%  m: filter f to only maintain the azimuthal mode number of m
%  l_max: keep the Bessel series up to l_max

global out2d

f_filtered = zeros(size(f));
for l = 1: l_max
	[coeff, f_ml] = fbt_ml(f, m, l);
	f_filtered = f_filtered + f_ml;
end

% The field outside the cylinder is considered to be m=0 
%  but does not belong to any l
out3d = out2d & true(size(f));
if m == 0
	f_filtered(out3d) = f(out3d);
end
