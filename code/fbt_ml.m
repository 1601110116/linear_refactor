function [coeff, f_ml]  = fbt_ml(f, m, l)
% This function calculates the coefficient of a Fourier-Bessel component
%  should be used after build_grid_2d.m
%  all variables here are normalized
%
% Input:
%  f: a scalar field in the original Cartesian grid,
%	  either be 2D(x-y) or 3D(x-y-z)
%  m: azimuthal mode number
%  l: l-th Bessel series
% Output:
%  coeff: the Fourier-Bessel coeffcient of f
%  f_ml: back-transformed m-l component of f

global tht2d r2d alpha delta_x r_max in2d out2d fourier_cores ...
	bessel_cores

% l-th positive zero of Bessel function of order m 
%  scaled by the invert of radius
lambda_ml = alpha(m+1, l) / r_max;
% J_m(\lambda_{ml} r)
%bessel_core = besselj(m, lambda_ml * r2d);
bessel_core = bessel_cores{m+1, l};
% e^{-im\theta}
%fourier_core = exp(-1j * m * tht2d);
fourier_core = fourier_cores{m+1};
% The factor in front of the integration
fac = 1 / (pi * r_max^2 * (besselj(m+1, alpha(m+1, l))).^2);
% For the integration to get the coefficients, we set the variable outside
%  the cylinder to be zero to exclude their contributions
f_in = f .* in2d;
coeff = fac * delta_x^2 * sum(sum(f_in .* bessel_core .* fourier_core, 1), 2);

% coeff is only able to reproduce the data inside the cylinder
% f_ml outside the cylinder is considered to be not belong to any l
if m == 0
	f_ml = coeff .* bessel_core;
else
	f_ml = 2 * real(coeff .* conj(fourier_core) .* bessel_core);
end
out3d = out2d & true(size(f));
f_ml(out3d) = 0;

