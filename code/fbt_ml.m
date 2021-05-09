function [coeff, f_ml]  = fbt_ml(f, m, l)
% This function calculates the coefficient of a Fourier-Bessel component
%  should be used after build_grid_2d.m
%  all variables here are normalized
%
% Input:
%  f: a 2D field in the original Cartesian grid
%  m: azimuthal mode number
%  l: l-th Bessel series
% Output:
%  coeff: the Fourier-Bessel coeffcient of f
%  f_ml: back-transformed m-l component of f

global tht2d r2d alpha dx

% l-th positive zero of Bessel function of order m 
%  scaled by the invert of radius
lambda_ml = alpha(m+1, l) / radius;
% J_m(\lambda_{ml} r)
bessel_core = besselj(m, lambda_ml * r2d);
% e^{-im\theta}
fourier_core = exp(-1j * m * tht2d);
% The factor in front of the integration
fac = 1 / (pi * radius^2 * (besselj(m+1, alpha(m+1, l))).^2);
coeff = fac * dx^2 * sum(sum(f .* bessel_core .* fourier_core)
if m == 0
	f_ml = coeff * bessel_core;
else
	f_ml = 2 * real(coeff * conj(fourier_core) .* bessel_core);
end

