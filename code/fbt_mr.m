function f_mr = fbt_mr(f, r_array, m, l_max)
% This function calculates the azimuthally Fourier-transformed radial
%   distribution of a single m component

% Input:
%  f: a 2D field in the original Cartesian grid
%  r_array: an array of radial positions, normalized
%  m: the azimuthal mode number to analyse
%  l_max: keep the radial Bessel series up to and include l_max

global r_max alpha

if max(r_array) > r_max
	warning(['Input r_array exceeds the radius r_max = ', num2str(r_max)]);
end

f_mr = zeros(size(r_array));
for l = 1: l_max
	% l-th positive zero of Bessel function of order m
	%  scaled by the invert of radius
	lambda_ml = alpha(m+1, l) / r_max;
	[coeff, f_ml] = fbt_ml(f, m, l);
	% calculate the Bessel component at r_array
	f_mr = f_mr + besselj(m, lambda_ml * r_array);
end
