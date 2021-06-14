% This script is the 2-dimensional version of build_grid.m
%   usually used before 2D Fourier-Bessel decomposition

% For these global variables
global x2d y2d tht2d r2d alpha delta_x r_max in2d out2d ...
	bessel_cores fourier_cores

%---input---
m_max = 20;
l_max = 30;
cores_file = fullfile('data', 'fb_cores.mat');
%-----------
load(fullfile(code_path, 'alpha.mat'));
[x2d, y2d] = ndgrid(x, x);
[tht2d, r2d] = cart2pol(x2d, y2d);
delta_x = dx;
r_max = radius;
in2d = r2d < radius;
out2d = r2d >= radius;
need_core_gen = 1;
if exist(cores_file, 'file')
	load(cores_file, 'l_up_to', 'm_up_to');
	if l_up_to >= l_max & m_up_to >= m_max
		need_core_gen = 0;
	end
end
if need_core_gen
	disp('calculating Fourier-Bessel integral kernels');
	bessel_cores = cell(m_max+1, l_max);
	fourier_cores = cell(m_max+1, 1);
	for m = 0: m_max
	    disp(['m = ', num2str(m), ' of ', num2str(m_max)]);
	    for l = 1: l_max
	        lambda_ml = alpha(m+1, l) / radius;
	        bessel_cores{m+1, l} = besselj(m, lambda_ml * r2d);
	    end
	    fourier_cores{m+1} = exp(-1j * m * tht2d);
	end
	l_up_to = l_max;  m_up_to = m_max;
	save(cores_file, 'bessel_cores', 'fourier_cores', ...
		'l_up_to', 'm_up_to');
else
	load(cores_file);
end
