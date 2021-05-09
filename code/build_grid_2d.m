% This script is the 2-dimensional version of build_grid.m
%   usually used before 2D Fourier-Bessel decomposition

% For these global variables
global x2d y2d tht2d r2d alpha delta_x r_max in2d out2d

load(fullfile(code_path, 'alpha.mat'));
[x2d, y2d] = ndgrid(x, x);
[tht2d, r2d] = cart2pol(x2d, y2d);
delta_x = dx;
r_max = radius;
in2d = r2d < radius;
out2d = r2d >= radius;
