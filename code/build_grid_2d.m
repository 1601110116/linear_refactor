% This script is the 2-dimensional version of build_grid.m
%   usually used before 2D Fourier-Bessel decomposition

% For these global variables
global x2d y2d tht2d r2d alpha dx

load(fullfile(code_path, 'alpha.mat'));
[x2d, y2d] = ndgrid(x, x);
[tht2d, r2d] = cart2pol(x2d, y2d);
