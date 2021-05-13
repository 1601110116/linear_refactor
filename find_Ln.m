% This script draws different profiles with different Ln
%  to find the Ln that is close to the experimental density profile

clear all;  close all;
load('parameters.mat');
addpath(code_path);

%---input---
idiag = 70;
l_max = 30;
font_size = 15;
fig_position = [50, 50, 300, 300];
line_width = 2;
%-----------

data_file = fullfile('data', ['dat', sprintf('%4.4d', idiag)]);
load(data_file, 'den');

% find the profile using interpolation
build_grid;
denp = field2pol(den);
den_profile_interp = squeeze(zonal_average(denp));

% find the profile using Fourier-Bessel decomposition
build_grid_2d;
den_profile_fbt = fbt_mr(den, r, 0, l_max);

fig = figure();
rX = rhos0 * r;
plot(rX, den_profile_interp, 'LineWidth', line_width);  hold on;
plot(rX, den_profile_fbt, 'LineWidth', line_width);  hold off;
legend('interp', 'FBt');
set(gca, 'fontSize', font_size);
