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

Ln = 3; %1.85;  % cm
N0 = 15e12;  % cm^{-3}
%-----------

data_file = fullfile('data', ['dat', sprintf('%4.4d', idiag)]);
load(data_file, 'den');

% find the profile using interpolation
build_grid;
denp = field2pol(den);
den_profile_interp = denref * squeeze(zonal_average(denp));

% find the profile using Fourier-Bessel decomposition
build_grid_2d;
den_profile_fbt = fbt_mr(den, r, 0, l_max);
den_profile_fbt = denref * mean(den_profile_fbt, 3);

% Canonical profile, n=N0exp(-0.5*(r/Ln)^2)
Ln = Ln * cm;
den_profile_canonical = N0 * exp(-0.5*(r/Ln).^2);

fig = figure();
rX = rhos0 * r;
plot(rX, den_profile_interp, 'r-', 'LineWidth', line_width);  hold on;
plot(rX, den_profile_fbt, 'k+')
plot(rX, den_profile_canonical, 'b--', 'LineWidth', line_width);  hold off;
xline(r_max * rhos0, 'g--', 'LineWidth', line_width)
legend('interp', 'FBt', ...
	['N0=', sprintf('%.2e', N0), ',Ln=', sprintf('%.2f', Ln/cm)]);
xlabel('r (cm)');
ylabel('$$n\ \mathrm{\left(cm^{-3}\right)}$$', 'interpreter', 'latex');
set(gca, 'fontSize', font_size);

