
clear all;  close all;

load('parameters.mat');
addpath(code_path);

%---input---
idiag = 71;
zdiag = 12;
ydiag = length(xX) / 2;
l_max = 30;
m = 7;
fig_position = [50, 50, 1600, 600];
font_size = 15;
filter_n1 = 1;  % the parallel harmonic to keep. Keep all if <0
%-----------

build_grid_2d;
data_file = fullfile('data', ['dat', sprintf('%4.4d', idiag)]);
load(data_file);

if filter_n1 >= 0
	den_n = filter_z(den, filter_n1);
	phi_n = filter_z(phi, filter_n1);
else
	den_n = den;
	phi_n = phi;
end

den_n = den_n * denref;
phi_n = phi_n * Tref;

den_m = zeros(size(den));
phi_m = zeros(size(phi));

for iz = 2: nz+1
	den_m(:, :, iz) = filter_m(den_n(:, :, iz), m, l_max);
	phi_m(:, :, iz) = filter_m(phi_n(:, :, iz), m, l_max);
end
den_m = zbcs(den_m);
phi_m = zbcs(phi_m);

fig = figure();
set(fig, 'Position', fig_position);

subplot(2, 4, 1);
pcolor(rhos0*y2d, rhos0*x2d, den_n(:, :, zdiag));
colormap jet;  colorbar;  shading interp;
xlabel('y (cm)');
ylabel('x (cm)');
title('$$n\ \mathrm{\left(cm^{-3}\right)}$$', 'interpreter', 'latex');
set(gca, 'FontSize', font_size)

subplot(2, 4, 2);
pcolor(rhos0*y2d, rhos0*x2d, den_m(:, :, zdiag));
colormap jet;  colorbar;  shading interp;
xlabel('y (cm)');
ylabel('x (cm)');
title(['$$n_', num2str(m), '\ \mathrm{\left(cm^{-3}\right)}$$'], 'interpreter', 'latex');
set(gca, 'FontSize', font_size)

subplot(2, 4, 3);
pcolor(rhos0*y2d, rhos0*x2d, phi_n(:, :, zdiag));
colormap jet;  colorbar;  shading interp;
xlabel('y (cm)');
ylabel('x (cm)');
title('$$\phi\ \mathrm{\left(V\right)}$$', 'interpreter', 'latex');
set(gca, 'FontSize', font_size)

subplot(2, 4, 4);
pcolor(rhos0*y2d, rhos0*x2d, phi_m(:, :, zdiag));
colormap jet;  colorbar;  shading interp;
xlabel('y (cm)');
ylabel('x (cm)');
title(['$$\phi_', num2str(m), '\ \mathrm{\left(V\right)}$$'], 'interpreter', 'latex');
set(gca, 'FontSize', font_size)

subplot(2, 4, 5);
pcolor(zX, xX, squeeze(den_n(:, ydiag, :)));
colormap jet;  colorbar;  shading interp;
xlabel('z (cm)');
ylabel('x (cm)');
title('$$n\ \mathrm{\left(cm^{-3}\right)}$$', 'interpreter', 'latex');
set(gca, 'FontSize', font_size)

subplot(2, 4, 6);
pcolor(zX, xX, squeeze(den_m(:, ydiag, :)));
colormap jet;  colorbar;  shading interp;
xlabel('z (cm)');
ylabel('x (cm)');
title(['$$n_', num2str(m), '\ \mathrm{\left(cm^{-3}\right)}$$'], 'interpreter', 'latex');
set(gca, 'FontSize', font_size)

subplot(2, 4, 7);
pcolor(zX, xX, squeeze(phi_n(:, ydiag, :)));
colormap jet;  colorbar;  shading interp;
xlabel('z (cm)');
ylabel('x (cm)');
title('$$\phi\ \mathrm{\left(V\right)}$$', 'interpreter', 'latex');
set(gca, 'FontSize', font_size)

subplot(2, 4, 8);
pcolor(zX, xX, squeeze(phi_m(:, ydiag, :)));
colormap jet;  colorbar;  shading interp;
xlabel('z (cm)');
ylabel('x (cm)');
title(['$$\phi_', num2str(m), '\ \mathrm{\left(V\right)}$$'], 'interpreter', 'latex');
set(gca, 'FontSize', font_size)

print(gcf, '-dpng', ['test_m', num2str(m), 't', sprintf('%4.4d', idiag), '.png']);
