
clear all;  close all;

load('parameters.mat');
addpath(code_path);

%---input---
idiag = 71;
zdiag = 12;
l_max = 30;
m = 3;
fig_position = [50, 50, 1200, 900];
font_size = 15;
%-----------

build_grid_2d;
data_file = fullfile('data', ['dat', sprintf('%4.4d', idiag)]);
load(data_file);

den_before = denref * squeeze(den(:, :, zdiag));
phi_before = Tref * squeeze(phi(:, :, zdiag));

den_m3 = filter_m(den_before, m, l_max);
phi_m3 = filter_m(phi_before, m, l_max);

fig = figure();
set(fig, 'Position', fig_position);

subplot(2, 2, 1);
pcolor(rhos0*y2d, rhos0*x2d, den_before);
colormap jet;  colorbar;  shading interp;
xlabel('y (cm)');
ylabel('x (cm)');
title('$$n\ \mathrm{\left(cm^{-3}\right)}$$ before', 'interpreter', 'latex');
set(gca, 'FontSize', font_size)

subplot(2, 2, 2);
pcolor(rhos0*y2d, rhos0*x2d, den_m3);
colormap jet;  colorbar;  shading interp;
xlabel('y (cm)');
ylabel('x (cm)');
title('$$n_3\ \mathrm{\left(cm^{-3}\right)}$$', 'interpreter', 'latex');
set(gca, 'FontSize', font_size)

subplot(2, 2, 3);
pcolor(rhos0*y2d, rhos0*x2d, phi_before);
colormap jet;  colorbar;  shading interp;
xlabel('y (cm)');
ylabel('x (cm)');
title('$$\phi\ \mathrm{\left(V\right)}$$ before', 'interpreter', 'latex');
set(gca, 'FontSize', font_size)

subplot(2, 2, 4);
pcolor(rhos0*y2d, rhos0*x2d, phi_m3);
colormap jet;  colorbar;  shading interp;
xlabel('y (cm)');
ylabel('x (cm)');
title('$$\phi_3\ \mathrm{\left(V\right)}$$', 'interpreter', 'latex');
set(gca, 'FontSize', font_size)

print(gcf, '-dpng', 'test_m3.png');
