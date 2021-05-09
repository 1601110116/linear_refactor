
clear all;  close all;

load('parameters.mat');
addpath(code_path);

%---input---
idiag = 1212;
zdiag = 5;
font_size = 15;
fig_position = [50, 50, 1200, 500];
l_max = 30;
m_max = 20;
%-----------

build_grid_2d;
data_file = fullfile('data', ['dat', sprintf('%4.4d', idiag)]);
load(data_file);

phi_before = Tref * squeeze(phi(:, :, zdiag));
%phi_before = denref * squeeze(den(:, :, zdiag));
phi_after = zeros(size(phi_before));
for m = 0: 10
	phi_after = phi_after + filter_m(phi_before, m, l_max);
end

fig = figure();
set(fig, 'Position', fig_position);
sgtitle('$$\phi$$ (V)', 'interpreter', 'latex');

subplot(1, 2, 1);
pcolor(y2d*rhos0, x2d*rhos0, phi_before);
colormap jet;  colorbar;  shading interp;
xlabel('y (cm)');
ylabel('x (cm)');
title('before');
set(gca, 'FontSize', font_size);

subplot(1, 2, 2);
pcolor(y2d*rhos0, x2d*rhos0, phi_after);
colormap jet;  colorbar;  shading interp;
xlabel('y (cm)');
ylabel('x (cm)');
title('after');
set(gca, 'FontSize', font_size);


print(gcf, '-dpng', 'test_fb.png');
