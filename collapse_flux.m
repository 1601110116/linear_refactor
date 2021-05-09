% This script draws the particle flux during two typical collapses

clear;  close all;
global den Te pe vi jz ve phi vEx vEy vdex vdey dt inv_nustar calc ddx

%%
load('parameters.mat');
addpath(code_path);

%---input---
start_diag1 = 1135;
start_diag2 = 1437; %1241;
end_diag1 = 1157;
end_diag2 = 1459; %1263;
r_diag = 1.5;  % cm
fig_position = [50, 50, 800, 500];  % [left, bottom, width, height];
y_tick = 1e16*(-1:6);
y_lim = 1e16*[-1.5, 6.5];
x_tick1 = 35.5:0.2:36.3;
x_tick2 = 38.8:0.2:39.6;
line_width = 2;
marker_size = 10;
font_size = 20;
no_size = 25;
dst_path = 'Documentation/KH/';  % specify where the .eps output should be copied to
%-----------
%%

build_grid;
rX = rhos0 * r;
nr_diag = find(rX >= r_diag, 1, 'first');
jz = zeros(nx+2, nx+2, nz+2);
ve = jz;  vEx = jz;  vEy = jz;  vdex = jz;  vdey = jz; 
generate_constants(height, radius, dx, dz, nx, nz, dt, ... 
    rdif, dif_perp_in, dif_perp_out, dif_z_in, dif_z_out, ... 
	rconduct, conduct_z_in, conduct_z_out, viscosity, ... 
    den_damp, momentum_damp, Te_tanhsrc_max, Te_tanhsrc_radius, ... 
    Te_tanhsrc_incline, Te_gausssrc_magnitude, Te_gausssrc_sigma, ... 
    den_tanhsrc_max, den_tanhsrc_radius, den_tanhsrc_incline, ... 
    den_gausssrc_magnitude, den_gausssrc_sigma);
ndiag = end_diag1 - start_diag1 + 1;
flux_t1 = zeros(1, ndiag);
flux_t2 = zeros(1, ndiag);

for idiag = start_diag1: end_diag1
	data_file = fullfile('data', ['dat', sprintf('%4.4d', idiag)]);
	load(data_file);
	sdata();
	[vErp, vEthtp] = vec2pol(vEx, vEy);
	denp = field2pol(den);
	flux_t1(idiag-start_diag1+1) = squeeze(zonal_average(denp(nr_diag, :, :) .* vErp(nr_diag, :, :)));
end
for idiag = start_diag2: end_diag2
	data_file = fullfile('data', ['dat', sprintf('%4.4d', idiag)]);
	load(data_file);
	sdata();
	[vErp, vEthtp] = vec2pol(vEx, vEy);
	denp = field2pol(den);
	flux_t2(idiag-start_diag2+1) = squeeze(zonal_average(denp(nr_diag, :, :) .* vErp(nr_diag, :, :)));
end
%%
close all;
tX1_ms = 1e3 * t0 * dt * nt_per_diagnose * (start_diag1: end_diag1);
tX2_ms = 1e3 * t0 * dt * nt_per_diagnose * (start_diag2: end_diag2);
fig = figure;
set(fig, 'Position', fig_position);

subplot(1,2,1);
plot(tX1_ms, denref*cs0*flux_t1, 'b+', 'LineWidth', line_width, 'MarkerSize', marker_size);
xlabel('$$t$$/ms', 'interpreter', 'latex');
ylabel('$$\left<nv_{E,r}\right>/\rm\left(cm^{-3}\cdot cm/s\right)$$', 'interpreter', 'latex')
title('With','FontWeight', 'normal');
x_lim1 = [x_tick1(1), x_tick1(end)];
set(gca, 'XTick', x_tick1);
set(gca, 'XLim', x_lim1);
set(gca, 'YTick', y_tick);
set(gca, 'YLim', y_lim);
text(x_lim1(1), y_lim(end), ['\fontsize{', num2str(no_size), '}(a)'], ...
	'VerticalAlignment', 'top', 'color', 'black');
set(gca, 'FontSize', font_size);

subplot(1,2,2);
plot(tX2_ms, denref*cs0*flux_t2, 'b+', 'LineWidth', line_width, 'MarkerSize', marker_size);
title('W/O', 'FontWeight', 'normal');
xlabel('$$t$$/ms', 'interpreter', 'latex');
ylabel('$$\left<nv_{E,r}\right>/\rm\left(cm^{-3}\cdot cm/s\right)$$', 'interpreter', 'latex')
x_lim2 = [x_tick2(1), x_tick2(end)];
%set(gca, 'XTick', x_tick2);
%set(gca, 'XLim', x_lim2);
set(gca, 'YTick', y_tick);
set(gca, 'YLim', y_lim);
text(x_lim2(1), y_lim(end), ['\fontsize{', num2str(no_size), '}(b)'], ...
	'VerticalAlignment', 'top', 'color', 'black');
set(gca, 'FontSize', font_size);

print(gcf, '-dpng', 'collapse_flux.png');
print(gcf, '-depsc', 'collapse_flux.eps');
if exist(dst_path, 'dir')
	copyfile('collapse_flux.eps', dst_path);
else
	disp('Destination path do not exist');
end
