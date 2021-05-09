% This script draws the theta-direction EXB and electron diamagnetic velocities 
%  versus time
% Copy this script to the folder of the data and then run it

clear;  close all;
global den Te pe vi jz ve phi vEx vEy vdex vdey dt inv_nustar calc ddx

%%
load('parameters.mat');

addpath(code_path)
last_file = get_last_file('data');
last_diag = str2num(last_file(end-7: end-4));


%---input---
start_diag = 549;
end_diag = last_diag;
r_diag = 3.6;  % vEtht at r=r_diag cm is diagnosed
r_start = 0.5;  % the maximum value is taken from [r_start, r_end] cm
r_end = 9.5;
font_size = 13;
fig_position = [50, 50, 700, 900];  % [left, bottom, width, height]
line_width = 1;
no_size = 18;  % font size of the number of each subplot
legend_font_size = 13;
x_tick = 0: 5: floor(end_diag*dt*t0*nt_per_diagnose*1e3);
flow_y_tick = 0:4e4:20e4;
vorticity_y_tick = 0:0.5e5:4e5;
dst_path = 'Documentation/KH/';  % specify where the .eps output should be copied to
%-----------
%%

build_grid;
rX = rhos0*r;
nr_diag = find(rX >= r_diag, 1, 'first');
nr_start = find(rX >= r_start, 1, 'first');
if nr_start <= 2
	errpr(['r_start must be greater than ', num2str(2*dr*rhos0), ...
		' for the current grid']);
end
nr_end = find (rX <= r_end, 1, 'last');
jz = zeros(nx+2, nx+2, nz+2);
ve = jz;  vEx = jz;  vEy = jz;  vdex = jz;  vdey = jz;
generate_constants(height, radius, dx, dz, nx, nz, dt, ... 
    rdif, dif_perp_in, dif_perp_out, dif_z_in, dif_z_out, ... 
    rconduct, conduct_z_in, conduct_z_out, viscosity, ... 
    den_damp, momentum_damp, Te_tanhsrc_max, Te_tanhsrc_radius, ... 
    Te_tanhsrc_incline, Te_gausssrc_magnitude, Te_gausssrc_sigma, ... 
    den_tanhsrc_max, den_tanhsrc_radius, den_tanhsrc_incline, ... 
    den_gausssrc_magnitude, den_gausssrc_sigma);
n_diag = end_diag - start_diag + 1;

ZF_rdiag_t = zeros(1, n_diag);
ZF_max_t = ZF_rdiag_t;
ZF_rshear_max_t = ZF_rdiag_t;
ZF_vorticity_max_t = ZF_rdiag_t;
vEtht_rshear_max_t = ZF_rdiag_t;
vEtht_vorticity_max_t = ZF_rdiag_t;
w_max_t = ZF_rdiag_t;


for idiag = start_diag: end_diag
	disp(['getting zonal averaged velocities: step ', num2str(idiag), ' of ', ...
		num2str(end_diag)]);
	data_file = fullfile('data', ['dat', sprintf('%4.4d', idiag)]);
	load(data_file);
	sdata();

	[vErp, vEthtp] = vec2pol(vEx, vEy);

	ZF_rdiag = squeeze(zonal_average(vEthtp(nr_diag, :, :)));
	ZF_rdiag_t(idiag-start_diag+1) = ZF_rdiag;

	ZF = squeeze(zonal_average(vEthtp(nr_start-1: nr_end+1, :, :)));
	ZF_max_t(idiag-start_diag+1) = max(ZF(2:end-1));

	ZF_rshear = 1/(2*dr) .* (ZF(3:end) - ZF(1:end-2));
	ZF_rshear_max_t(idiag-start_diag+1) = max(ZF_rshear);
	r_tmp = reshape(r(nr_start-1: nr_end+1), [], 1);
	rZF = r_tmp .* ZF;
	ZF_vorticity = 1/(2*dr)./r_tmp(2: end-1) .* (rZF(3:end) - rZF(1:end-2));
	ZF_vorticity_max_t(idiag-start_diag+1) = max(ZF_vorticity);

	vEtht_rshear = 1/(2*dr) .* (vEthtp(nr_start+1: nr_end+1, :, :) ...
		- vEthtp(nr_start-1: nr_end-1, : ,:));
	vEtht_rshear_max_t(idiag-start_diag+1) = squeeze( ...
		max(max(max(vEtht_rshear, [], 3), [], 2), [], 1));
	r_tmp = repmat(r_tmp, 1, ntht, nz+2);
	rvEtht = r_tmp .* vEthtp(nr_start-1: nr_end+1, :, :);
	vEtht_vorticity = 1/(2*dr)./r_tmp(2: end-1, :, :) ...
		.* (rvEtht(3:end, :, :) - rvEtht(1:end-2, : ,:));
	vEtht_vorticity_max_t(idiag-start_diag+1) = squeeze( ...
		max(max(max(vEtht_vorticity, [], 3), [], 2), [], 1));


	wp = field2pol(w);
	w_max = squeeze(max(max(max(wp(nr_start: nr_end,:,:), [], 3), [], 2), [], 1));
	w_max_t(idiag-start_diag+1) = w_max;
end
%%
close;
fig = figure;
tX_ms = 1e3 * t0 * dt * nt_per_diagnose * (start_diag: end_diag);

subplot(3,1,1);
plot(tX_ms, cs0*ZF_rdiag_t, 'b-', 'LineWidth', line_width);  hold on;
plot(tX_ms, cs0*ZF_max_t, 'r-', 'LineWidth', line_width);  hold off;
ylabel('Flow $$\left(\rm{cm/s}\right)$$', 'interpreter', 'latex');
%title('Flow', 'interpreter', 'latex');
set(gca, 'XTick', x_tick);
x_lim = get(gca, 'XLim');
y_tick = flow_y_tick;set(gca, 'XLim', x_lim);
y_lim = [y_tick(1), y_tick(end)];
set(gca, 'YTick', y_tick);
set(gca, 'YLim', y_lim);
%text(x_lim(end), y_lim(end), ['\fontsize{', num2str(no_size), '}r=3.6 cm'], ...
%	'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'color', 'black');
text(x_lim(1), y_lim(1), ['\fontsize{', num2str(no_size), '}(a)'], ...
	'VerticalAlignment', 'bottom', 'color', 'black');
set(gca, 'fontSize', font_size);
lgd1 = legend('$$\left<v_{E,\theta}\right>\Big|_{r=3.6\rm{cm}}$$', ...
	'$$\left\{\left<v_{E,\theta}\right>\right\}_{max}$$');
set(lgd1, 'interpreter', 'latex');
set(lgd1, 'FontSize', legend_font_size);
set(lgd1, 'Location', 'NorthEast');
set(lgd1, 'NumColumns', 2);
set(lgd1, 'Box', 'off');

subplot(3,1,2);
plot(tX_ms, cs0/rhos0*ZF_rshear_max_t, 'b-', 'LineWidth', line_width);  hold on;set(gca, 'XLim', x_lim);
plot(tX_ms, cs0/rhos0*w_max_t, 'g-', 'LineWidth', line_width);  hold on;
plot(tX_ms, cs0/rhos0*vEtht_rshear_max_t, 'r-', 'LineWidth', line_width);  hold off;
ylabel('$$w$$ components $$\left(\rm{s^{-1}}\right)$$', 'interpreter', 'latex');
%title('Vorticity from $$\hat{\bf r\it}\times\partial_{r}\bf v\it _{E}$$', 'interpreter', 'latex');
set(gca, 'XTick', x_tick);
set(gca, 'XLim', x_lim);
y_tick = vorticity_y_tick;
y_lim = [y_tick(1), y_tick(end)];
set(gca, 'YTick', y_tick);
set(gca, 'YLim', y_lim);
text(x_lim(1), y_lim(1), ['\fontsize{', num2str(no_size), '}(b)'], ...
	'VerticalAlignment', 'bottom', 'color', 'black');
set(gca, 'fontSize', font_size);
lgd2 = legend('$$\rm Max \it \left\{\partial_{r}\left<v_{E,\theta}\right>\right\}$$', ...
	'$$\rm Max \it \left\{w\right\}$$', ...
    '$$\rm Max \it \left\{\partial_{r}v_{E,\theta}\right\}$$');
set(lgd2, 'interpreter', 'latex');
set(lgd2, 'FontSize', legend_font_size);
set(lgd2, 'Location', 'NorthWest');
set(lgd2, 'NumColumns', 2);
set(lgd2, 'Box', 'off');

subplot(3,1,3);
plot(tX_ms, cs0/rhos0*ZF_vorticity_max_t, 'b-', 'LineWidth', line_width);  hold on;
plot(tX_ms, cs0/rhos0*w_max_t, 'g-', 'LineWidth', line_width);
plot(tX_ms, cs0/rhos0*vEtht_vorticity_max_t, 'r-', 'LineWidth', line_width); hold off;
xlabel('$$t$$/ms', 'interpreter', 'latex');
ylabel('$$w$$ components $$\left(\rm{s^{-1}}\right)$$', 'interpreter', 'latex');
%title('Vorticity from $$\nabla\times\left(v_{E,\theta}\hat{\bf \theta\it}\right)$$', ...
%	'interpreter', 'latex');
set(gca, 'XTick', x_tick);
set(gca, 'XLim', x_lim);
y_tick = vorticity_y_tick;
y_lim = [y_tick(1), y_tick(end)];
set(gca, 'YTick', y_tick);
set(gca, 'YLim', y_lim);
text(x_lim(1), y_lim(1), ['\fontsize{', num2str(no_size), '}(c)'], ...
	'VerticalAlignment', 'bottom', 'color', 'black');
set(gca, 'fontSize', font_size);
lgd3 = legend('$$\rm Max\left\{\frac{1}{\it r}\it\partial_{r}\left<rv_{E,\theta}\right>\right\}$$', ...
	'$$\rm Max \it \left\{w\right\}$$', ...
	'$$\rm Max \left\{\frac{1}{\it r}\it\partial_{r}\left(rv_{E,\theta}\right)\right\}$$');
set(lgd3, 'interpreter', 'latex');
set(lgd3, 'FontSize', legend_font_size);
set(lgd3, 'Location', 'NorthWest');
set(lgd3, 'NumColumns', 2);
set(lgd3, 'Box', 'off');

set(fig, 'Position', fig_position);
png_file = fullfile('figs', 'w_contributions_short.png');
eps_file = fullfile('figs', 'w_contributions_short.eps');
print(gcf, '-dpng', png_file);
print(gcf, '-depsc', eps_file);
if exist(dst_path, 'dir')
	copyfile(eps_file, dst_path);
else
	disp('Destination path do not exist');
end
