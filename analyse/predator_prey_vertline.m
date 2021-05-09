% This script draws the zonal-averaged auto-spectrum versus m and r

clear;  close all;
global den Te pe vi jz ve phi vEx vEy vdex vdey dt inv_nustar calc ddx

load('parameters.mat');

addpath(code_path);
last_file = get_last_file('data');
last_diag = str2num(last_file(end-7: end-4));
save('parameters.mat', 'last_diagnose', '-append');

%%
%---input---
start_diag = 549;
end_diag = last_diag;
r_dst = 3;  % cm, used to separate shear flow growth phase from collapse phase by vertical lines
line_width = 1.2;  % width of the vertical lines
r_start = 0.5;  % cm
r_end = 9.5;  % cm
font_size = 20; 
fig_position = [50, 50, 1000, 1400]; % [left, bottom, width, height]
no_size = 20;  % font size of the number of each subplot
% for the area of autopower(not normalized) lower than 10^autopower_clim_low, we always, 
%  we always use the color for 10^autopower_clim_low
autopower_clim_low = 16;  
int_flux_clim_low = 15;  % 10^int_flux_clim_low
int_vE2_clim_low = 8;
x_tick = ceil(start_diag*dt*t0*nt_per_diagnose*1e3): 2: floor(end_diag*dt*t0*nt_per_diagnose*1e3);
y_tick = 1:2:9;
height_incre = 0.008;  % heigh increment of each axes with respect to default value by subplot()
width_incre = 0.08;
left_shift = 0.008;
dst_path = 'Documentation/KHv2/';  % specify where the .eps output should be copied to
%------------
%%

build_grid;
rX = rhos0 * r;
nr_dst = find(rX >= r_dst, 1, 'first');
nr_start = find(rX >= r_start, 1, 'first');
nr_end = find(rX <= r_end, 1, 'last');
nrdiag = nr_end - nr_start + 1;

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
m1_autopower_t = zeros(nrdiag, n_diag);
m3_autopower_t = zeros(nrdiag, n_diag);
ZF_t = zeros(nrdiag, n_diag);
int_den_t = zeros(nrdiag, n_diag);
int_flux_t = zeros(nrdiag, n_diag);
int_vE2_t = zeros(nrdiag, n_diag);  % kinetic energy

for idiag = start_diag: end_diag
	disp(['analysing: step ', num2str(idiag), ' of ', num2str(end_diag)]);
	data_file = fullfile('data', ['dat', sprintf('%4.4d', idiag)]);
	load(data_file);
	sdata();
	[vErp, vEthtp] = vec2pol(vEx, vEy);
	ZF = squeeze(zonal_average(vEthtp(nr_start: nr_end, :, :)));
	ZF_t(:, idiag-start_diag+1) = ZF;

	denp = field2pol(den);
	den_pert = zonal_pert(denp);
	for iz = 2: nz+1
		den_m = fft(den_pert(nr_start: nr_end, 1:end-1, iz), [], 2);
		I1 = 0.5 * (abs(den_m(:, 2)) / (ntht/2)).^2;
		I3 = 0.5 * (abs(den_m(:, 4)) / (ntht/2)).^2;
		m1_autopower_t(:, idiag-start_diag+1) = m1_autopower_t(:, idiag-start_diag+1) + I1;
		m3_autopower_t(:, idiag-start_diag+1) = m3_autopower_t(:, idiag-start_diag+1) + I3;
	end
	int_den = 2*pi * reshape(r, [], 1) .* zonal_average(denp);
	int_den_t(:, idiag-start_diag+1) = int_den(nr_start: nr_end);
	int_flux = 2*pi * reshape(r, [], 1) .* zonal_average(denp .* vErp);
	int_flux_t(:, idiag-start_diag+1) = int_flux(nr_start: nr_end);
	int_vE2 = 2*pi * reshape(r, [], 1) .* zonal_average(vErp.^2 + vEthtp.^2);
	int_vE2_t(:, idiag-start_diag+1) = int_vE2(nr_start: nr_end);
end
m1_autopower_t = m1_autopower_t ./ nz;
m3_autopower_t = m3_autopower_t ./ nz;
tX = t0 * dt * nt_per_diagnose * (start_diag: end_diag);
tX_ms = 1e3 * tX;
ZF_dst = ZF_t(nr_dst, :);
peak_idiag = find(diff(sign(diff(ZF_dst))) == -2) + 1;
trough_idiag = find(diff(sign(diff(ZF_dst))) == 2) + 1;
peak_tX_ms = tX_ms(peak_idiag);
trough_tX_ms = tX_ms(trough_idiag);

%%
close
fig = figure;
set(fig, 'Position', fig_position);


subplot(5,1,1);
init_position = get(gca, 'Position');
gca_left = init_position(1) - width_incre/2 - left_shift;
gca_width = init_position(3) + width_incre;
gca_height = init_position(4) + height_incre;
gca_position = [gca_left, init_position(2), gca_width, gca_height];
set(gca, 'Position', gca_position);
pcolor(tX_ms, rX(nr_start: nr_end), log10(cs0*ZF_t));  hold on;
colormap jet;  colorbar;  shading interp;
title('$$\lg\left<v_{E,\theta}/\left(\rm{cm/s}\right)\right>$$', 'interpreter', 'latex');
set(gca, 'xticklabel', []);
ylabel('$$r$$/cm', 'interpreter', 'latex');
set(gca, 'XTick', x_tick);
set(gca, 'YTick', y_tick);
x_lim = get(gca, 'XLim');
y_lim = get(gca, 'YLim');
text(x_lim(1), y_lim(end), ['\fontsize{', num2str(no_size), '}(a)'], ...
	'VerticalAlignment', 'top', 'color', 'white');
set(gca, 'FontSize', font_size);
for iline = 1: length(peak_idiag)
	plot([peak_tX_ms(iline), peak_tX_ms(iline)], y_lim, 'w-', 'LineWidth', line_width);  hold on;
end
for iline = 1: length(trough_idiag)
	plot([trough_tX_ms(iline), trough_tX_ms(iline)], y_lim, 'k-', 'LineWidth', line_width);  hold on;
end
hold off

m3_clim_high = max(max(max(log10(denref^2*m3_autopower_t), [], 3), [], 2), [], 1);
m1_clim_high = max(max(max(log10(denref^2*m1_autopower_t), [], 3), [], 2), [], 1);
autopower_clim = [autopower_clim_low, max(m3_clim_high, m1_clim_high)];

subplot(5,1,2);
init_position = get(gca, 'Position');
gca_position = [gca_left, init_position(2), gca_width, gca_height];
set(gca, 'Position', gca_position);
pcolor(tX_ms, rX(nr_start: nr_end), log10(denref^2*m3_autopower_t));  hold on;
colormap jet;  colorbar;  shading interp;
set(gca, 'CLim', autopower_clim);
title('$$\lg\left<\tilde{n}_{3}^{2}/\rm{cm^{-3}}\right>$$', 'interpreter', 'latex');
set(gca, 'xticklabel', []);
ylabel('$$r$$/cm', 'interpreter', 'latex');
set(gca, 'XTick', x_tick);
set(gca, 'YTick', y_tick);
x_lim = get(gca, 'XLim');
y_lim = get(gca, 'YLim');
text(x_lim(1), y_lim(end), ['\fontsize{', num2str(no_size), '}(b)'], ...
	'VerticalAlignment', 'top', 'color', 'white');
set(gca, 'FontSize', font_size);
for iline = 1: length(peak_idiag)
	plot([peak_tX_ms(iline), peak_tX_ms(iline)], y_lim, 'w-', 'LineWidth', line_width);  hold on;
end
for iline = 1: length(trough_idiag)
	plot([trough_tX_ms(iline), trough_tX_ms(iline)], y_lim, 'k-', 'LineWidth', line_width);  hold on;
end
hold off

subplot(5,1,3);
init_position = get(gca, 'Position');
gca_position = [gca_left, init_position(2), gca_width, gca_height];
set(gca, 'Position', gca_position);
%pcolor(tX_ms, rX(nr_start: nr_end), log10(denref^2*m1_autopower_t.*...
%	repmat(rX(nr_start:nr_end)'.^2, 1, n_diag)));
pcolor(tX_ms, rX(nr_start: nr_end), log10(denref^2*m1_autopower_t));  hold on;
colormap jet;  colorbar;  shading interp;
set(gca, 'CLim', autopower_clim);
title('$$\lg\left<\tilde{n}_{1}^{2}/\rm{cm^{-3}}\right>$$', 'interpreter', 'latex');
set(gca, 'xticklabel', []);
ylabel('$$r$$/cm', 'interpreter', 'latex');
set(gca, 'XTick', x_tick);
set(gca, 'YTick', y_tick);
x_lim = get(gca, 'XLim');
y_lim = get(gca, 'YLim');
text(x_lim(1), y_lim(end), ['\fontsize{', num2str(no_size), '}(c)'], ...
	'VerticalAlignment', 'top', 'color', 'white');
set(gca, 'FontSize', font_size);
for iline = 1: length(peak_idiag)
	plot([peak_tX_ms(iline), peak_tX_ms(iline)], y_lim, 'w-', 'LineWidth', line_width);  hold on;
end
for iline = 1: length(trough_idiag)
	plot([trough_tX_ms(iline), trough_tX_ms(iline)], y_lim, 'k-', 'LineWidth', line_width);  hold on;
end
hold off

%subplot(6,1,5);
%int_vE2_tX = rhos0*cs0^2*int_vE2_t;
%int_vE2_tX(int_vE2_tX < (10^int_vE2_clim_low)) = 10^int_vE2_clim_low;
%pcolor(tX_ms, rX(nr_start: nr_end), log10(int_vE2_tX));
%colormap jet;  colorbar;  shading interp;
%title(['$$\lg\left<2\pi r\left(v_{E,r}^{2}+v_{E,\theta}^{2}\right)', ...
%	'/\left(\rm{cm^{2}/s^{2}\cdot cm}\right)\right>$$'], ...
%	'interpreter', 'latex');
%xlabel('$$t$$/ms', 'interpreter', 'latex');
%ylabel('$$r$$/cm', 'interpreter', 'latex');
%set(gca, 'XTick', x_tick);
%set(gca, 'YTick', y_tick);
%x_lim = get(gca, 'XLim');
%y_lim = get(gca, 'YLim');
%text(x_lim(1), y_lim(end), ['\fontsize{', num2str(no_size), '}(f)'], ...
%	'VerticalAlignment', 'top', 'color', 'white');
%set(gca, 'FontSize', font_size);
%

subplot(5,1,4)
init_position = get(gca, 'Position');
gca_position = [gca_left, init_position(2), gca_width, gca_height];
set(gca, 'Position', gca_position);
pcolor(tX_ms, rX(nr_start: nr_end), rhos0*denref*int_den_t);  hold on;
colormap jet;  colorbar;  shading interp;
title('$$\left<2\pi r\cdot n/\rm{cm^{-3}}\right>$$', 'interpreter', 'latex');
set(gca, 'xticklabel', []);
ylabel('$$r$$/cm', 'interpreter', 'latex')
set(gca, 'XTick', x_tick);
set(gca, 'YTick', y_tick);
x_lim = get(gca, 'XLim');
y_lim = get(gca, 'YLim');
text(x_lim(1), y_lim(end), ['\fontsize{', num2str(no_size), '}(d)'], ...
	'VerticalAlignment', 'top', 'color', 'white');
set(gca, 'FontSize', font_size);
for iline = 1: length(peak_idiag)
	plot([peak_tX_ms(iline), peak_tX_ms(iline)], y_lim, 'w-', 'LineWidth', line_width);  hold on;
end
for iline = 1: length(trough_idiag)
	plot([trough_tX_ms(iline), trough_tX_ms(iline)], y_lim, 'k-', 'LineWidth', line_width);  hold on;
end
hold off


subplot(5,1,5)
init_position = get(gca, 'Position');
gca_position = [gca_left, init_position(2), gca_width, gca_height];
set(gca, 'Position', gca_position);
positive_int_flux_tX = rhos0*denref*cs0*int_flux_t;
positive_int_flux_tX(positive_int_flux_tX <= 0) = nan;
lg_int_flux_tX = log10(positive_int_flux_tX);
pcolor(tX_ms, rX(nr_start: nr_end), lg_int_flux_tX);  hold on;
colormap jet;  colorbar;  shading interp;
tmp_clim = get(gca, 'CLim');
set(gca, 'CLim', [int_flux_clim_low, tmp_clim(end)]);
title('$$\lg\left<2\pi r\cdot nv_{E,r}/\left(\rm{cm^{-3}\cdot cm/s\cdot cm}\right)\right>$$', ...
	'interpreter', 'latex')
%set(gca, 'xticklabel', []);
xlabel('$$t$$/ms', 'interpreter', 'latex');
ylabel('$$r$$/cm', 'interpreter', 'latex')
set(gca, 'XTick', x_tick);
set(gca, 'YTick', y_tick);
x_lim = get(gca, 'XLim');
y_lim = get(gca, 'YLim');
text(x_lim(1), y_lim(end), ['\fontsize{', num2str(no_size), '}(e)'], ...
	'VerticalAlignment', 'top', 'color', 'black');
set(gca, 'FontSize', font_size);
for iline = 1: length(peak_idiag)
	plot([peak_tX_ms(iline), peak_tX_ms(iline)], y_lim, 'w-', 'LineWidth', line_width);  hold on;
end
for iline = 1: length(trough_idiag)
	plot([trough_tX_ms(iline), trough_tX_ms(iline)], y_lim, 'k-', 'LineWidth', line_width);  hold on;
end
hold off

%subplot(6,1,6,margins)
%negative_int_flux_tX = rhos0*denref*cs0*int_flux_t;
%negative_int_flux_tX(negative_int_flux_tX >= 0) = nan;
%pcolor(tX_ms, rX(nr_start: nr_end), negative_int_flux_tX);
%colormap jet;  colorbar;  shading interp;
%title('$$2\pi r\left<nv_{E,r}\right>/\left(\rm{cm^{-3}\cdot cm/s\cdot cm}\right)$$', ...
%	'interpreter', 'latex')
%xlabel('$$t$$/ms', 'interpreter', 'latex');
%ylabel('$$r$$/cm', 'interpreter', 'latex');
%set(gca, 'XTick', x_tick);
%set(gca, 'YTick', y_tick);
%x_lim = get(gca, 'XLim');
%y_lim = get(gca, 'YLim');
%text(x_lim(1), y_lim(end), ['\fontsize{', num2str(no_size), '}(f)'], ...
%	'VerticalAlignment', 'top', 'color', 'black');
%set(gca, 'FontSize', font_size);
%
pic_name = fullfile('figs', 'ppv2');
print(fig, '-dpng', [pic_name, '.png']);
print(fig, '-depsc', [pic_name, '.eps']);
if exist(dst_path, 'dir')
	copyfile([pic_name, '.eps'], dst_path);
else
	disp('Destination path do not exist');
end
