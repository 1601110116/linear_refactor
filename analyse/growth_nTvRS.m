% This scpript plots the growth of ZF in a shear flow growth phase

clear;  close all;
global phi vx vy

%%
load('parameters.mat');
addpath(code_path);

%---input---
start_diag = 1226; %1035;
end_diag = 1246;  %1045;
r_start = 0.4;  % cm
r_end = 7.6;  % cm
font_size = 15;
legend_font_size = 13;
fig_position = [50, 50, 1000, 600];
cbar_position = [0.8305, 0.1371,0.03,0.7881];
line_width = 2;
%-----------

build_grid;
rX = rhos0 * r;
nr_start = find(rX >= r_start, 1, 'first');
if nr_start < 2
    error(['r_start must be greater than ', num2str(dr*rhos0), ...
        ' for the current grid']);
end
nr_end = find(rX <= r_end, 1, 'last');
nrdiag = nr_end - nr_start + 1;
n_diag = end_diag - start_diag + 1;
den_t = zeros(nrdiag, n_diag);
Te_t = zeros(nrdiag, n_diag);
zf_t = zeros(nrdiag, n_diag);
RS_t = zeros(nrdiag, n_diag);

%%
vx = zeros(nx+2, nx+2, nz+2);
vy = zeros(nx+2, nx+2, nz+2);
for idiag = start_diag: end_diag
	disp(['getting flow, step ', num2str(idiag), ' of ', ...
		num2str(end_diag)]);
	data_file = fullfile('data', ['dat', sprintf('%4.4d', idiag)]);
	load(data_file, 'phi', 'den', 'Te');

	denp = field2pol(den);
	den_t(:, idiag - start_diag + 1) = ...
		squeeze(zonal_average(denp(nr_start: nr_end, :, :)));
	Tep = field2pol(Te);
	Te_t(:, idiag - start_diag + 1) = ...
		squeeze(zonal_average(Tep(nr_start: nr_end, :, :)));
	vx(2:end-1, 2:end-1, :) = 1/(2*dx) * ...
		(phi(2:end-1, 1:end-2, :) - phi(2:end-1, 3:end, :));
	vy(2:end-1, 2:end-1, :) = 1/(2*dx) * ...
		(phi(3:end, 2:end-1, :) - phi(1:end-2, 2:end-1, :));
	[vr, vtht] = vec2pol(vx, vy);
	zf_t(:, idiag - start_diag + 1) = ...
		squeeze(zonal_average(vtht(nr_start: nr_end, :, :)));
	RS = squeeze(zonal_average(vtht(nr_start: nr_end, :, :) ...
		.* vr(nr_start: nr_end, :, :)));
	RS_t(:, idiag - start_diag + 1) = RS;
end

%%
close;
dtX = 1e3 * t0 * dt * nt_per_diagnose;
tX = dtX * (start_diag: end_diag);
fig = figure('name', 'flow growth');
r_axis = rX(nr_start: nr_end);
line_colors = jet(n_diag);

subplot('Position', [0.0650 0.5700 0.3000 0.3850]);
for iplot = 1: n_diag
	plot(r_axis, denref * den_t(:, iplot), 'Color', line_colors(iplot, :), ...
		'lineWidth', line_width);  hold on;
end
hold off;
grid on;
ylabel('$$\left<n\right>\ \rm\left(cm^{-3}\right)$$', 'interpreter', 'latex');
x_lim = get(gca, 'XLim');
y_lim = [0, 12e12];
text(x_lim(1), y_lim(end), ['\fontsize{', num2str(font_size), '}(a)'], ...
	'VerticalAlignment', 'top', 'color', 'black');
set(gca, 'YLim', y_lim)
set(gca, 'xticklabel', [])
set(gca, 'fontSize', font_size);

subplot('Position', [0.4500 0.5700 0.3000, 0.3850]);
for iplot = 1: n_diag
	plot(r_axis, Tref * Te_t(:, iplot), 'Color', line_colors(iplot, :), ...
		'lineWidth', line_width);  hold on;
end
hold off;
grid on;
ylabel('$$\left<T_e\right>\ \rm\left(eV\right)$$', 'interpreter', 'latex');
x_lim = get(gca, 'XLim');
y_lim = [0, 3.5];
text(x_lim(1), y_lim(end), ['\fontsize{', num2str(font_size), '}(b)'], ...
	'VerticalAlignment', 'top', 'color', 'black');
set(gca, 'YLim', y_lim)
set(gca, 'xticklabel', [])
set(gca, 'fontSize', font_size);

subplot('Position', [0.0650 0.1100 0.3000 0.3850]);
for iplot = 1: n_diag
	plot(r_axis, cs0 * zf_t(:, iplot), 'Color', line_colors(iplot, :), ...
		'lineWidth', line_width);  hold on;
end
hold off;
grid on;
xlabel('$$r$$ (cm)', 'interpreter', 'latex')
ylabel('$$\left<v_{E,\theta}\right>$$ (cm/s)', 'interpreter', 'latex');
x_lim = get(gca, 'XLim');
y_lim = get(gca, 'YLim');
text(x_lim(1), y_lim(end), ['\fontsize{', num2str(font_size), '}(c)'], ...
	'VerticalAlignment', 'top', 'color', 'black');
set(gca, 'YLim', y_lim)
set(gca, 'fontSize', font_size);

subplot('Position', [0.4500 0.1100 0.3000 0.3850]);
for iplot = 1: n_diag
	plot(r_axis, cs0^2 * RS_t(:, iplot), 'Color', line_colors(iplot, :), ...
		'lineWidth', line_width);  hold on;
end
hold off;
grid on;
xlabel('$$r$$ (cm)', 'interpreter', 'latex')
ylabel('$$\left<v_{E,r}v_{E,\theta}\right>\ \rm\left(cm^2/s^2\right)$$', ...
	'interpreter', 'latex')
x_lim = get(gca, 'XLim');
y_lim = get(gca, 'YLim');
text(x_lim(1), y_lim(end), ['\fontsize{', num2str(font_size), '}(d)'], ...
	'VerticalAlignment', 'top', 'color', 'black');
set(gca, 'YLim', y_lim)
set(gca, 'fontSize', font_size);

subplot('Position', [0.79, 0.1, 0.18, 0.85]);
set(gca, 'Visible', 'off');
colormap(line_colors);  cbar = colorbar;
set(gca, 'CLim', [tX(1) - 0.5*dtX, tX(end) + 0.5*dtX])
set(cbar.Label, 'String', '$$t$$ (ms)');
set(cbar.Label, 'interpreter', 'latex');
set(cbar.Label, 'fontSize', font_size);
set(cbar, 'Location', 'west');
%set(cbar.Ticks, tX);
cbar.Ticks = tX(1:2:end);
set(gca, 'fontSize', font_size);
set(fig, 'Position', fig_position);
pic_name = fullfile('figs', ['t=', sprintf('%4.4d', start_diag), '-', ...
	sprintf('%4.4d', end_diag), '_growth_nTvRS']);
print(gcf, '-dpng', [pic_name, '.png']);
print(gcf, '-depsc', [pic_name, '.eps']);
