% This script is a 2D version of vrvz_movie.m, showing the evolution of vrvz and vi

clear;  close all;
global vx vy phi

load('parameters.mat');
addpath(code_path);

%%
%---input---
start_diag = 100;%1216; %1020;
end_diag = 110;%1236; %1037;
% for the vr-vz 2D plot
rdiag = 3;
% for the vi profile
r_start = 0.2;
r_end = 7.5;
zdiag = ceil(length(zX) / 2);  % middle x-y plane
font_size = 20;
fig_position = [50, 50, 700, 1000];
line_width = 1;
%-----------

build_grid;
rX = rhos0 * r;
irdiag = min(find(rX > rdiag));
nr_start = find(rX >= r_start, 1, 'first');
nr_end = find(rX <= r_end, 1, 'last');
nrdiag = nr_end - nr_start + 1;
n_diag = end_diag - start_diag + 1;
%%

vx = zeros(nx+2, nx+2, nz+2);
vy = zeros(nx+2, nx+2, nz+2);
vr_pert_t = zeros(ntht, n_diag);
vi_pert_t = vr_pert_t;
vi_mean_t = zeros(nrdiag, n_diag);

for idiag = start_diag: end_diag
	disp(['step ', num2str(idiag), ' of ', num2str(end_diag)]);
	data_file = fullfile('data', ['dat', sprintf('%4.4d', idiag)]);
	load(data_file, 'phi', 'vi');
	vx(2:end-1, 2:end-1, :) = 1/(2*dx) * ...
		(phi(2:end-1, 1:end-2, :) - phi(2:end-1, 3:end, :));
	vy(2:end-1, 2:end-1, :) = 1/(2*dx) * ...
        (phi(3:end, 2:end-1, :) - phi(1:end-2, 2:end-1, :));
	[vr, vtht] = vec2pol(vx, vy);
	vip = field2pol(vi);
	vr_pert = get_pert(vr);
	vi_pert = get_pert(vip);
	vi_mean = get_average(vip);
	vr_pert_t(:, idiag - start_diag + 1) = squeeze(vr_pert(irdiag, :, zdiag));
	vi_pert_t(:, idiag - start_diag + 1) = squeeze(vi_pert(irdiag, :, zdiag));
	vi_mean_t(:, idiag - start_diag + 1) = squeeze(vi_mean(nr_start: nr_end, 1, zdiag));
end

%%
close;
dtX = 1e3 * t0 * dt * nt_per_diagnose;
tX = dtX * (start_diag: end_diag);
fig = figure('name', 'vi growth');
r_axis = rX(nr_start: nr_end);
line_colors = jet(n_diag);

subplot('Position', [0.1800 0.5700 0.5500 0.3700]);
for iplot = 1:n_diag
	plot(r_axis, cs0 * vi_mean_t(:, iplot), 'Color', line_colors(iplot, :), ...
		'lineWidth', line_width);  hold on;
end
hold off;
grid on;
xlabel('$$r$$ (cm)', 'interpreter', 'latex')
ylabel('$$v_{\parallel i}\ \rm \left(cm/s\right)$$', 'interpreter', 'latex');
x_lim = get(gca, 'XLim');
y_lim = get(gca, 'YLim');
text(x_lim(1), y_lim(end), ['\fontsize{', num2str(font_size), '}(a)'], ...
	'VerticalAlignment', 'top', 'color', 'black');
set(gca, 'YLim', y_lim);
set(gca, 'fontSize', font_size);

subplot('Position', [0.1800 0.1100, 0.5500, 0.3700])
for iplot = 1: n_diag
	plot(cs0 * vi_pert_t(:, iplot), cs0 * vr_pert_t(:, iplot), '+', 'Color', ...
		line_colors(iplot, :), 'lineWidth', line_width);  hold on;
end
hold off;
grid on;
xlabel('$$\tilde{v}_{\parallel i}\ \rm \left(cm/s\right)$$', 'interpreter', 'latex');
ylabel('$$\tilde{v}_{E,r}\ \rm \left(cm/s\right)$$', 'interpreter', 'latex');
x_lim = get(gca, 'XLim');
y_lim = get(gca, 'YLim');
text(x_lim(1), y_lim(end), ['\fontsize{', num2str(font_size), '}(b)'], ... 
    'VerticalAlignment', 'top', 'color', 'black');
set(gca, 'YLim', y_lim);
set(gca, 'XLim', x_lim);
set(gca, 'fontSize', font_size);

subplot('Position', [0.7500 0.1100, 0.2000, 0.8000]);
set(gca, 'Visible', 'off');
colormap(line_colors);  cbar = colorbar;
set(gca, 'CLim', [tX(1) - 0.5*dtX, tX(end) + 0.5*dtX])
set(cbar.Label, 'String', '$$t$$ (ms)');
set(cbar.Label, 'interpreter', 'latex');
set(cbar.Label, 'fontSize', font_size);
set(cbar, 'Location', 'west');
cbar.Ticks = tX(1:2:end);
set(gca, 'fontSize', font_size);

set(fig, 'Position', fig_position);
pic_name = fullfile('figs', ['t=', sprintf('%4.4d', start_diag), '-', ...
	sprintf('%4.4d', end_diag), '_r=', num2str(rX(irdiag)), '_vrvz_2d']);
print(gcf, '-dpng', [pic_name, '.png']);
print(gcf, '-depsc', [pic_name, '.eps']);


