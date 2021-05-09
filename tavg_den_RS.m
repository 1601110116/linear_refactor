% This script draws the azimuthal Renolds' stress using plasma potential and floating potential

clear;  close all;
global Te phi

load ('parameters.mat');
addpath(code_path);
last_file = get_last_file('data');
last_diag = str2num(last_file(end-7: end-4));
%%
%---input---
start_diag = 549;
end_diag = last_diag;
r_start = 0.4;
r_end = 9;
Lambda = 4.68;  % phi_fl = phi - Lambda*Te

fig_position = [50, 50, 1500, 700];  % [left, bottom, width, height]
x_tick = 0:6;
dst_path = 'Documentation/ZFv2/';
font_size = 19;
legend_font_size = 19;
no_size = 25;
line_width = 2;
%------------

build_grid;
rX = rhos0 * r;
nr_start = find(rX >= r_start, 1, 'first');
if nr_start < 2
	error(['r_start must be greater than ', num2str(dr*rhos0), ...
		' for the current grid']);
end
nr_end = find(rX <= r_end, 1, 'last');
nrdiag = nr_end - nr_start + 1;
tX = t0 * dt * nt_per_diagnose * (start_diag: end_diag);

%%
RS_pl_tavg = zeros(nrdiag, 1);
RS_fl_tavg = zeros(nrdiag, 1);
vr2_pl_tavg = zeros(nrdiag, 1);
vr2_fl_tavg = zeros(nrdiag, 1);
den_tavg = zeros(nrdiag, 1);
n_diag = end_diag - start_diag + 1;
vx_pl = zeros(nx+2, nx+2, nz+2);
vy_pl = zeros(nx+2, nx+2, nz+2);
vx_fl = zeros(nx+2, nx+2, nz+2);
vy_fl = zeros(nx+2, nx+2, nz+2);
disp(['t = ', num2str(tX(1)), ' ~ ', num2str(tX(end)), ' s'])

for idiag = start_diag: end_diag
	disp(['getting Renolds stress, step ', num2str(idiag), ' of ', ...
		num2str(end_diag)]);
	data_file = fullfile('data', ['dat', sprintf('%4.4d', idiag)]);
	load(data_file, 'Te', 'phi', 'den');
	denp = field2pol(den);
	den_tavg = den_tavg + squeeze(zonal_average(denp(nr_start: nr_end, :, :)));

	vx_pl(2:end-1, 2:end-1, 2:end-1) = 1/(2*dx) * ...
		(phi(2:end-1, 1:end-2, 2:end-1) - phi(2:end-1, 3:end, 2:end-1));
	vy_pl(2:end-1, 2:end-1, 2:end-1) = 1/(2*dx) * ...
		(phi(3:end, 2:end-1, 2:end-1) - phi(1:end-2, 2:end-1, 2:end-1));

	phi_fl = phi - Lambda * Te;
	vx_fl(2:end-1, 2:end-1, 2:end-1) = 1/(2*dx) * ...
		(phi_fl(2:end-1, 1:end-2, 2:end-1) - phi_fl(2:end-1, 3:end, 2:end-1));
	vy_fl(2:end-1, 2:end-1, 2:end-1) = 1/(2*dx) * ...
		(phi_fl(3:end, 2:end-1, 2:end-1) - phi_fl(1:end-2, 2:end-1, 2:end-1));
	[vr_pl, vtht_pl] = vec2pol(vx_pl, vy_pl);
	[vr_fl, vtht_fl] = vec2pol(vx_fl, vy_fl);
	
	RS_pl = zonal_average(vr_pl(nr_start: nr_end, :, :) .* vtht_pl(nr_start: nr_end, :, :));
	RS_fl = zonal_average(vr_fl(nr_start: nr_end, :, :) .* vtht_fl(nr_start: nr_end, :, :));
	RS_pl_tavg = RS_pl_tavg + squeeze(RS_pl);
	RS_fl_tavg = RS_fl_tavg + squeeze(RS_fl);

%	vr2_pl = zonal_average(vr_pl(nr_start: nr_end, :, :) .^ 2);
%	vr2_fl = zonal_average(vr_fl(nr_start: nr_end, :, :) .^ 2);
%	vr2_pl_tavg = vr2_pl_tavg + squeeze(vr2_pl);
%	vr2_fl_tavg = vr2_fl_tavg + squeeze(vr2_fl);
end	
den_tavg = den_tavg / n_diag;
RS_pl_tavg = RS_pl_tavg / n_diag;
RS_fl_tavg = RS_fl_tavg / n_diag;
%vr2_pl_tavg = vr2_pl_tavg / n_diag;
%vr2_fl_tavg = vr2_fl_tavg / n_diag;

%%
close;
fig = figure('name', 'Renolds stress');  
r_axis = rX(nr_start: nr_end);
xlim = [x_tick(1), x_tick(end)];


subplot(1, 2, 1);
plot(r_axis, denref * den_tavg, 'k-', 'LineWidth', line_width);
grid on;
xlabel('$$r$$/cm', 'interpreter', 'latex');
ylabel('$$\left<n\right>\ \rm\left(cm^{-3}\right)$$', ...
	'interpreter', 'latex');
set(gca, 'XLim', xlim);
set(gca, 'XTick', x_tick);
set(gca, 'FontSize', font_size);
ylim = get(gca, 'YLim');
txt = text(xlim(1), ylim(end), '(a)', 'VerticalAlignment', 'top');
set(txt, 'FontSize', no_size);
set(gca, 'YLim', ylim);

subplot(1,2,2);
plot(r_axis, 100*cs0^2*RS_pl_tavg, 'b-', 'LineWidth', line_width);  hold on;
plot(r_axis, cs0^2*RS_fl_tavg, 'r-', 'LineWidth', line_width);
grid on;
xlabel('$$r$$/cm', 'interpreter', 'latex');
ylabel(['$$\left<v_{E,r}v_{E,\theta}\right>\ ', ...
    '\left(\rm cm^{2}/s^{2}\right)$$'], ...
	'interpreter', 'latex');
set(gca, 'XLim', xlim);
set(gca, 'XTick', x_tick);
set(gca, 'FontSize', font_size);
lgd = legend('use $$10*\phi$$', 'use $$\phi_{f}$$');
legend('boxoff');
set(lgd, 'interpreter', 'latex');
set(lgd, 'FontSize', legend_font_size);
set(lgd, 'Location', 'NorthEast');
ylim = get(gca, 'YLim');
%txt = text(xlim(1), ylim(end), '(a)', 'VerticalAlignment', 'top');
txt = text(xlim(1), 2.5e9, '(b)', 'VerticalAlignment', 'top');
set(txt , 'FontSize', no_size);

%subplot(1,2,2);
%plot(r_axis, 25*cs0^2*vr2_pl_tavg, 'b-', 'LineWidth', line_width);  hold on;
%plot(r_axis, cs0^2*vr2_fl_tavg, 'r-', 'LineWidth', line_width);
%grid on;
%xlabel('$$r$$/cm', 'interpreter', 'latex');
%ylabel('$$\left<v_{E,r}^{2}\right>\ \rm \left(cm^{2}/s^{2}\right)$$', ...
%	'interpreter', 'latex');
%set(gca, 'XLim', xlim);
%set(gca, 'XTick', x_tick);
%set(gca, 'FontSize', font_size);
%lgd = legend('use $$5*\phi$$', 'use $$\phi_{fl}$$');
%legend('boxoff');
%set(lgd, 'interpreter', 'latex');
%set(lgd, 'FontSize', legend_font_size);
%set(lgd, 'Location', 'NorthEast');
%ylim = get(gca, 'YLim');
%set(gca, 'YLim', ylim);
%txt = text(xlim(1), ylim(end), '(b)', 'VerticalAlignment', 'top');
%set(txt , 'FontSize', no_size);

set(fig, 'Position', fig_position);
disp(['t = ', num2str(tX(1)), ' ~ ', num2str(tX(end)), ' s'])
pic_name = fullfile('figs', ...
	['t=', num2str(start_diag), '-', num2str(end_diag), 'tavg_den_RS_f']);
print(gcf, '-dpng', [pic_name, '.png']);
print(gcf, '-depsc', [pic_name, '.eps']);
if exist(dst_path, 'dir')
	copyfile([pic_name, '.eps'], dst_path);
else
	disp('Warning: Destination path do not exist');
end
