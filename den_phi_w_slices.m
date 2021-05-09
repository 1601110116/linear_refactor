%  This script draws 3D structures of density, potential and vorticity by slices

clear;  close all;

%%
load('parameters.mat');
addpath(code_path);


%---input---
%  save_eps = 1: figures are saved as .eps files and copied to dst_path
%  save_eps = 0; figures are saved as .png files and copied to slices/
save_eps = 0;
dst_path = 'Documentation/KH/';

start_diag = 1541; %1549; %1138; %1146;  % 514; %506
end_diag = 1548; %1551; %1140; %1148;
nz_per_slice = 3;  % draw a slice every nz_per_diag x-y planes;

fig_position = [50, 50, 1200, 1500];  % [left, bottom, width, height];
view_angle = [13, 30];  % [horizontal rotation, vertical elevation]
font_size = 20;
no_size = 25;
title_size = 27;
line_width = 3;
cbar_right_shift = 0.02;
x_tick = zX(2:nz_per_slice:end-1);
y_tick = -10: 5: 10;
z_tick = -10: 5: 10;
%--------------

x_tick = roundn(x_tick, -1);
if ~exist('3D_structures', 'dir')
	mkdir('3D_structures');
end
[zX3d, yX3d, xX3d] = meshgrid(zX(2:end-1), xX(2:end-1), xX(2:end-1));
rX3d = sqrt(xX3d.^2 + yX3d.^2);
out3d = rX3d > (radius*rhos0);
xslice = zX(2: nz_per_slice: end-1);  % zX of the x-y planes to diagnose
x_lim = [x_tick(1), x_tick(end)];
y_lim = [y_tick(1), y_tick(end)];
z_lim = [z_tick(1), z_tick(end)];


for idiag = start_diag: end_diag
	close;
	data_file = fullfile('data', ['dat', sprintf('%4.4d', idiag)]);
	load(data_file, 'den', 'phi', 'w');
	tX_ms = 1e3 * t0 * dt * nt_per_diagnose * idiag;
	disp(['painting slices, step', num2str(idiag), ' of ', num2str(idiag), ', t = ', num2str(tX_ms, 6), ' ms']);
	fig = figure;
	set(fig, 'Position', fig_position);
	subplot(3,1,1);
	den_disp = permute(den(2:end-1, 2:end-1, 2:end-1), [2, 3, 1]);
	den_disp(out3d) = nan;
	plot3([x_lim(1), x_lim(end)], [0, 0], [0, 0], 'w-', 'LineWidth', line_width);  
	hold on;
	slice(zX3d, yX3d, xX3d, denref*den_disp, xslice, [], [], 'nearest');
	colormap jet;  cbar=colorbar;  shading interp;
	hold off;
	view(view_angle);
	cbar_position = get(cbar, 'Position');  %[left, bottom, width, height]
	cbar_position(1) = cbar_position(1) + cbar_right_shift;
	gca_position = get(gca, 'Position');
	set(cbar, 'Position', cbar_position);
	set(gca, 'Position', gca_position);
	set(gca, 'XTick', x_tick);
	set(gca, 'YTick', y_tick);
	set(gca, 'ZTick', z_tick);
	set(gca, 'XLim', x_lim);
	set(gca, 'YLim', y_lim);
	set(gca, 'ZLim', z_lim);
	set(gca, 'XDir', 'reverse');
	xlabel('$$z$$/cm', 'interpreter', 'latex');
	ylabel('$$y$$/cm', 'interpreter', 'latex');
	zlabel('$$x$$/cm', 'interpreter', 'latex');
	set(gca, 'FontSize', font_size);
	title('$$n/\rm cm^{-3}$$', 'interpreter', 'latex', 'FontSize', title_size);
	text(x_lim(end), y_lim(end), z_lim(end), ['\fontsize{', num2str(no_size), '}(a)'], ...
		'VerticalAlignment', 'top', 'color', 'black');
    grid on;

	subplot(3,1,2);
	phi_disp = permute(phi(2:end-1, 2:end-1, 2:end-1), [2, 3, 1]);
	phi_disp(out3d) = nan;
	plot3([x_lim(1), x_lim(end)], [0, 0], [0, 0], 'w-', 'LineWidth', line_width);  
	hold on;
	slice(zX3d, yX3d, xX3d, Tref*phi_disp, xslice, [], [], 'nearest');
	colormap jet;  cbar=colorbar;  shading interp;
	hold off;
	view(view_angle);
	cbar_position = get(cbar, 'Position');  %[left, bottom, width, height]
	cbar_position(1) = cbar_position(1) + cbar_right_shift;
	gca_position = get(gca, 'Position');
	set(cbar, 'Position', cbar_position);
	set(gca, 'Position', gca_position);
	set(gca, 'XTick', x_tick);
	set(gca, 'YTick', y_tick);
	set(gca, 'ZTick', z_tick);
	set(gca, 'XLim', x_lim);
	set(gca, 'YLim', y_lim);
	set(gca, 'ZLim', z_lim);
	set(gca, 'XDir', 'reverse');
	xlabel('$$z$$/cm', 'interpreter', 'latex');
	ylabel('$$y$$/cm', 'interpreter', 'latex');
	zlabel('$$x$$/cm', 'interpreter', 'latex');
	set(gca, 'FontSize', font_size);
	title('$$\phi/\rm V$$', 'interpreter', 'latex', 'FontSize', title_size);
	text(x_lim(end), y_lim(end), z_lim(end), ['\fontsize{', num2str(no_size), '}(b)'], ...
		'VerticalAlignment', 'top', 'color', 'black');
	grid on;

	subplot(3,1,3);
	w_disp = permute(w(2:end-1, 2:end-1, 2:end-1), [2, 3, 1]);
	w_disp(out3d) = nan;
	plot3([x_lim(1), x_lim(end)], [0, 0], [0, 0], 'w-', 'LineWidth', line_width);  
	hold on;
	slice(zX3d, yX3d, xX3d, cs0/rhos0*w_disp, xslice, [], [], 'nearest');
	colormap jet;  cbar=colorbar;  shading interp;
	hold off;
	view(view_angle);
	cbar_position = get(cbar, 'Position');  %[left, bottom, width, height]
	cbar_position(1) = cbar_position(1) + cbar_right_shift;
	gca_position = get(gca, 'Position');
	set(cbar, 'Position', cbar_position);
	set(gca, 'Position', gca_position);
	set(gca, 'XTick', x_tick);
	set(gca, 'YTick', y_tick);
	set(gca, 'ZTick', z_tick);
	set(gca, 'XLim', x_lim);
	set(gca, 'YLim', y_lim);
	set(gca, 'ZLim', z_lim);
	set(gca, 'XDir', 'reverse');
	xlabel('$$z$$/cm', 'interpreter', 'latex');
	ylabel('$$y$$/cm', 'interpreter', 'latex');
	zlabel('$$x$$/cm', 'interpreter', 'latex');
	set(gca, 'FontSize', font_size);
	title('$$w/\rm s^{-1}$$', 'interpreter', 'latex', 'FontSize', title_size);
	text(x_lim(end), y_lim(end), z_lim(end), ['\fontsize{', num2str(no_size), '}(c)'], ...
		'VerticalAlignment', 'top', 'color', 'black');
	grid on;

	if save_eps == 1
		print(gcf, '-depsc', fullfile('figs', ['slices_t=', sprintf('%4.4d', idiag), '.eps']));
    elseif save_eps == 0
		print(gcf, '-dpng', fullfile('figs', ['3D_structures/3D_t=', sprintf('%4.4d', idiag), '.png']));
	else
		error('the input save_eps has to be 0 or 1');
	end
end


