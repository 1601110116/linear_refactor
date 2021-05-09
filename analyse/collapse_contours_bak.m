%  This script shows the den and phi contours in the collapse phase

clear;  close all;
global den phi

%%
load('parameters.mat');
addpath(code_path);


%---input---
start_diag = 499;  %499 818
% pictures of one certain field is listed num_rows*num_columns
num_columns = 6;
num_rows = 4;
zdiag = ceil(length(zX)/2);  % middle x-y plane
font_size = 15;
no_size = 12;
cbar_label_size = 18;
fig_width = 800;
margin = 0.0025;
x_tick = -10:5:10;
cbar_width = 0.02;
dst_path = 'Documentation/KH/';
%-----------
%%

total_num_rows = num_rows * 3;  % three fields are displayed, i.e. den, phi and w
margins = [margin, total_num_rows/num_columns * margin];  % [vertical, horizontal]
% fig_position = [left, bottom, width, height];
fig_position = [50, 50, fig_width, total_num_rows/num_columns * fig_width];
% the total number of diagnoses is the column number times the row number
n_diag = num_columns * num_rows - 1;  % the last gca is used to draw colobar
end_diag = start_diag + n_diag - 1;
[xX2d, yX2d] = ndgrid(xX(2:end-1), xX(2:end-1));  % excludes boundaries
rX2d = sqrt(xX2d.^2 + yX2d.^2);
out2d = rX2d > (radius*rhos0);  % get the indices of the outside grid points
show_ylabel = 1:num_columns:n_diag;  % index of gca that should show yticklabel and ylabel
show_xlabel = num_columns*(num_rows-1)+1:n_diag;
init_time = 1e3 * t0 * dt * nt_per_diagnose * start_diag;
final_time = 1e3 * t0 * dt * nt_per_diagnose * end_diag;
time_interval = 1e3 * t0 * dt * nt_per_diagnose;
disp(['time of the first picture: ', num2str(init_time), ' ms']);
disp(['time of the last picture: ', num2str(final_time), ' ms']);
disp(['time interval between adjacent pictures: ', num2str(time_interval), ' ms']);

% find the CLim of den, phi and w 
load(['dat', sprintf('%4.4d', start_diag)], 'den', 'phi', 'w')
tmp = floor(nx/2);
den_clim = [0, denref*den(tmp, tmp, zdiag)];
phi_clim = Tref * phi(tmp, tmp, zdiag) .* ones(2, 1);
w_clim = cs0/rhos0 * w(tmp, tmp, zdiag) .* ones(2, 1);
for idiag = start_diag: end_diag
	data_file = ['dat', sprintf('%4.4d', idiag)];
	load(data_file, 'den', 'phi', 'w');
	den_max = denref * max(max(den(:, :, zdiag), [], 1), [], 2);
	phi_min = Tref * min(min(phi(:, :, zdiag), [], 1), [], 2);
	phi_max = Tref * max(max(phi(:, :, zdiag), [], 1), [], 2);
	w_min = cs0/rhos0 * min(min(w(:, :, zdiag), [], 1), [], 2);
	w_max = cs0/rhos0 * max(max(w(:, :, zdiag), [], 1), [], 2);
	if den_max > den_clim(end)
		den_clim(end) = den_max;
	end
	if phi_min < phi_clim(1)
		phi_clim(1) = phi_min;
	end
	if phi_max > phi_clim(2)
		phi_clim(2) = phi_max;
	end
	if w_min < w_clim(1)
		w_clim(1) = w_min;
	end
	if w_max > w_clim(end)
		w_clim(end) = w_max;
	end
end

fig = figure;

for idiag = start_diag: end_diag
   data_file = ['dat', sprintf('%4.4d', idiag)];
   load(data_file, 'den', 'phi', 'w');
   
   igca = idiag-start_diag+1;  % indiex of the current gca
   subplot_tight(3*num_rows, num_columns, igca, margins);
   den_disp = den(2:end-1, 2:end-1, zdiag);
   den_disp(out2d) = nan;
   pcolor(yX2d, xX2d, denref*den_disp);
   colormap jet;  shading interp;
   set(gca, 'CLim', den_clim);
   set(gca, 'XTick', x_tick);
   set(gca, 'YTick', x_tick);
   set(gca, 'xticklabel', []);
   set(gca, 'yticklabel', []);
   text(x_tick(1), x_tick(end), ['\fontsize{', num2str(no_size), '}', num2str(igca)], ...
	   'VerticalAlignment', 'top', 'color', 'black');
   
   subplot_tight(3*num_rows, num_columns, igca+n_diag+1, margins);
   phi_disp = phi(2:end-1, 2:end-1, zdiag);
   phi_disp(out2d) = nan;
   pcolor(yX2d, xX2d, Tref*phi_disp);
   colormap jet;  shading interp;
   set(gca, 'CLim', phi_clim);
   set(gca, 'XTick', x_tick);
   set(gca, 'YTick', x_tick);
   set(gca, 'xticklabel', []);
   set(gca, 'yticklabel', []);
   text(x_tick(1), x_tick(end), ['\fontsize{', num2str(no_size), '}', num2str(igca)], ...
	   'VerticalAlignment', 'top', 'color', 'black');

   subplot_tight(3*num_rows, num_columns, igca+2*(n_diag+1), margins);
   w_disp = w(2:end-1, 2:end-1, zdiag);
   w_disp(out2d) = nan;
   pcolor(yX2d, xX2d, cs0/rhos0*w_disp);
   colormap jet;  shading interp;
   set(gca, 'CLim', w_clim);
   set(gca, 'XTick', x_tick);
   set(gca, 'YTick', x_tick);
   set(gca, 'xticklabel', []);
   set(gca, 'yticklabel', []);
   text(x_tick(1), x_tick(end), ['\fontsize{', num2str(no_size), '}', num2str(igca)], ...
	   'VerticalAlignment', 'top', 'color', 'black');

end
subplot_tight(3*num_rows, num_columns, n_diag+1, margins);
set(gca, 'Visible', 'off');
colormap jet;  cbar = colorbar;
set(cbar, 'AxisLocation', 'in');
set(gca, 'CLim', 1e-12*den_clim);
set(cbar.Label, 'String', '$$n/10^{12}\rm{cm^{-3}}$$');
set(cbar.Label, 'interpreter', 'latex');
cbar_position = get(cbar, 'Position');
cbar_position(3) = cbar_width;
cbar_position(1) = cbar_position(1) + cbar_width;
set(cbar, 'Position', cbar_position);
set(gca, 'fontSize', font_size);
set(cbar.Label, 'fontSize', cbar_label_size);

subplot_tight(3*num_rows, num_columns, 2*(n_diag+1), margins);
set(gca, 'Visible', 'off');
colormap jet;  cbar = colorbar;
set(cbar, 'AxisLocation', 'in');
set(gca, 'CLim', phi_clim);
set(cbar.Label, 'String', '$$\phi/\rm{V}$$');
set(cbar.Label, 'interpreter', 'latex');
cbar_position = get(cbar, 'Position');
cbar_position(3) = cbar_width;
cbar_position(1) = cbar_position(1) + cbar_width;
set(cbar, 'Position', cbar_position);
set(gca, 'fontSize', font_size);
set(cbar.Label, 'FontSize', cbar_label_size);

subplot_tight(3*num_rows, num_columns, 3*(n_diag+1), margins);
set(gca, 'Visible', 'off');
colormap jet;  cbar = colorbar;
set(cbar, 'AxisLocation', 'in');
set(gca, 'CLim', 1e-4*w_clim);
set(cbar.Label, 'String', '$$w/10^{4}\rm{s^{-1}}$$');
set(cbar.Label, 'interpreter', 'latex');
cbar_position = get(cbar, 'Position');
cbar_position(3) = cbar_width;
cbar_position(1) = cbar_position(1) + cbar_width;
set(cbar, 'Position', cbar_position);
set(gca, 'fontSize', font_size);
set(cbar.Label, 'FontSize', cbar_label_size);

set(fig, 'Position', fig_position);
print(gcf, '-dpng', ['collapse_from_t=', sprintf('%4.4d', start_diag), '.png']);
print(gcf, '-depsc', ['collapse_from_t=', sprintf('%4.4d', start_diag), '.eps']);
if exist(dst_path, 'dir')
	copyfile(['collapse_from_t=', sprintf('%4.4d', start_diag), '.eps'], dst_path);
else
	disp('Destination path do not exist');
end
