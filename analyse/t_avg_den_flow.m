% This script draws the time-averaged azimuthal EXB flow

clear;  close all;
global Te phi

load ('parameters.mat');
addpath(code_path);
last_file = get_last_file('./');
last_diag = str2num(last_file(end-7: end-4));
%%
%---input---
start_diag = 600;
end_diag = last_diag;
r_start = 0.4;
r_end = 9;

fig_position = [50, 50, 1300, 500];  % [left, bottom, width, height]
x_tick = 1:1:7.5;
x_lim = [0, 7.5];
line_width = 2;
dst_path = 'Documentation/ZF/';
font_size = 24;
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

%%
zf_tavg = zeros(nrdiag, 1);
den_tavg = zeros(nrdiag, 1);
n_diag = end_diag - start_diag + 1;
vx = zeros(nx+2, nx+2, nz+2);
vy = zeros(nx+2, nx+2, nz+2);

for idiag = start_diag: end_diag
	disp(['getting Renolds stress, step ', num2str(idiag), ' of ', ...
		num2str(end_diag)]);
	data_file = ['dat', sprintf('%4.4d', idiag)];
	load(data_file, 'phi', 'den');

	vx(2:end-1, 2:end-1, 2:end-1) = 1/(2*dx) * ...
		(phi(2:end-1, 1:end-2, 2:end-1) - phi(2:end-1, 3:end, 2:end-1));
	vy(2:end-1, 2:end-1, 2:end-1) = 1/(2*dx) * ...
		(phi(3:end, 2:end-1, 2:end-1) - phi(1:end-2, 2:end-1, 2:end-1));

	[vr, vtht] = vec2pol(vx, vy);
	
	zf_tavg = zf_tavg + zonal_average(vtht(nr_start: nr_end, : ,:));
	denp = field2pol(den);
	den_tavg =  den_tavg + zonal_average(denp(nr_start: nr_end, :, :));
end	
zf_tavg = zf_tavg / n_diag;
den_tavg = den_tavg / n_diag;
%%
close;
fig = figure('name', 'time averaged ZF');  
r_axis = rX(nr_start: nr_end);

subplot(1,2,1)
plot(r_axis, denref * 1e-12 * den_tavg, 'b-', 'lineWidth', line_width)
xlabel('$$r$$/cm', 'interpreter', 'latex');
ylabel('$$\left<n\right>\left(10^{12}\rm{cm^-3}\right)$$', 'interpreter', 'latex');
set(gca, 'XTick', x_tick);
set(gca, 'XLim', x_lim);
set(gca, 'FontSize', font_size);

subplot(1,2,2)
plot(r_axis, 1e-2*cs0*zf_tavg, 'b-', 'lineWidth', line_width);
xlabel('$$r$$/cm', 'interpreter', 'latex');
ylabel('$$\left<v_{E,\theta}\right>\left(\rm{m/s}\right)$$', 'interpreter', 'latex');
set(gca, 'XTick', x_tick);
set(gca, 'XLim', x_lim);
set(gca, 'FontSize', font_size);

set(fig, 'Position', fig_position);
pic_name = ['t=', num2str(start_diag), '-', num2str(end_diag), '_ZF'];
print(gcf, '-dpng', [pic_name, '.png']);
print(gcf, '-depsc', [pic_name, '.eps']);
if exist(dst_path, 'dir')
	copyfile([pic_name, '.eps'], dst_path);
else
	disp('Warning: Destination path do not exist');
end
