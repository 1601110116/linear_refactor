% This script draws the azimuthal Renolds' stress using plasma potential and floating potential

clear;  close all;
global Te phi

load ('parameters.mat');
addpath(code_path);
last_file = get_last_file('./');
last_diag = str2num(last_file(end-7: end-4));
%%
%---input---
start_diag = 303;
end_diag = last_diag;
r_start = 0.4;
r_end = 9;
font_size = 20;
Lambda = 3;  % phi_floating = phi - Lambda*Te

fig_position = [50, 50, 1000, 800];  % [left, bottom, width, height]
x_tick = 1:2:9;
y_tick = -1e9:0.5e9:2e9;
dst_path = 'Documentation/KH/';
font_size = 24;
legend_font_size = 24;
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
Renolds_stress_plasma_tavg = zeros(nrdiag, 1);
Renolds_stress_floating_tavg = zeros(nrdiag, 1);
zf_plasma_tavg = zeros(nrdiag, 1);
zf_floating_tavg = zeros(nrdiag, 1);
n_diag = end_diag - start_diag + 1;
vx_plasma = zeros(nx+2, nx+2, nz+2);
vy_plasma = zeros(nx+2, nx+2, nz+2);
vx_floating = zeros(nx+2, nx+2, nz+2);
vy_floating = zeros(nx+2, nx+2, nz+2);

for idiag = start_diag: end_diag
	disp(['getting Renolds stress, step ', num2str(idiag), ' of ', ...
		num2str(end_diag)]);
	data_file = ['dat', sprintf('%4.4d', idiag)];
	load(data_file);

	vx_plasma(2:end-1, 2:end-1, 2:end-1) = 1/(2*dx) * ...
		(phi(2:end-1, 1:end-2, 2:end-1) - phi(2:end-1, 3:end, 2:end-1));
	vy_plasma(2:end-1, 2:end-1, 2:end-1) = 1/(2*dx) * ...
		(phi(3:end, 2:end-1, 2:end-1) - phi(1:end-2, 2:end-1, 2:end-1));

	phi_floating = phi - Lambda * Te;
	vx_floating(2:end-1, 2:end-1, 2:end-1) = 1/(2*dx) * ...
		(phi_floating(2:end-1, 1:end-2, 2:end-1) - phi_floating(2:end-1, 3:end, 2:end-1));
	vy_floating(2:end-1, 2:end-1, 2:end-1) = 1/(2*dx) * ...
		(phi_floating(3:end, 2:end-1, 2:end-1) - phi_floating(1:end-2, 2:end-1, 2:end-1));
	[vr_plasma, vtht_plasma] = vec2pol(vx_plasma, vy_plasma);
	[vr_floating, vtht_floating] = vec2pol(vx_floating, vy_floating);
	
	Renolds_stress_plasma = zonal_average(vr_plasma(nr_start: nr_end, :, :) .* vtht_plasma(nr_start: nr_end, :, :));
	Renolds_stress_floating = zonal_average(vr_floating(nr_start: nr_end, :, :) .* vtht_floating(nr_start: nr_end, :, :));
	Renolds_stress_plasma_tavg = Renolds_stress_plasma_tavg + squeeze(Renolds_stress_plasma);
	Renolds_stress_floating_tavg = Renolds_stress_floating_tavg + squeeze(Renolds_stress_floating);
	zf_plasma_tavg = zf_plasma_tavg + zonal_average(vtht_plasma(nr_start: nr_end, : ,:));
	zf_floating_tavg = zf_floating_tavg + zonal_average(vtht_floating(nr_start: nr_end, :, :));
end	
Renolds_stress_plasma_tavg = Renolds_stress_plasma_tavg / n_diag;
Renolds_stress_floating_tavg = Renolds_stress_floating_tavg / n_diag;
zf_plasma_tavg = zf_plasma_tavg / n_diag;
zf_floating_tavg = zf_floating_tavg / n_diag;

%%
close;
fig = figure('name', 'Renolds stress');  
r_axis = rX(nr_start: nr_end);
subplot(1, 2, 1);
plot(r_axis, 100*cs0^2*Renolds_stress_plasma_tavg, 'b-');  hold on;
plot(r_axis, cs0^2*Renolds_stress_floating_tavg, 'r-');
set(gca, 'XTick', x_tick);  set(gca, 'YTick', y_tick);
set(gca, 'FontSize', font_size);
lgd = legend('$$\left<v_{r,pl}v_{\theta,pl}\right>$$', '$$\left<v_{r,fl}v_{\theta,fl}\right>$$');
set(lgd, 'interpreter', 'latex');
set(lgd, 'FontSize', legend_font_size);

subplot(1,2,2);
plot(r_axis, cs0*zf_plasma_tavg, 'b-');  hold on;
plot(r_axis, cs0*zf_floating_tavg, 'r-'); 
set(gca, 'XTick', x_tick);
set(gca, 'FontSize', font_size);
lgd = legend('$$\left<v_{\theta,pl}\right>$$', '$$\left<v_{\theta,fl}\right>$$');
set(lgd, 'interpreter', 'latex');
set(lgd, 'FontSize', legend_font_size);

set(fig, 'Position', fig_position);
pic_name = ['t=', num2str(start_diag), '-', num2str(end_diag), '_Renolds_stress'];
print(gcf, '-dpng', [pic_name, '.png']);
print(gcf, '-depsc', [pic_name, '.eps']);
if exist(dst_path, 'dir')
	copyfile([pic_name, '.eps'], dst_path);
else
	disp('Warning: Destination path do not exist');
end
