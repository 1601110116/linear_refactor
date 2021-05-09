% This script draws time-averaged profiles of zonal flow drives 
%  versus time1
% Copy this script to the folder of the data and then run it

clear;  close all;
global den Te pe vi jz ve phi vEx vEy dt inv_nustar

%%
%---input---
r_start = 0.1;  % cm
r_end = 10;  % cm
interval_nr = 3;  % for every interval_nr radial points, only one point of data is used
init_diag = 814;
font_size = 30;
fig_size = [50, 50, 1000, 1400]; %[left bottom width height]
x_lim = [25.5, 25.9];
x_tick = 25.5:0.1:25.9;
y_tick = 0:3:9;
vE_z_tick = 0:3:9;
vi_z_tick = -50:25:50;
view_angle = [-45, 45];  % [horizontal rotation, vertical elevation]
dst_path = '../../Documentation/KH';  % specify where the .eps output should be copied to
%-----------
%%

load('../../parameters.mat');
init_time = t0 * dt * nt_per_diagnose * init_diag;
load('parameters.mat');
%code_path='../../code';
addpath(code_path);
last_file = get_last_file('./');
last_diag = str2num(last_file(end-7: end-4));

build_grid;
rX = rhos0 * r;
nr_start = find(rX >= r_start, 1, 'first');
nr_end = find(rX <= r_end, 1, 'last');
nrdiag = nr_end - nr_start + 1;

%%
jz = zeros (nx+2, nx+2, nz+2);
ve = jz;  vEx = jz;  vEy = jz;
generate_constants(height, radius, dx, dz, nx, nz, dt, ... 
    rdif, dif_perp_in, dif_perp_out, dif_z_in, dif_z_out, ... 
    rconduct, conduct_z_in, conduct_z_out, viscosity, ... 
    den_damp, momentum_damp, Te_tanhsrc_max, Te_tanhsrc_radius, ... 
    Te_tanhsrc_incline, Te_gausssrc_magnitude, Te_gausssrc_sigma, ... 
    den_tanhsrc_max, den_tanhsrc_radius, den_tanhsrc_incline, ... 
    den_gausssrc_magnitude, den_gausssrc_sigma);

vE_profile_t = zeros(nrdiag, last_diag+1);
vi_profile_t = vE_profile_t;

for idiag = 0: last_diag
	disp(['getting zf profile: step ', num2str(idiag), ' of ', ...
		num2str(last_diag)]);
	data_file = ['dat', sprintf('%4.4d', idiag)];
	load(data_file);
	sdata();

	[vErp, vEthtp] = vec2pol(vEx, vEy);
	vip = field2pol(vi);

	vE_profile = squeeze(zonal_average(vEthtp(nr_start: nr_end, :, :)));
	vE_profile_t(:, idiag+1) = vE_profile;
	vi_profile = squeeze(zonal_average(vip(nr_start: nr_end, :, :)));
	vi_profile_t(:, idiag+1) = vi_profile;
end
%%
close;
fig = figure('name', 'ZF drives');
set(fig, 'Position', fig_size);
r_axis = rX(nr_start: nr_end);
normalize_factor = cs0/(t0*rhos0);
%tX = t0 * dt * nt_per_diagnose * (start_diag: end_diag);
tX = init_time + t0 * dt * nt_per_diagnose * (0:last_diag);

subplot(2,1,1);
s1 = surf(1e3*tX, r_axis(1:interval_nr:end), 1e-4*cs0*vE_profile_t(1:interval_nr:end,:));
set(s1, 'FaceAlpha', 0.2, 'EdgeColor', 'b', 'FaceColor', 'b');
title('zonal flow profile', 'FontWeight', 'normal');
xlabel('$$t$$/ms', 'interpreter', 'latex');
ylabel('$$r$$/cm', 'interpreter', 'latex');
zlabel('$$\left<v_{E,\theta}\right>/\left(10^{4}\rm{cm/s}\right)$$', 'interpreter', 'latex');
set(gca, 'XLim', x_lim);
set(gca, 'XTick', x_tick);
set(gca, 'YTick', y_tick);
set(gca, 'ZTick', vE_z_tick);
view(view_angle);
set(gca, 'fontSize', font_size);

subplot(2,1,2);
s2 = surf(1e3*tX, r_axis(1:interval_nr:end), 1e-2*cs0*vi_profile_t(1:interval_nr:end,:));
set(s2, 'FaceAlpha', 0.2, 'EdgeColor', 'g', 'FaceColor', 'g');
title('parallel flow profile', 'FontWeight', 'normal');
xlabel('$$t$$/ms', 'interpreter', 'latex');
ylabel('$$r$$/cm', 'interpreter', 'latex');
zlabel('$$\left<v_{\parallel i}\right>/\left(\rm{m/s}\right)$$', 'interpreter', 'latex');
set(gca, 'XLim', x_lim);
set(gca, 'XTick', x_tick);
set(gca, 'YTick', y_tick);
set(gca, 'ZTick', vi_z_tick);
view(view_angle);
set(gca, 'fontSize', font_size);

print(gcf, '-dpng', 'flow_growth.png');
print(gcf, '-depsc', 'flow_growth.eps');

if exist(dst_path)
	copyfile('flow_growth.eps', dst_path);
else
	disp('Destination path do not exist');
end
