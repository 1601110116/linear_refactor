% This script draws time-averaged profiles of zonal flow drives 
%  versus time1
% Copy this script to the folder of the data and then run it

clear;  close all;
global den Te pe vi jz ve phi vEx vEy dt inv_nustar

load('parameters.mat');
addpath(code_path);
last_file = get_last_file('data');
last_diag = str2num(last_file(end-7: end-4));

%---input---
start_diag = 1035;
end_diag = 1045;
r_start = 0.4;  % cm
r_end = 6;  % cm
font_size = 15;
legend_font_size = 13;
fig_size = [50, 50, 1500, 1500]; %[left bottom width height]
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

zf_shear_t = zeros(nrdiag, n_diag);
zf_shear_polar_t = zeros(nrdiag, n_diag);
perp_drive_ex_t = zf_shear_t;  para_drive_ex_t = zf_shear_t;
sum_t = zeros(nrdiag-4, n_diag);

for idiag = start_diag: end_diag
	disp(['getting zf drives: step ', num2str(idiag), ' of ', ...
		num2str(end_diag)]);
	data_file = fullfile('data', ['dat', sprintf('%4.4d', idiag)]);
	load(data_file);
	sdata();

	denp = field2pol(den);
	jzp = field2pol(jz);
	[vErp, vEthtp] = vec2pol(vEx, vEy);

	% radial derivative of zonal flow
	zf_shear = 1/(2*dr) .* squeeze( ... 
		zonal_average(vEthtp(nr_start+1: nr_end+1, :, :)) ... 
		- zonal_average(vEthtp(nr_start-1: nr_end-1, :, :)));
	zf_shear_t(:, idiag-start_diag+1) = zf_shear;

	% the additional term in the laplacian operator is added to zf shear
	r_tmp = reshape(r(nr_start-1: nr_end+1), nrdiag+2, 1);
	zf_shear_polar = 1/(2*dr) ./r_tmp(2:end-1) .* squeeze( ...
		r_tmp(3:end) .* zonal_average(vEthtp(nr_start+1: nr_end+1, :, :)) ...
		- r_tmp(1:end-2) .* zonal_average(vEthtp(nr_start-1: nr_end-1, :, :)));
	zf_shear_polar_t(:, idiag-start_diag+1) = zf_shear_polar;
	
	% exact perpendicular drive
	r_tmp = repmat(reshape(r(nr_start-2: nr_end+2), nrdiag+4, 1), 1, ntht, nz+2);
	tmp1 = r_tmp .* zonal_pert(vEthtp(nr_start-2: nr_end+2, :, :));
	tmp1 = 1/(2*dr) * (tmp1(3:end, :, :) - tmp1(1:end-2, :, :));
	tmp1 = 1./r_tmp(2:end-1, :, :) .* tmp1;
	tmp1 = squeeze(zonal_average( ...
		zonal_pert(vErp(nr_start-1: nr_end+1, :, :)) .* tmp1));
	perp_drive_ex = -radius/(2*dr) * (tmp1(3:end) - tmp1(1:end-2));
	perp_drive_ex_t(:, idiag-start_diag+1) = perp_drive_ex;

	% exact parallel drive
	tmp1 = zeros(nrdiag, ntht, nz+2);
	tmp1(:, :, 2:end-1) = (radius/height) / (2*dz) ...
		* (zonal_pert(jzp(nr_start: nr_end, :, 3:end)) ...
	   	- zonal_pert(jzp(nr_start: nr_end, :, 1:end-2)));
	para_drive_ex = squeeze(zonal_average(tmp1./denp(nr_start: nr_end, :, :)));
	para_drive_ex_t(:, idiag-start_diag+1) = para_drive_ex;

	% viscosity contribution
	r_tmp = reshape(r(nr_start+1: nr_end-1), [], 1);
	tmp1 = 1/(2*dr)*r_tmp .* (zf_shear_polar(3:end) - zf_shear_polar(1:end-2));
	shear_diff = viscosity./(2*dr)./r_tmp(2:end-1).*(tmp1(3:end) - tmp1(1:end-2));
	shear_damp = -momentum_damp * zf_shear_polar;
	sum_t(:, idiag-start_diag+1) = perp_drive_ex(3:end-2) + para_drive_ex(3:end-2) ...
		+ shear_diff + shear_damp(3:end-2);
end
%%
close;
dtS_minus_sum = 1/(2*dt*nt_per_diagnose) ...
	* (zf_shear_polar_t(3:end-2, 3:end) - zf_shear_polar_t(3:end-2, 1:end-2)) - sum_t(:, 2:end-1);
figure('name', 'ZF drives');
set(gcf, 'Position', fig_size);
r_axis = rX(nr_start: nr_end);
normalize_factor = cs0/(t0*rhos0);
tX = t0 * dt * nt_per_diagnose * (start_diag: end_diag);

subplot(2,1,1);
s1 = surf(1e3*tX, r_axis, cs0/rhos0*zf_shear_t);  hold on;
s2 = surf(1e3*tX, r_axis, cs0/rhos0*zf_shear_polar_t);
set(s1, 'FaceAlpha', 0.2, 'EdgeColor', 'c', 'FaceColor', 'c');
set(s2, 'FaceAlpha', 0.2, 'EdgeColor', 'g', 'FaceColor', 'g');
title('zonal flow shear');
xlabel('t/ms');
ylabel('r/cm');
zlabel('shear$$/\rm{s^{-1}}$$', 'interpreter', 'latex');
set(gca, 'fontSize', font_size);
lgd = legend('$$\frac{\partial}{\partial r}\left<v_{E,\theta}\right>$$', ...
	['$$\frac{1}{r}\frac{\partial}{\partial r}\left(r\left<v_{E,\theta}', ...
	'\right>\right)$$'], 'Location', 'NorthEast');
set(lgd, 'interpreter', 'latex');
set(lgd, 'FontSize', legend_font_size);

subplot(2,1,2);
s3 = surf(1e3*tX, r_axis, normalize_factor*perp_drive_ex_t);  hold on;
s4 = surf(1e3*tX, r_axis, normalize_factor*para_drive_ex_t);
s5 = surf(1e3*tX(2:end-1), r_axis(3:end-2), normalize_factor*dtS_minus_sum);
set(s3, 'FaceAlpha', 0.2, 'EdgeColor', 'r', 'FaceColor', 'r');
set(s4, 'FaceAlpha', 0.2, 'EdgeColor', 'b', 'FaceColor', 'b');
set(s5, 'FaceAlpha', 0.3, 'EdgeColor', 'k', 'FaceColor', 'k');
xlabel('t/ms');
ylabel('r/cm');
zlabel('shear\ drive$$/\rm{s^{-2}}$$', 'interpreter', 'latex');
title('drives of zonal flow shear');
set(gca, 'fontSize', font_size);

lgd = legend('$$F_{\perp}^{\dagger}$$', '$$F_{\parallel}^{\dagger}$$', ...
	'$$\rm{d}_{t}S^{\dagger}-\rm{sum}$$');
%lgd = legend('$$F_{\perp}^{\dagger}$$', '$$F_{\parallel}^{\dagger}$$');

set(lgd, 'interpreter', 'latex');
set(lgd, 'FontSize', legend_font_size);

file_name = fullfile('figs', ['t=', num2str(start_diag), '-',num2str(end_diag), ...
	'flow_shear_drives_t']);
print(gcf, '-dpng', [file_name, '.png'])
print(gcf, '-depsc', [file_name, '.eps'])
