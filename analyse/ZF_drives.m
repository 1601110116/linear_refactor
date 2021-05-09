% This script draws time-averaged profiles of zonal flow drives using data 
%  between start_diag and end_diag in its directory
% Copy this script to the folder of the data and then run it

clear;  close all;
global den Te pe vi jz ve phi vEx vEy dt inv_nustar

load('parameters.mat');
addpath(code_path);
last_file = get_last_file('./');
last_diag = str2num(last_file(end-7: end-4));

%---input---
start_diag = 814;
end_diag = 814;
r_start = 0.4;  % cm
r_end = 9.5;  % cm
font_size = 15;
legend_font_size = 13;
fig_size = [50, 50, 800, 1000]; %[left bottom width height]
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
zf_tavg = zeros(nrdiag, 1);  zf_shear_tavg = zf_tavg;
zf_shear_polar_tavg = zf_tavg;
perp_drive_tavg = zeros(nrdiag, 1);  para_drive_tavg = perp_drive_tavg;
perp_drive_ex_tavg = perp_drive_tavg;
para_drive_ex_tavg = para_drive_tavg;
shear_diff_tavg = zeros(nrdiag-4, 1);  shear_damp_tavg = zeros(nrdiag, 1);

n_diag = end_diag - start_diag + 1;

for idiag = start_diag: end_diag
	disp(['getting zf drives: step ', num2str(idiag), ' of ', ...
		num2str(end_diag)]);
	data_file = ['dat', sprintf('%4.4d', idiag)];
	load(data_file);
	sdata();

	Tep = field2pol(Te);  denp = field2pol(den);
	phip = field2pol(phi);  jzp = field2pol(jz);
	[vErp, vEthtp] = vec2pol(vEx, vEy);

	% zonal flow
	zf = squeeze(zonal_average(vEthtp(nr_start:nr_end, :, :)));
	zf_tavg = zf_tavg + zf;

	% radial derivative of zonal flow
	zf_shear = 1/(2*dr) .* squeeze( ...
		zonal_average(vEthtp(nr_start+1: nr_end+1, :, :)) ...
		- zonal_average(vEthtp(nr_start-1: nr_end-1, :, :)));
	zf_shear_tavg = zf_shear_tavg + zf_shear;
	
	% the additional term in the laplacian operator is added to zf shear
	r_tmp = reshape(r(nr_start-1: nr_end+1), nrdiag+2, 1);
	zf_shear_polar = 1/(2*dr) ./r_tmp(2:end-1) .* squeeze( ...
		r_tmp(3:end) .* zonal_average(vEthtp(nr_start+1: nr_end+1, :, :)) ...
		- r_tmp(1:end-2) .* zonal_average(vEthtp(nr_start-1: nr_end-1, :, :)));
	zf_shear_polar_tavg = zf_shear_polar_tavg + zf_shear_polar;
	tmp1 = 1/(2*dr)*r_tmp(3:end-2) .* (zf_shear_polar(3:end) - zf_shear_polar(1:end-2));
	r_tmp = r_tmp(4:end-3);
	difr = dif_perp_in * (r_tmp<rdif) + dif_perp_out * (r_tmp>=rdif);
	shear_diff = difr./(2*dr)./r_tmp.*(tmp1(3:end) - tmp1(1:end-2));
	shear_diff_tavg = shear_diff_tavg + shear_diff;
	shear_damp = -momentum_damp * zf_shear_polar;
	shear_damp_tavg = shear_damp_tavg + shear_damp;

	
	% perpendicular shear drive
	Renolds_stress = squeeze(zonal_average( ...
		zonal_pert(vErp(nr_start-2: nr_end+2, :, :)) ...
		.* zonal_pert(vEthtp(nr_start-2: nr_end+2, :, :))));
	perp_drive = -radius/(4*dr^2) * (Renolds_stress(5:end) ...
		- 2*Renolds_stress(3:end-2) + Renolds_stress(1:end-4));
%	perp_drive = -radius/(dr^2) * (Renolds_stress(3:end) ...
%		- 2*Renolds_stress(2:end-1) + Renolds_stress(1:end-2));
	perp_drive_tavg = perp_drive_tavg + perp_drive;

	% parallel shear drive
	tmp1 = zonal_pert(phip) ./ Tep;
	tmp1 = tmp1(nr_start: nr_end, :, :);
	tmp2 = zonal_pert(denp) ./ repmat(zonal_average(denp), 1, ...
		size(denp, 2), size(denp, 3));
	tmp2 = tmp2(nr_start: nr_end, :, :);
	tmp1 = tmp1 - tmp2;
	tmp1(:, :, 2:end-1) = (radius/height)^2 / dz^2 * (tmp1(:, :, 3:end) ...
		- 2*tmp1(:, :, 2:end-1) + tmp1(:, :, 1:end-2));
	tmp1 = inv_nustar * zonal_average(Tep(nr_start: nr_end, :, :)) ...
        .* zonal_average(tmp1 .* tmp2);
	para_drive = squeeze(tmp1);
	para_drive_tavg = para_drive_tavg + para_drive;

	% rigorous perpendicular drive
	r_tmp = repmat(reshape(r(nr_start-2: nr_end+2), nrdiag+4, 1), 1, ntht, nz+2);
	tmp1 = r_tmp .* zonal_pert(vEthtp(nr_start-2: nr_end+2, :, :));
	tmp1 = 1/(2*dr) * (tmp1(3:end, :, :) - tmp1(1:end-2, :, :));
	tmp1 = 1./r_tmp(2:end-1, :, :) .* tmp1;
	tmp1 = squeeze(zonal_average( ...
		zonal_pert(vErp(nr_start-1: nr_end+1, :, :)) .* tmp1));
	perp_drive_ex = -radius/(2*dr) * (tmp1(3:end) - tmp1(1:end-2));
	perp_drive_ex_tavg = perp_drive_ex_tavg + perp_drive_ex;

	% rigorous parallel drive
	tmp1 = zeros(nrdiag, ntht, nz+2);
	tmp1(:, :, 2:end-1) = (radius/height) / (2*dz) ...
		* (zonal_pert(jzp(nr_start: nr_end, :, 3:end)) ...
	   	- zonal_pert(jzp(nr_start: nr_end, :, 1:end-2)));
	para_drive_ex = squeeze(zonal_average(tmp1./denp(nr_start: nr_end, :, :)));
	para_drive_ex_tavg = para_drive_ex_tavg + para_drive_ex;
end
zf_tavg = zf_tavg / n_diag;
zf_shear_tavg = zf_shear_tavg / n_diag;
zf_shear_polar_tavg = zf_shear_polar_tavg / n_diag;
shear_diff_tavg = shear_diff_tavg / n_diag;
shear_damp_tavg = shear_damp_tavg / n_diag;
perp_drive_tavg = perp_drive_tavg / n_diag;
para_drive_tavg = para_drive_tavg / n_diag;
perp_drive_ex_tavg = perp_drive_ex_tavg / n_diag;
para_drive_ex_tavg = para_drive_ex_tavg / n_diag;

%%
close;
figure('name', 'ZF drives');
set(gcf, 'Position', fig_size);
r_axis = rX(nr_start: nr_end);
normalize_factor = cs0/(t0*rhos0);

subplot(3,1,1);
plot(r_axis, cs0*zf_tavg, 'LineWidth', line_width);
title('zonal flow');
ylabel('$$\left<v_{E,\theta}\right>/\left(\rm{cm/s}\right)$$', ...
	'interpreter', 'latex');
set(gca, 'fontSize', font_size);

subplot(3,1,2);
plot(r_axis, cs0/t0*zf_shear_tavg, 'r-', 'LineWidth', line_width);  hold on;
plot(r_axis, cs0/t0*zf_shear_polar_tavg, 'b-', 'LineWidth', line_width);
title('zonal flow shear');
ylabel('shear$$/\rm{s^{-1}}$$', 'interpreter', 'latex');
set(gca, 'fontSize', font_size);
lgd = legend('$$\frac{\partial}{\partial r}\left<v_{E,\theta}\right>$$', ...
	['$$\frac{1}{r}\frac{\partial}{\partial r}\left(r\left<v_{E,\theta}', ...
	'\right>\right)$$'], 'Location', 'NorthEast');
set(lgd, 'interpreter', 'latex');
set(lgd, 'Box', 'off');
set(lgd, 'FontSize', legend_font_size);

subplot(3,1,3);
plot(r_axis, cs0/(t0*rhos0)*perp_drive_tavg, 'r-', 'LineWidth', line_width);
hold on;
plot(r_axis, cs0/(t0*rhos0)*para_drive_tavg, 'b-', 'LineWidth', line_width);
hold on;
plot(r_axis, cs0/(t0*rhos0)*perp_drive_ex_tavg, 'r+', 'LineWidth', line_width);
hold on;
plot(r_axis, cs0/(t0*rhos0)*para_drive_ex_tavg, 'b+', 'LineWidth', line_width);
hold on;
plot(r_axis(3:end-2), cs0/(t0*rhos0)*shear_diff_tavg, 'c-', 'LineWidth', line_width);
hold on;
plot(r_axis, cs0/(t0*rhos0)*shear_damp_tavg, 'm-', 'LineWidth', line_width);
hold on;
plot(r_axis, cs0/(t0*rhos0)*(perp_drive_ex_tavg + para_drive_ex_tavg), 'g-', ...
	'LineWidth', line_width);
hold on;
plot(r_axis(3:end-2), cs0/(t0*rhos0)*(perp_drive_ex_tavg(3:end-2) ...
	+ para_drive_ex_tavg(3:end-2) + shear_diff_tavg + shear_damp_tavg(3:end-2)), ...
	'k-', 'LineWidth', line_width);
xlabel('r/cm');
ylabel('shear\ drive$$/\rm{s^{-2}}$$', 'interpreter', 'latex');
title('drives of zonal flow shear');
set(gca, 'fontSize', font_size);

%lgd = legend( ...
%	['$$-\frac{\partial^{2}}{\partial r^{2}}\left<\delta v_{E,r}', ...
%	'\delta v_{E,\theta}\right>$$'], ...
%	['$$\frac{\Omega_{i}}{m_{e}\nu_{e,i}}\left<T_{e}\right>\left<', ...
%	'\frac{\delta n}{\left<n\right>}\frac{\partial^2}{\partial z^{2}}', ...
%	'\left(\frac{\delta\phi}{T_{e}}-\frac{\delta n}{\left<n\right>}', ...
%	'\right)\right>$$'], 'accurate perp drive', ...
%	'accurate para drive', 'sum of accurate drives', 'Location', 'NorthEast');
lgd = legend('$$F_{\perp}$$', '$$F_{\parallel}$$', ...
	'$$F_{\perp}^{\dagger}$$', '$$F_{\parallel}^{\dagger}$$', ...
	['$$D_{\perp}\frac{1}{r}\frac{\partial}{\partial r}', ...
	'\left(r\frac{\partial}{\partial r}S^{\dagger}\right)$$'], ...
	'$$-\gamma S^{\dagger}$$', '$$\tilde{n}$$ dirve', 'Location', 'NorthEast');
set(lgd, 'NumColumns', 2);
set(lgd, 'Box', 'off');
set(lgd, 'interpreter', 'latex');
set(lgd, 'FontSize', legend_font_size);

print(gcf, '-dpng', ['t=', num2str(start_diag), '-', num2str(end_diag), ...
	'flow_shear_drives.png']);
print(gcf, '-depsc', ['t=', num2str(start_diag), '-', num2str(end_diag), ...
	'flow_shear_drives.eps']);
