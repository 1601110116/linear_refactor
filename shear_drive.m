% This script draws time-averaged profiles of zonal flow drives 
%  versus time1
% Copy this script to the folder of the data and then run it

clear;  close all;
global den Te pe vi jz ve phi vEx vEy dt inv_nustar

load('parameters.mat');
addpath(code_path);
last_file = get_last_file('data');
last_diag = str2num(last_file(end-7: end-4));
%%
%---input---
start_diag = 1226;  % the figure starts from start_diag+1 and ends at end_diag-1
end_diag = 1246; %last_diag;  
r_start = 0.4;  % cm
r_end = 6;  % cm
font_size = 15;
legend_font_size = 14;
no_size = 16;
fig_size = [50, 50, 700, 900]; %[left bottom width height]
view_angle = [45, 15];  % [horizontal rotation, vertical elevation]
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
nustar_zonal_drive_t = zf_shear_t;
den_acc_drive_t = zf_shear_t;
ZF_t = zf_shear_t;
RS_t = zf_shear_t;

for idiag = start_diag: end_diag
	disp(['getting zf drives: step ', num2str(idiag), ' of ', ...
		num2str(end_diag)]);
	data_file = fullfile('data', ['dat', sprintf('%4.4d', idiag)]);
	load(data_file);
	sdata();

	denp = field2pol(den);
	jzp = field2pol(jz);
	Tep = field2pol(Te);
	phip = field2pol(phi);
	[vErp, vEthtp] = vec2pol(vEx, vEy);

	% zonal flow
	ZF = zonal_average(vEthtp(nr_start: nr_end, :, :));
	ZF_t(:, idiag-start_diag+1) = ZF;

	% Reynolds' stress
	RS = zonal_average(vEthtp(nr_start: nr_end, :, :) ...
		.* vErp(nr_start: nr_end, :, :));
	RS_t(:, idiag-start_diag+1) = RS;

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

	% parallel drive neglecting parallel Te/nustar perturbations
	if local_nustar
		inv_nustar_local = inv_nustar * Tep(nr_start: nr_end, :, :).^1.5 ...
			./ denp(nr_start: nr_end, :, :);
	else
		inv_nustar_local = inv_nustar * ones(nrdiag, ntht, nz+2);
    end
%	inv_nustar_zonal = squeeze(zonal_average(inv_nustar_local));
%	pep = denp .* Tep;
%	tmp1 = zeros(nrdiag, ntht, nz+2);  tmp2 = tmp1;
%	tmp1(:, :, 2:end-1) = radius/height/(2*dz)*(pep(nr_start: nr_end, :, 3:end) ...
%		- pep(nr_start: nr_end, : ,1:end-2));
%	tmp2(:, :, 2:end-1) = radius/height/(2*dz)*denp(nr_start: nr_end, :, 2:end-1) ...
%		.* (phip(nr_start: nr_end, :, 3:end) - phip(nr_start: nr_end, :, 1:end-2));
%	jz1 = inv_nustar_zonal .* (tmp1 - tmp2);
%	jz1 = zbcs(jz1);
%	tmp3 = zeros(nrdiag, ntht, nz+2);
%	tmp3(:, :, 2:end-1) = radius/height/(2*dz)*(jz1(:, :, 3:end) - jz1(:, :, 1:end-2));
%	nustar_zonal_drive = squeeze(zonal_average(tmp3./denp(nr_start: nr_end, :, :)));
%	nustar_zonal_drive_t(:, idiag-start_diag+1) = nustar_zonal_drive;
	tmp1 = zeros(nrdiag, ntht, nz+2);  tmp2 = tmp1;
	tmp1(:, :, 2:end-1) = radius/height/(2*dz) ...
		* (denp(nr_start: nr_end, :, 3:end) - denp(nr_start: nr_end, :, 1:end-2));
	tmp1 = zbcs(tmp1);
	tmp2(:, :, 2:end-1) = radius/height/(2*dz) ...
		* (phip(nr_start: nr_end, :, 3:end) - phip(nr_start: nr_end, :, 1:end-2));
	tmp2 = zbcs(tmp2);
	tmp2 = tmp2 .* denp(nr_start: nr_end, :, :) ./ Tep(nr_start: nr_end, :, :);
	tmp1 = tmp1 - tmp2;
	tmp2(:, :, 2:end-1) = radius/height/(2*dz) * (tmp1(:, :, 3:end) - tmp1(:, :, 1:end-2));
	tmp2 = zbcs(tmp2);
	t_nu_zonal = zonal_average(Tep(nr_start: nr_end, :, :) .* inv_nustar_local);
	nustar_zonal_drive = t_nu_zonal .* zonal_average(tmp2./denp(nr_start: nr_end, :, :));
	nustar_zonal_drive_t(:, idiag-start_diag+1) = nustar_zonal_drive;

%	pep = denp .* Tep;
%	tmp1 = zeros(nrdiag, ntht, nz+2);  tmp2 = tmp1;
%	tmp1(:, :, 2:end-1) = radius/height/(2*dz)*(pep(nr_start: nr_end, :, 3:end) ...
%		- pep(nr_start: nr_end, : ,1:end-2));
%	tmp2(:, :, 2:end-1) = radius/height/(2*dz)*denp(nr_start: nr_end, :, 2:end-1) ...
%		.* (phip(nr_start: nr_end, :, 3:end) - phip(nr_start: nr_end, :, 1:end-2));
%	jz1 = inv_nustar_local .* (tmp1 - tmp2);
%	jz1 = zbcs(jz1);
%	tmp3 = zeros(nrdiag, ntht, nz+2);
%	tmp3(:, :, 2:end-1) = radius/height/(2*dz)*(jz1(:, :, 3:end) - jz1(:, :, 1:end-2));
%	nustar_zonal_drive = squeeze(zonal_average(tmp3./denp(nr_start: nr_end, :, :)));
%	para_drive_ex_t(:, idiag-start_diag+1) = nustar_zonal_drive;
	

	% further neglect parallel phi and Te perturbations
	gradparn = zeros(nrdiag, ntht, nz+2);
	gradparn(:, :, 2:end-1) = radius/height/(2*dz) ...
		* (denp(nr_start: nr_end, :, 3:end) - denp(nr_start: nr_end, :, 1:end-2));
	gradparn = zbcs(gradparn);
	gradparn_norm = gradparn ./ denp(nr_start: nr_end, :, :);
	gradparn_norm_zonal = squeeze(zonal_average(gradparn_norm.^2));
	Te_inv_nustar_zonal = squeeze(zonal_average(Tep(nr_start: nr_end, :, :) ...
		.* inv_nustar_local));
	den_acc_drive = Te_inv_nustar_zonal .* gradparn_norm_zonal;
	den_acc_drive_t(:, idiag-start_diag+1) = den_acc_drive;
end
%%
close all;
%dtS_minus_sum = 1/(2*dt*nt_per_diagnose) ...
%	* (zf_shear_polar_t(3:end-2, 3:end) - zf_shear_polar_t(3:end-2, 1:end-2)) - sum_t(:, 2:end-1);
%dtS = 1/(2*dt*nt_per_diagnose) ...
%	* (zf_shear_polar_t(3:end-2, 3:end) - zf_shear_polar_t(3:end-2, 1:end-2));
r_axis = rX(nr_start: nr_end);
normalize_factor = cs0/(t0*rhos0);
tX = t0 * dt * nt_per_diagnose * ((start_diag+1): (end_diag-1));

%figure;
%surf(1e3*tX, r_axis(3:end-2), normalize_factor*sum_t(:,2:end-1));

figure('name', 'ZF drives');
set(gcf, 'Position', fig_size);

subplot(2,1,1);
s1 = surf(1e3*tX, r_axis, cs0/rhos0*zf_shear_t(:, 2:end-1));  hold on;
s2 = surf(1e3*tX, r_axis, cs0/rhos0*zf_shear_polar_t(:, 2:end-1));
set(s1, 'FaceAlpha', 0.2, 'EdgeColor', 'g', 'FaceColor', 'g');
set(s2, 'FaceAlpha', 0.2, 'EdgeColor', 'r', 'FaceColor', 'r');
%title('zonal flow shear');
xlabel('t/ms');
ylabel('r/cm');
zlabel('shear$$\rm \left(s^{-1}\right)$$', 'interpreter', 'latex');
set(gca, 'XLim', [1e3*tX(1), 1e3*tX(end)]);
set(gca, 'YLim', [r_axis(1), r_axis(end)]);
set(gca, 'fontSize', font_size);
%lgd = legend('$$\frac{\partial}{\partial r}\left<v_{E,\theta}\right>$$', ...
%	['$$\frac{1}{r}\frac{\partial}{\partial r}\left(r\left<v_{E,\theta}', ...
%	'\right>\right)$$'], 'Location', 'NorthEast');
lgd = legend('$$S$$', '$$S^{\dagger}$$');
set(lgd, 'interpreter', 'latex');
set(lgd, 'FontSize', legend_font_size);
set(lgd, 'Location', 'NorthWest');
view(view_angle);
tmp = get(gca, 'Position');
axis_width = tmp(3);
tmp = get(lgd, 'Position');
lgd_left = tmp(1);
xlim = get(gca, 'XLim');  ylim = get(gca, 'YLim'); zlim = get(gca, 'ZLim');
txt = text(xlim(2), ylim(2), zlim(end), '(a)', ...
	'HorizontalAlignment', 'right', 'VerticalAlignment', 'top'); 
set(txt, 'fontSize', no_size);

subplot(2,1,2);
s3 = surf(1e3*tX, r_axis, normalize_factor*para_drive_ex_t(:, 2:end-1));  hold on
s4 = surf(1e3*tX, r_axis, normalize_factor*nustar_zonal_drive_t(:, 2:end-1));  hold on;
s5 = surf(1e3*tX, r_axis, normalize_factor*den_acc_drive_t(:, 2:end-1));  hold on;
s6 = surf(1e3*tX, r_axis, normalize_factor*perp_drive_ex_t(:, 2:end-1)); 
set(s3, 'FaceAlpha', 0.2, 'EdgeColor', 'k', 'FaceColor', 'k');
set(s4, 'FaceAlpha', 0.2, 'EdgeColor', 'b', 'FaceColor', 'b');
set(s5, 'FaceAlpha', 0.2, 'EdgeColor', 'g', 'FaceColor', 'g');
set(s6, 'FaceAlpha', 0.2, 'EdgeColor', 'r', 'FaceColor', 'r');
xlabel('t/ms');
ylabel('r/cm');
zlabel('shear\ drive$$\rm \left(s^{-2}\right)$$', 'interpreter', 'latex');
set(gca, 'XLim', [1e3*tX(1), 1e3*tX(end)]);
set(gca, 'YLim', [r_axis(1), r_axis(end)]);
view(view_angle);
set(gca, 'fontSize', font_size);
lgd = legend('$$F_{\parallel}^{\dagger}$$', '$$F_{\parallel}$$', '$$G_{\parallel}$$', ...
	'$$F_{\perp}^{\dagger}$$');
set(lgd, 'interpreter', 'latex');
set(lgd, 'FontSize', legend_font_size);
set(lgd, 'Location', 'NorthWest');
xlim = get(gca, 'XLim');  ylim = get(gca, 'YLim'); zlim = get(gca, 'ZLim');
txt = text(xlim(2), ylim(2), zlim(end), '(b)', ...
	'HorizontalAlignment', 'right', 'VerticalAlignment', 'top'); 
set(txt, 'fontSize', no_size);



fig_name = fullfile('figs', ...
	['t=', num2str(start_diag), '-', num2str(end_diag), 'shear_drives']);
print(gcf, '-dpng', [fig_name, '.png']);
print(gcf, '-depsc', [fig_name, '.eps']);
