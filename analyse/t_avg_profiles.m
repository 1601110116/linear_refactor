% This script draws the time averaged profiles using data files between
%  start_diagnose and end_diagnose in its directory
% Copy this script to the folder of the data and then run it

clear;  close all;
global den Te pe vi jz ve phi vEx vEy vdex vdey dt inv_nustar ...
    calc ddx

load('parameters.mat');

addpath(code_path);
last_file = get_last_file('./');
last_diagnose = str2num(last_file(end-7: end-4));
save('parameters.mat', 'last_diagnose', '-append');

%---input---
start_diagnose = 491;
end_diagnose = 540;  %last_diagnose
zdiag = 8;
font_size = 15;
fig_size = [50, 50, 1000, 720];
%-----------

build_grid;
jz = zeros(nx+2, nx+2, nz+2);
ve = jz;  vEx = jz;  vEy = jz;  vdex = jz;  vdey = jz;
generate_constants(height, radius, dx, dz, nx, nz, dt, ... 
    rdif, dif_perp_in, dif_perp_out, dif_z_in, dif_z_out, ... 
    rconduct, conduct_z_in, conduct_z_out, viscosity, ... 
    den_damp, momentum_damp, Te_tanhsrc_max, Te_tanhsrc_radius, ... 
    Te_tanhsrc_incline, Te_gausssrc_magnitude, Te_gausssrc_sigma, ... 
    den_tanhsrc_max, den_tanhsrc_radius, den_tanhsrc_incline, ... 
    den_gausssrc_magnitude, den_gausssrc_sigma);
Te_tavg = zeros(1, nr);
den_tavg = Te_tavg;  vi_tavg = Te_tavg;  phi_tavg = Te_tavg;
vEtht_tavg = Te_tavg;  vdetht_tavg = Te_tavg;
radial_particle_flux_tavg = Te_tavg;
axial_Reynolds_stress_tavg = Te_tavg;
azimuthal_Reynolds_stress_tavg = Te_tavg;

n_diag = end_diagnose - start_diagnose + 1;

for idiagnose = start_diagnose: end_diagnose
	close;
	disp(['getting profiles: step ', num2str(idiagnose), ' of', ...
		num2str(end_diagnose)]);
	data_file = ['dat', sprintf('%4.4d', idiagnose)];
	load(data_file);
	sdata();
	tX = idiagnose * dt * nt_per_diagnose * t0;

	Tep = field2pol(Te);  denp = field2pol(den);  vip = field2pol(vi);
	phip = field2pol(phi);
	pe = Te .* den;
	vdex(2:end-1, 2:end-1, 2:end-1) = calc .* ddx .* (pe(2:end-1, 3:end, 2:end-1) ...
		- pe(2:end-1, 1:end-2, 2:end-1)) ./ den(2:end-1, 2:end-1, 2:end-1);
	vdey(2:end-1, 2:end-1, 2:end-1) = calc .* ddx .* (pe(1:end-2, 2:end-1, 2:end-1) ...
		- pe(3:end, 2:end-1, 2:end-1)) ./ den(2:end-1, 2:end-1, 2:end-1);
	[vErp, vEthtp] = vec2pol(vEx, vEy);
	[vderp, vdethtp] = vec2pol(vdex, vdey);

	Te_mean = get_average(Tep);
	den_mean = get_average(denp);
	vi_mean = get_average(vip);
	phi_mean = get_average(phip);
	vEr_mean = get_average(vErp);
	vEtht_mean = get_average(vEthtp);
	vder_mean = get_average(vderp);
	vdetht_mean = get_average(vdethtp);

	radial_particle_flux = get_average(get_pert(vErp) .* get_pert(denp));
	axial_Reynolds_stress = get_average(get_pert(vErp) .* get_pert(vip));
	azimuthal_Reynolds_stress = get_average(get_pert(vErp) .* get_pert(vEthtp));

	Te_tavg = Te_tavg + squeeze(Te_mean(:, 1, zdiag));
	den_tavg = den_tavg + squeeze(den_mean(:, 1, zdiag));
	vi_tavg = vi_tavg + squeeze(vi_mean(:, 1, zdiag));
	phi_tavg = phi_tavg + squeeze(phi_mean(:, 1, zdiag));
	vEtht_tavg = vEtht_tavg + squeeze(vEtht_mean(:, 1, zdiag));
	vdetht_tavg = vdetht_tavg + squeeze(vdetht_mean(:, 1, zdiag));
	radial_particle_flux_tavg = radial_particle_flux_tavg + squeeze(radial_particle_flux(:, 1, zdiag));
	axial_Reynolds_stress_tavg = axial_Reynolds_stress_tavg + squeeze(axial_Reynolds_stress(:, 1, zdiag));
	azimuthal_Reynolds_stress_tavg = azimuthal_Reynolds_stress_tavg + squeeze(azimuthal_Reynolds_stress(:, 1, zdiag));
end
Te_tavg = Te_tavg / n_diag;
den_tavg = den_tavg / n_diag;
vi_tavg = vi_tavg / n_diag;
phi_tavg = phi_tavg / n_diag;
vEtht_tavg = vEtht_tavg / n_diag;
vdetht_tavg = vdetht_tavg / n_diag;
radial_particle_flux_tavg = radial_particle_flux_tavg / n_diag;
axial_Reynolds_stress_tavg = axial_Reynolds_stress_tavg / n_diag;
azimuthal_Reynolds_stress_tavg = azimuthal_Reynolds_stress_tavg / n_diag;

fig = figure;  set(gcf, 'Position', get(0, 'ScreenSize'));

subplot(3,3,1);  plot(rhos0*r, Tref*Te_tavg);
xlabel('r/cm', 'interpreter', 'latex');  
ylabel('$$T_{e}/eV$$', 'interpreter', 'latex');
set(gca, 'fontSize', font_size);

subplot(3,3,2);  plot(rhos0*r, denref*den_tavg);
xlabel('r/cm', 'interpreter', 'latex');  
ylabel('$$n/cm^{-3}$$', 'interpreter', 'latex');
set(gca, 'fontSize', font_size);

subplot(3,3,3);  plot(rhos0*r, cs0*vi_tavg);
xlabel('r/cm', 'interpreter', 'latex');  
ylabel('$$v_{\parallel i}/\left(cm/s\right)$$', 'interpreter', 'latex');
set(gca, 'fontSize', font_size);

subplot(3,3,4);  plot(rhos0*r, Tref*phi_tavg);
xlabel('r/cm', 'interpreter', 'latex');  
ylabel('$$\phi/V$$', 'interpreter', 'latex');
set(gca, 'fontSize', font_size);

subplot(3,3,5);  plot(rhos0*r, denref*cs0*radial_particle_flux_tavg);
xlabel('r/cm', 'interpreter', 'latex');  
ylabel('$$\left<\tilde{n}\tilde{v_{E,r}}\right>/\left(cm^{-2}s^{-1}\right)$$', ...
	'interpreter', 'latex');
set(gca, 'fontSize', font_size);

subplot(3,3,6);  plot(rhos0*r, cs0*vEtht_tavg);
xlabel('r/cm', 'interpreter', 'latex');  
ylabel('$$v_{E,\theta}/\left(cm/s\right)$$', 'interpreter', 'latex');
set(gca, 'fontSize', font_size);

subplot(3,3,7);  plot(rhos0*r, cs0*vdetht_tavg);
xlabel('r/cm', 'interpreter', 'latex');  
ylabel('$$v_{de,\theta}/\left(cm/s\right)$$', 'interpreter', 'latex');
set(gca, 'fontSize', font_size);

subplot(3,3,8);  plot(rhos0*r, cs0^2*axial_Reynolds_stress_tavg);
xlabel('r/cm', 'interpreter', 'latex');  
ylabel('$$\left<\tilde{v_{E,r}}\tilde{v_{\parallel i}}\right>/\left(cm^{2}/s^{2}\right)$$', ...
	'interpreter', 'latex');
set(gca, 'fontSize', font_size);

subplot(3,3,9);  plot(rhos0*r, cs0^2*azimuthal_Reynolds_stress_tavg);
xlabel('r/cm', 'interpreter', 'latex');  
ylabel('$$\left<\tilde{v_{E,r}}\tilde{v_{E,\theta}}\right>/\left(cm^{2}/s^{2}\right)$$', 'interpreter', 'latex');
set(gca, 'fontSize', font_size);

set(fig, 'Position', fig_size);

print(gcf, '-dpng', ['zdiag=', num2str(zdiag), '_t=', ...
	num2str(start_diagnose), '-', num2str(end_diagnose), '_avg_profiles.png']);


