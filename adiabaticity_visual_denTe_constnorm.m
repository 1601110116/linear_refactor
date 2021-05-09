% This script draws the fields using all dat****.mat files in its directory
% Copy this script to the folder of the data and then run it

clear;  close all;
global den Te phi

%%
load('parameters.mat');
addpath(code_path);
last_file = get_last_file('data');
last_diag = str2num(last_file(end-7: end-4));

%---input---
start_diag = 1035;
end_diag = 1045;
r_diag = 3; %cm
font_size = 25;
fig_position = [50, 50, 1000, 500];
%-----------
%%

build_grid;
rX = rhos0 * r;
nr_diag = find(rX >=r_diag, 1, 'first');
heightX = height * rhos0;
if ~exist(fullfile('figs', 'adiabaticity_diag'), 'dir')
	mkdir(fullfile('figs', 'adiabaticity_diag'));
end

for idiag = start_diag: end_diag
	close;
	disp(['visualizing den and phi at r = ', num2str(r_diag), ' cm. step ', ...
		num2str(idiag), ' of ', num2str(end_diag)]);
	data_file = fullfile('data', ['dat', sprintf('%4.4d', idiag)]);
	load(data_file, 'den', 'Te', 'phi');

	Tep = field2pol(Te);  denp = field2pol(den);  phip = field2pol(phi);
	den_zonal = zonal_average(denp(nr_diag, :, :));
	Te_zonal = zonal_average(Tep(nr_diag, :, :));
	den_pert_norm = zonal_pert(denp(nr_diag, :, :)) ./ zonal_average(denp(nr_diag, :, :));
	phi_pert_norm = zonal_pert(phip(nr_diag, :, :)) ./ zonal_average(Tep(nr_diag, :, :));
	Te_pert_norm = zonal_pert(Tep(nr_diag, :, :)) ./ zonal_average(Tep(nr_diag, : ,:));
	fig = figure;
	set(fig, 'Position', fig_position);

	subplot(1,3,1);
	pcolor(zX, tht, squeeze(den_pert_norm));
	colormap jet;  colorbar;  shading interp;
	title('$$\tilde{n}/n$$', 'interpreter', 'latex');
	xlabel('$$z$$/cm', 'interpreter', 'latex');
	ylabel('$$\theta$$/rad', 'interpreter', 'latex');
	set(gca, 'FontSize', font_size);

	subplot(1,3,2);
	pcolor(zX, tht, squeeze(phi_pert_norm));
	colormap jet;  colorbar;  shading interp;
	title('$$e\tilde{\phi}/T_{e}$$', 'interpreter', 'latex');
	xlabel('$$z$$/cm', 'interpreter', 'latex');
	ylabel('$$\theta$$/rad', 'interpreter', 'latex');
	set(gca, 'FontSize', font_size);

	subplot(1,3,3);
	pcolor(zX, tht, squeeze(Te_pert_norm));
	colormap jet;  colorbar;  shading interp;
	title('$$\tilde{T_e}/T_e$$', 'interpreter', 'latex');
	xlabel('$$z$$/cm', 'interpreter', 'latex');
	ylabel('$$\theta$$/rad', 'interpreter', 'latex');
	set(gca, 'FontSize', font_size);

	fig_name = fullfile('figs', 'adiabaticity_diag', ...
		['r=', num2str(r_diag), 'cm_den_phi_Te_pert_constnorm', ...
		sprintf('%4.4d', idiag), '.png']);
	print(gcf, '-dpng', fig_name);
end
