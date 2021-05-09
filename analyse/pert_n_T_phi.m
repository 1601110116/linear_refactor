% This script draws the perturbed density, temperature and electric potential in its directory
%  'p' means polar and 'c' means cartecian
% Copy this script to the folder of the data and then run it

clear;  close all;
load('parameters.mat');
addpath(code_path);
last_file = get_last_file('./');
last_diagnose = str2num(last_file(end-7: end-4));
save('parameters.mat', 'last_diagnose', '-append');
build_grid;


rdiag = floor(nr/3);
zdiag = 9;


for idiagnose = 266: 266
	close;
	disp(['visualizing perturbed fields: step ', num2str(idiagnose), ' of', ...
		num2str(last_diagnose)]);
	data_file = ['dat', sprintf('%4.4d', idiagnose)];
	load(data_file, 'den', 'Te', 'phi');
	tX = idiagnose * dt * nt_per_diagnose * t0;

	denp = interp3(yc, xc, zc, den, yp, xp, zp);
	Tep = interp3(yc, xc, zc, Te, yp, xp, zp);
	phip = interp3(yc, xc, zc, phi, yp, xp, zp);

%	den_pert = denp - repmat(mean(denp(:,1:end-1,:), 2), 1, ntht, 1);
%	Te_pert = Tep - repmat(mean(Tep(:,1:end-1,:), 2), 1, ntht, 1);
%	phi_pert = phip - repmat(mean(phip(:,1:end-1,:), 2), 1, ntht, 1);
	den_pert = get_pert(denp);
	Te_pert = get_pert(Tep);
	phi_pert = get_pert(phip);

	figure;  set(gcf, 'Position', get(0, 'ScreenSize'));
	subplot(2,3,1);  pcolor(rhos0*yp(:,:,zdiag), rhos0*xp(:,:,zdiag), denref*den_pert(:, :, zdiag));
	colormap jet;  colorbar;  shading interp;
	title('$$\delta n/cm^{-3}$$', 'interpreter', 'latex');
	xlabel('y/cm');  ylabel('x/cm');
	
	subplot(2,3,2);  pcolor(rhos0*yp(:,:,zdiag), rhos0*xp(:,:,zdiag), Tref*Te_pert(:, :, zdiag));
	colormap jet;  colorbar;  shading interp;
	title('$$\delta T_{e}/eV$$', 'interpreter', 'latex');
	xlabel('y/cm');  ylabel('x/cm');

	subplot(2,3,3);  pcolor(rhos0*yp(:,:,zdiag), rhos0*xp(:,:,zdiag), Tref*phi_pert(:, :, zdiag));
	colormap jet;  colorbar;  shading interp;
	title('$$\delta \phi/V$$', 'interpreter', 'latex');
	xlabel('y/cm');  ylabel('x/cm');

	subplot(2,3,4);  pcolor(height*rhos0*squeeze(zp(rdiag, :, :)), squeeze(thtp(rdiag, :, :)), ...
		denref*squeeze(den_pert(rdiag, :, :)));
	colormap jet;  colorbar;  shading interp;
	title('$$\delta n/cm^{-3}$$', 'interpreter', 'latex');
	xlabel('z/cm');  ylabel('$$\theta/rad$$', 'interpreter', 'latex');

	subplot(2,3,5);  pcolor(height*rhos0*squeeze(zp(rdiag, :, :)), squeeze(thtp(rdiag, :, :)), ...
		Tref*squeeze(Te_pert(rdiag, :, :)));
	colormap jet;  colorbar;  shading interp;
	title('$$\delta Te/eV$$', 'interpreter', 'latex');
	xlabel('z/cm');  ylabel('$$\theta/rad$$', 'interpreter', 'latex');

	subplot(2,3,6);  pcolor(height*rhos0*squeeze(zp(rdiag, :, :)), squeeze(thtp(rdiag, :, :)), ...
		Tref*squeeze(phi_pert(rdiag, :, :)));
	colormap jet;  colorbar;  shading interp;
	title('$$\delta \phi/V$$', 'interpreter', 'latex');
	xlabel(['z/cm', sprintf('\ntime = '), num2str(tX), sprintf(' s\nidiagnose = '), num2str(idiagnose)]);
   	ylabel('$$\theta/rad$$', 'interpreter', 'latex');

	print(gcf, '-dpng', ['zdiag=', num2str(zdiag), 'pert_n_T_phi', sprintf('%4.4d', idiagnose), '.png']);
end
