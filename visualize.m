% This script draws the fields using all dat****.mat files in its directory
% Copy this script to the folder of the data and then run it

clear;  close all;
global den Te vi jz ve phi vEx vEy vdex vdey dt inv_nustar

load('parameters.mat');
addpath(code_path);
last_file = get_last_file('data/');
last_diagnose = str2num(last_file(end-7: end-4));
save('parameters.mat', 'last_diagnose', '-append');

%---input---
ydiag = length(xX) / 2;
zdiag = ceil(length(zX)/2);
x_tick = -10: 5: 10;
start_diag = 299;
end_diag = 300;
fig_position = [50, 50, 1500, 1000];
%-----------


jz = zeros(nx+2, nx+2, nz+2);
ve = jz;  vEx = jz;  vEy = jz;  vdex = jz;  vdey = jz;
generate_constants(height, radius, dx, dz, nx, nz, dt, ... 
	    rdif, dif_perp_in, dif_perp_out, dif_z_in, dif_z_out, ... 
	    rconduct, conduct_z_in, conduct_z_out, viscosity, ... 
	    den_damp, momentum_damp, Te_tanhsrc_max, Te_tanhsrc_radius, ... 
	    Te_tanhsrc_incline, Te_gausssrc_magnitude, Te_gausssrc_sigma, ... 
	    den_tanhsrc_max, den_tanhsrc_radius, den_tanhsrc_incline, ... 
	    den_gausssrc_magnitude, den_gausssrc_sigma);

for idiag = start_diag: end_diag
	close;
	disp(['visualizing fields: step ', num2str(idiag), ' of', ...
		num2str(end_diag)]);
	data_file = fullfile('data', ['dat', sprintf('%4.4d', idiag)]);
	load(data_file);
	sdata();
	tX = idiag * dt * nt_per_diagnose * t0;
	figure('name', pwd, 'Visible', 'off');  
    %set(gcf, 'Position', get(0, 'ScreenSize'));
    set(gcf, 'Position', fig_position);

	subplot(4,4,1);  pcolor(xX, xX, denref*den(:, :, zdiag));
	colormap jet;  colorbar;  shading interp;
	title('$$n/cm^-3$$', 'interpreter', 'latex');
	xlabel('y/cm');  ylabel('x/cm');
	set(gca, 'XTick', x_tick);
	set(gca, 'YTick', x_tick);

	subplot(4,4,2);  pcolor(xX, xX, Tref*Te(:,:,zdiag));
	colormap jet;  colorbar;  shading interp;
	title('$$T_{e}/eV$$', 'interpreter', 'latex');
	xlabel('y/cm');  ylabel('x/cm');
	set(gca, 'XTick', x_tick);
	set(gca, 'YTick', x_tick);

	subplot(4,4,3);  pcolor(xX, xX, cs0*vi(:,:,zdiag));
	colormap jet;  colorbar;  shading interp;
	title('$$v_{\parallel i}/\left(cm/s\right)$$', 'interpreter', 'latex');
	xlabel('y/cm');  ylabel('x/cm');
	set(gca, 'XTick', x_tick);
	set(gca, 'YTick', x_tick);

	subplot(4,4,4);  pcolor(xX, xX, cs0/rhos0*w(:,:,zdiag));
	colormap jet;  colorbar;  shading interp;
	title('$$w/s^{-1}$$', 'interpreter', 'latex');
	xlabel('y/cm');  ylabel('x/cm');
	set(gca, 'XTick', x_tick);
	set(gca, 'YTick', x_tick);

	subplot(4,4,5);  pcolor(xX, xX, cs0*jz(:,:,zdiag));
	colormap jet;  colorbar;  shading interp;
	title('$$j_{\parallel}/\left(cm/s\right)$$', 'interpreter', 'latex');
	xlabel('y/cm');  ylabel('x/cm');
	set(gca, 'XTick', x_tick);
	set(gca, 'YTick', x_tick);

	subplot(4,4,6);  pcolor(xX, xX, cs0*ve(:,:,zdiag));
	colormap jet;  colorbar;  shading interp;
	title('$$v_{\parallel e}/\left(cm/s\right)$$', 'interpreter', 'latex');
	xlabel('y/cm');  ylabel('x/cm');
	set(gca, 'XTick', x_tick);
	set(gca, 'YTick', x_tick);

	subplot(4,4,7);  pcolor(xX, xX, Tref*phi(:,:,zdiag));
	colormap jet;  colorbar;  shading interp;
	title('$$\phi/V$$', 'interpreter', 'latex');
	xlabel('y/cm');  ylabel('x/cm');
	set(gca, 'XTick', x_tick);
	set(gca, 'YTick', x_tick);


	subplot(4,4,9);  pcolor(zX, xX, squeeze(denref*den(:,ydiag,:)));
	colormap jet;  colorbar;  shading interp;
	title('$$n/cm^-3$$', 'interpreter', 'latex');
	xlabel('z/cm');  ylabel('x/cm');
%	set(gca, 'XTick', x_tick);
	set(gca, 'YTick', x_tick);

	subplot(4,4,10);  pcolor(zX, xX, squeeze(Tref*Te(:,ydiag,:)));
	colormap jet;  colorbar;  shading interp;
	title('$$T_{e}/eV$$', 'interpreter', 'latex');
	xlabel('z/cm');  ylabel('x/cm');
	set(gca, 'XTick', x_tick);
	set(gca, 'YTick', x_tick);

	subplot(4,4,11);  pcolor(zX, xX, cs0*squeeze(vi(:,ydiag,:)));
	colormap jet;  colorbar;  shading interp;
	title('$$v_{\parallel i}/\left(cm/s\right)$$', 'interpreter', 'latex');
	xlabel('z/cm');  ylabel('x/cm');
	set(gca, 'XTick', x_tick);
	set(gca, 'YTick', x_tick);

	subplot(4,4,12);  pcolor(zX, xX, cs0/rhos0*squeeze(w(:,ydiag,:)));
	colormap jet;  colorbar;  shading interp;
	title('$$w/s^{-1}$$', 'interpreter', 'latex');
	xlabel('z/cm');  ylabel('x/cm');
	set(gca, 'XTick', x_tick);
	set(gca, 'YTick', x_tick);

	subplot(4,4,13);  pcolor(zX, xX, cs0*squeeze(jz(:,ydiag,:)));
	colormap jet;  colorbar;  shading interp;
	title('$$j_{\parallel}/\left(cm/s\right)$$', 'interpreter', 'latex');
	xlabel('z/cm');  ylabel('x/cm');
	set(gca, 'XTick', x_tick);
	set(gca, 'YTick', x_tick);

	subplot(4,4,14);  pcolor(zX, xX, cs0*squeeze(ve(:,ydiag,:)));
	colormap jet;  colorbar;  shading interp;
	title('$$v_{\parallel e}/\left(cm/s\right)$$', 'interpreter', 'latex');
	xlabel('z/cm');  ylabel('x/cm');
	set(gca, 'XTick', x_tick);
	set(gca, 'YTick', x_tick);

	subplot(4,4,15);  pcolor(zX, xX, Tref*squeeze(phi(:,ydiag,:)));
	colormap jet;  colorbar;  shading interp;
	title('$$\phi/V$$', 'interpreter', 'latex');
	xlabel('z/cm');  ylabel('x/cm');
	set(gca, 'XTick', x_tick);
	set(gca, 'YTick', x_tick);

	subplot(4,4,16);  xlabel(['time = ', num2str(tX), ...
		sprintf(' s\nidiag = '), num2str(idiag)]);
    
	fig_file = fullfile('diagnose', ['zdiag=', num2str(zdiag), 'fields', ...
		sprintf('%4.4d', idiag), '.png']);
    print(gcf, '-dpng', fig_file)
end
