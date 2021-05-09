% This script makes a movie to show the fields between start_diagnose 
%  and end_diagnose in its directory
% Copy this script to the folder of the data and then run it

clear;  close all;
global den Te vi jz ve phi vEx vEy vdex vdey dt inv_nustar

load('parameters.mat');


addpath(code_path);
last_file = get_last_file('./');
last_diagnose = str2num(last_file(end-7: end-4));
save('parameters.mat', 'last_diagnose', '-append');

%---input---
start_diagnose = 1;
end_diagnose = last_diagnose;
ydiag = length(xX) / 2;
zdiag = 9;
font_size = 9.8;
fig_size = [50, 50, 1500, 950];
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
writerObj = VideoWriter(['movie_t=', num2str(start_diagnose), ...
    '-', num2str(end_diagnose), '.mp4'], 'MPEG-4');
writerObj.FrameRate = 8;
open(writerObj);

for idiagnose = start_diagnose: end_diagnose
	disp(['making movie: step ', num2str(idiagnose), ' of', ...
		num2str(end_diagnose)]);
	data_file = ['dat', sprintf('%4.4d', idiagnose)];
	load(data_file);
	sdata();
	tX = idiagnose * dt * nt_per_diagnose * t0;
	fig = figure;  set(gcf, 'Position', fig_size);

	subplot(4,4,1);  pcolor(xX, xX, denref*den(:, :, zdiag));
	colormap jet;  colorbar;  shading interp;  set(gca, 'fontSize', 9);
	title('$$n/cm^-3$$', 'interpreter', 'latex');
	xlabel('y/cm');  ylabel('x/cm');
    set(gca, 'fontSize', font_size);

	subplot(4,4,2);  pcolor(xX, xX, Tref*Te(:,:,zdiag));
	colormap jet;  colorbar;  shading interp;  set(gca, 'fontSize', 9);
	title('$$T_{e}/eV$$', 'interpreter', 'latex');
	xlabel('y/cm');  ylabel('x/cm');
    set(gca, 'fontSize', font_size);


	subplot(4,4,3);  pcolor(xX, xX, cs0*vi(:,:,zdiag));
	colormap jet;  colorbar;  shading interp;  set(gca, 'fontSize', 9);
	title('$$v_{\parallel i}/\left(cm/s\right)$$', 'interpreter', 'latex');
	xlabel('y/cm');  ylabel('x/cm');
    set(gca, 'fontSize', font_size);

	subplot(4,4,4);  pcolor(xX, xX, Tref/rhos0^2*w(:,:,zdiag));
	colormap jet;  colorbar;  shading interp;  set(gca, 'fontSize', 9);
	title('$$w/\left(V*cm^{-2}\right)$$', 'interpreter', 'latex');
	xlabel('y/cm');  ylabel('x/cm');
    set(gca, 'fontSize', font_size);

	subplot(4,4,5);  pcolor(xX, xX, cs0*jz(:,:,zdiag));
	colormap jet;  colorbar;  shading interp;  set(gca, 'fontSize', 9);
	title('$$j_{\parallel}/\left(cm/s\right)$$', 'interpreter', 'latex');
	xlabel('y/cm');  ylabel('x/cm');
    set(gca, 'fontSize', font_size);

	subplot(4,4,6);  pcolor(xX, xX, cs0*ve(:,:,zdiag));
	colormap jet;  colorbar;  shading interp;  set(gca, 'fontSize', 9);
	title('$$v_{\parallel e}/\left(cm/s\right)$$', 'interpreter', 'latex');
	xlabel('y/cm');  ylabel('x/cm');
    set(gca, 'fontSize', font_size);

	subplot(4,4,7);  pcolor(xX, xX, Tref*phi(:,:,zdiag));
	colormap jet;  colorbar;  shading interp;  set(gca, 'fontSize', 9);
	title('$$\phi/V$$', 'interpreter', 'latex');
	xlabel('y/cm');  ylabel('x/cm');
    set(gca, 'fontSize', font_size);


	subplot(4,4,9);  pcolor(zX, xX, squeeze(denref*den(:,ydiag,:)));
	colormap jet;  colorbar;  shading interp;  set(gca, 'fontSize', 9);
	title('$$n/cm^-3$$', 'interpreter', 'latex');
	xlabel('z/cm');  ylabel('x/cm');
    set(gca, 'fontSize', font_size);

	subplot(4,4,10);  pcolor(zX, xX, squeeze(Tref*Te(:,ydiag,:)));
	colormap jet;  colorbar;  shading interp;  set(gca, 'fontSize', 9);
	title('$$T_{e}/eV$$', 'interpreter', 'latex');
	xlabel('z/cm');  ylabel('x/cm');
    set(gca, 'fontSize', font_size);

	subplot(4,4,11);  pcolor(zX, xX, cs0*squeeze(vi(:,ydiag,:)));
	colormap jet;  colorbar;  shading interp;  set(gca, 'fontSize', 9); grid on;
	title('$$v_{\parallel i}/\left(cm/s\right)$$', 'interpreter', 'latex');
	xlabel('z/cm');  ylabel('x/cm');
    set(gca, 'fontSize', font_size);

	subplot(4,4,12);  pcolor(zX, xX, Tref/rhos0^2*squeeze(w(:,ydiag,:)));
	colormap jet;  colorbar;  shading interp;  set(gca, 'fontSize', 9);
	title('$$w/\left(V*cm^{-2}\right)$$', 'interpreter', 'latex');
	xlabel('z/cm');  ylabel('x/cm');
    set(gca, 'fontSize', font_size);

	subplot(4,4,13);  pcolor(zX, xX, cs0*squeeze(jz(:,ydiag,:)));
	colormap jet;  colorbar;  shading interp;  set(gca, 'fontSize', 9);
	title('$$j_{\parallel}/\left(cm/s\right)$$', 'interpreter', 'latex');
	xlabel('z/cm');  ylabel('x/cm');
    set(gca, 'fontSize', font_size);

	subplot(4,4,14);  pcolor(zX, xX, cs0*squeeze(ve(:,ydiag,:)));
	colormap jet;  colorbar;  shading interp;  set(gca, 'fontSize', 9);
	title('$$v_{\parallel e}/\left(cm/s\right)$$', 'interpreter', 'latex');
	xlabel('z/cm');  ylabel('x/cm');
    set(gca, 'fontSize', font_size);

	subplot(4,4,15);  pcolor(zX, xX, Tref*squeeze(phi(:,ydiag,:)));
	colormap jet;  colorbar;  shading interp;  set(gca, 'fontSize', 9);
	title('$$\phi/V$$', 'interpreter', 'latex');
	xlabel('z/cm');  ylabel('x/cm');
    set(gca, 'fontSize', font_size);

	subplot(4,4,16);  xlabel(['time = ', num2str(tX), ...
		sprintf(' s\nidiagnose = '), num2str(idiagnose)]);
    set(gca, 'fontSize', font_size);
    
	pause(0.1);
	writeVideo(writerObj, getframe(fig));
	close;
end
close(writerObj);
