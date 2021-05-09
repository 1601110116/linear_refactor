% This script show the result same as vrvz.m, but in a movie 
%  rather than pictures

clear;  close all;
global den Te pe vi jz ve phi vEx vEy vdex vdey dt inv_nustar ...
   calc ddx	
load('parameters.mat');
addpath(code_path);
last_file = get_last_file('./');
last_diagnose = str2num(last_file(end-7: end-4));
save('parameters.mat', 'last_diagnose', '-append');

%---input---
rdiag = 3.2;  % cm
zdiag = 8;
start_diag = 1008;
end_diag = 1045;
font_size = 20;
fig_size = [50, 50, 1200, 500];
%-----------

build_grid;
rX = rhos0*r;
nrdiag = min(find(rX>rdiag));
jz = zeros(nx+2, nx+2, nz+2);
ve = jz;  vEx = jz;  vEy = jz;  vdex = jz;  vdey = jz; 
generate_constants(height, radius, dx, dz, nx, nz, dt, ...
	rdif, dif_perp_in, dif_perp_out, dif_z_in, dif_z_out, ...
	rconduct, conduct_z_in, conduct_z_out, viscosity, ...
	den_damp, momentum_damp, Te_tanhsrc_max, Te_tanhsrc_radius, ...
	Te_tanhsrc_incline, Te_gausssrc_magnitude, Te_gausssrc_sigma, ...
	den_tanhsrc_max, den_tanhsrc_radius, den_tanhsrc_incline, ...
	den_gausssrc_magnitude, den_gausssrc_sigma);
ndiag = end_diag - start_diag + 1;
writerObj = VideoWriter('vrvz.mp4', 'Motion JPEG AVI');
writerObj.FrameRate = 10;
open(writerObj);

for idiag = start_diag: end_diag
	close;
	disp(['visualizing vr vz: step ', num2str(idiag), ' of ', ...
		num2str(end_diag)]);
	data_file = ['dat', sprintf('%4.4d', idiag)];
	load(data_file);
	sdata();
	tX = idiag * dt * nt_per_diagnose * t0;
	vip = field2pol(vi);
	[vErp, vEthtp] = vec2pol(vEx, vEy);
	vEr_pert = get_pert(vErp);
	vi_pert = get_pert(vip);
	fig = figure('Visible', 'off');
	set(fig, 'Position', fig_size);
	subplot(1,2,1)
	plot(cs0*vi_pert(nrdiag, :, zdiag), cs0*vEr_pert(nrdiag, :, zdiag), 'b+');
	xlim([-0.8e4, 0.8e4]);  ylim([-2.8e4, 2.8e4]);
	xlabel('$$v_{\parallel i}/\left(cm/s\right)$$', 'interpreter', 'latex');
	ylabel('$$v_{E,r}/\left(cm/s\right)$$', 'interpreter', 'latex');
    title(['r = ', num2str(rdiag),' cm'])
	set(gca, 'fontSize', font_size);
    vi_mean = get_average(vip);
	subplot(1,2,2)
	plot(rhos0*r, cs0*squeeze(vi_mean(:, 1, zdiag)));
	ylim([-3.6e3, 3.6e3]);
    xlabel('$$r/cm$$', 'interpreter', 'latex');
    ylabel('$$v_{\parallel i}\left(cm/s\right)$$', 'interpreter', 'latex');
    title('$$v_{\parallel i}$$', 'interpreter', 'latex');
	set(gca, 'fontSize', font_size);
	writeVideo(writerObj, getframe(fig));
end	
close(writerObj);
