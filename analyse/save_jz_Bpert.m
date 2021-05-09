% This script calculate the perturbed magnetic field and save to B_pert/

clear;  close all;
global den Te vi jz ve phi vEx vEy vdex vdey dt inv_nustar calc outside ...
    ddx
if ~exist('B_pert')
	mkdir('B_pert');
end

%%
load('parameters.mat');
addpath(code_path);
last_file = get_last_file('./');
last_diag = str2num(last_file(end-7: end-4));

%---input---
start_diag = 0;
end_diag = last_diag;
%-----------
%%

jz = zeros(nx+2, nx+2, nz+2);
ve = jz;  vEx = jz;  vEy = jz;  vdex = jz;  vdey = jz;
Psi = jz;  Bx = jz;  By = jz;
generate_constants(height, radius, dx, dz, nx, nz, dt, ...
        rdif, dif_perp_in, dif_perp_out, dif_z_in, dif_z_out, ...
        rconduct, conduct_z_in, conduct_z_out, viscosity, ...
        den_damp, momentum_damp, Te_tanhsrc_max, Te_tanhsrc_radius, ...
        Te_tanhsrc_incline, Te_gausssrc_magnitude, Te_gausssrc_sigma, ...
        den_tanhsrc_max, den_tanhsrc_radius, den_tanhsrc_incline, ...
        den_gausssrc_magnitude, den_gausssrc_sigma);
tmp1 = 2*speye(nx) - diag(ones(1,nx-1), 1) - diag(ones(1, nx-1), -1);
tmp2 = -(dx^-2) * (kron(speye(nx), tmp1) + kron(tmp1, speye(nx)));
calc_col = reshape(calc(:,:,1), [], 1); 
lap = tmp2(calc_col, calc_col);
beta0 = 4.03e-11 * Tref * denref / B0^2;

for idiag = 0: end_diag
	disp(['calculating B_pert, step', num2str(idiag), ' of ', num2str(end_diag)]);
	data_file = ['dat', sprintf('%4.4d', idiag)];
	load(data_file);
	sdata();
    Psi_par = Psi(2:end-1, 2:end-1, :);
    jz_par = jz;
	parfor iz = 2:nz+1
		tmp_bigcol = beta0/2 * reshape(jz_par(2:nx+1, 2:nx+1, iz), nx*nx, 1);
		tmp_smallcol = tmp_bigcol(calc_col);
		tmp_smallcol = lap \ tmp_smallcol;
		tmp_bigcol(calc_col) = tmp_smallcol;
		Psi_par(:, :, iz) = reshape(tmp_bigcol, nx, nx);
    end
    Psi(2:end-1, 2:end-1, :) = Psi_par;
	Psi(outside) = 0;
	Psi = zbcs(Psi);
	Bx(2:end-1, 2:end-1, 2:end-1) = calc .* ddx .* (Psi(2:end-1, 1:end-2, 2:end-1) - Psi(2:end-1, 3:end, 2:end-1));
	By(2:end-1, 2:end-1, 2:end-1) = calc .* ddx .* (Psi(3:end, 2:end-1, 2:end-1) - Psi(1:end-2, 2:end-1, 2:end-1));
	Bx = zbcs(Bx);
	By = zbcs(By);
	dst_file = ['B_pert/jz_Bpert', sprintf('%4.4d', idiag)];
	save(dst_file, 'jz', 'Psi', 'Bx', 'By');
end
