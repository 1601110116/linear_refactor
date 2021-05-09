function initialize (simulate_mode, init_uniform, init_perturbation)
% This function gets the inittial values of all simulated quantities
global Te den Te_aux den_aux vi w vi_aux w_aux jz ve ...
	phi vEx vEy last_diagnose nx nz calc vdex vdey data_path
jz = zeros(nx+2, nx+2, nz+2);
ve = jz;  vEx = jz;  vEy = jz;  vdex = jz;  vdey = jz;  phi = jz;

if simulate_mode == 1
	Te = init_uniform .* ones(nx+2, nx+2, nz+2);
	den = init_uniform .* ones(nx+2, nx+2, nz+2);
	Te(2:end-1, 2:end-1, 2:end-1) = Te(2:end-1, 2:end-1, 2:end-1) + ...
		init_perturbation * calc .* (rand(nx, nx, nz) - 0.5);
	den(2:end-1, 2:end-1, 2:end-1) = den(2:end-1, 2:end-1, 2:end-1) + ...
		init_perturbation * calc .* (rand(nx, nx, nz) - 0.5);
% 	x3d = repmat(reshape(x, [], 1, 1), 1, nx+2, nz+2);
% 	y3d = repmat(reshape(x, 1, [], 1), nx+2, 1, nz+2);
% 	r3d = sqrt(x3d.^2 + y3d.^2);   
% 	Te(r3d<Te_tanhsrc_radius) = Te(r3d<Te_tanhsrc_radius) + 0.3;
% 	den(r3d<Te_tanhsrc_radius) = den(r3d<Te_tanhsrc_radius) + 0.3;
	w = zeros(nx+2, nx+2, nz+2);  vi = w;
%	w(2:end-1, 2:end-1, 2:end-1) = init_perturbation * calc .* (rand(nx,nx,nz)-0.5);
	sphi(nx, nz);
	sdata();
elseif simulate_mode == 2
	last_file = get_last_file(data_path);
	disp(['started from file ', last_file]);
	load(last_file)
	sdata();
elseif simulate_mode == 3
	last_file = get_last_file(fullfile('..', 'data'));
	disp(['started from ', last_file]);
	last_diagnose = 0;
	% copy the last data file in the previous directory to current
	% directory as rest.mat and then load it
	copyfile(last_file, 'rest.mat');
	load rest.mat
	sdata();
end
den_aux = den;  Te_aux = Te;
vi_aux = vi;      w_aux = w;

