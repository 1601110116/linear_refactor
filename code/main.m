global last_diagnose

read_parameters;

generate_constants(height, radius, dx, dz, nx, nz, dt, ...
	rdif, dif_perp_in, dif_perp_out, dif_z_in, dif_z_out, ...
	rconduct, conduct_z_in, conduct_z_out, viscosity, ...
	den_damp, momentum_damp, Te_tanhsrc_max, Te_tanhsrc_radius, ...
	Te_tanhsrc_incline, Te_gausssrc_magnitude, Te_gausssrc_sigma, ...
	den_tanhsrc_max, den_tanhsrc_radius, den_tanhsrc_incline, ...
	den_gausssrc_magnitude, den_gausssrc_sigma);

build_Poisson_coefficient_matrix(dx, nx, x_max);

initialize(simulate_mode, init_uniform, init_perturbation);

diagnose(last_diagnose, nt_per_diagnose);
diagnose_start = last_diagnose + 1;

for idiagnose = diagnose_start : ndiagnose
	disp(['step ', num2str(idiagnose), ' of ', num2str(ndiagnose)]);
	for it = 1: nt_per_diagnose
		f = 0.5;  fi = 0.5;  % constants for the time-advancing method
		sterms();
		sden(f, fi);
		sTe(f, fi);
		svi(f, fi);
		sw(f, fi);
		sphi(nx, nz);
		sdata();

		f = 1.0;  fi = 0.0;
		sterms();
		sden(f, fi);
		sTe(f, fi);
		svi(f, fi);
		sw(f, fi);
		sphi(nx, nz);
		sdata();
	end
	diagnose(idiagnose, nt_per_diagnose);
	last_diagnose = idiagnose;
	save('parameters.mat', 'last_diagnose', '-append');
end

