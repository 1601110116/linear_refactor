% This script reads input parameters for simulation, 
%  and saves parameters for further analysis.

global last_diagnose x z xX zX cs0 rhos0 t0 dx dz data_path ...
	inv_nustar ln_lambda max_difxdt min_difxdt cm second dt
if mod(nx, 2) ~= 0
	error('the input parameter nx has to be an even number');
end
if simulate_mode == 1
	if exist('../in_*.m', 'file')
		error('You are in a further/ directory. Simulate_mode 1 unavailable.');
	end
end
if simulate_mode == 1 || simulate_mode == 3
	last_diagnose = 0;
	if ~exist('further', 'dir')
		mkdir('further');
	end
	copyfile('in_*.m', 'further/');
	
	% calculate normalizing units
	omega_c0 = 9.58e3*B0/mu;  % ion gyrofrequency in Hz
	cs0 = 9.79e5*(Tref/mu)^0.5;  % ion sound speed in cm/s
	rhos0 = cs0/omega_c0;  % ion gyroradius in cm
	t0 = radius/cs0;  % t0 in second
	disp(['Ion sound speed is ', num2str(cs0), ' cm/s']);
	disp(['Gyroradius is ', num2str(rhos0), ' cm']);
	disp(['t0 = ', num2str(t0), ' s']);
	disp(['Time step is ', num2str(dt), ' s']);
	cm = 1/rhos0;  % centimeter in ion gyroradius
	second = 1/t0;  % second in t0;

	% Coulomb logarithm
	ln_lambda = 22.36 + 1.5*log(Tref) - 0.5*log(denref);

	% Classical diffusion
	if dif_mode == 1 || dif_mode == 2
		dif_perp_in = 1.652e-5 * denref * (Tref^-0.5) * (B0^-2) * ln_lambda;
		dif_perp_out = dif_perp_in;
		disp(['For denref and Tref: dif_perp = ', num2str(dif_perp_in), 'cm^2/s']);
	end
	if dif_mode == 1
		% classical diffusion with den=0.5*denref and Te=0.5*Tref
		dif_perp_in = sqrt(0.5) * dif_perp_in;
		dif_perp_out = sqrt(0.5) * dif_perp_out;
	end

	% inverseof (me/mi)\nu_e, nu_e is e-i collision frequency
	nustar = 2.906e-6 * ln_lambda * denref * Tref^-1.5 * t0/(mu*1836);
	inv_nustar = 1 / nustar;
	disp(['inv_nustar = ', num2str(inv_nustar)]);

	% nomalize
	radius = radius * cm;
	x_max = x_max * cm;
	height = height * cm;
	Te_tanhsrc_radius = Te_tanhsrc_radius * cm;
	Te_tanhsrc_incline = Te_tanhsrc_incline * cm;
	Te_gausssrc_sigma = Te_gausssrc_sigma * cm;
	den_tanhsrc_radius = den_tanhsrc_radius * cm;
	den_tanhsrc_incline = den_tanhsrc_incline * cm;
	den_gausssrc_sigma = den_gausssrc_sigma * cm;
	rdif = rdif * cm;
	dif_perp_in = dif_perp_in * cm^2 / second;
	dif_perp_out = dif_perp_out * cm^2 / second;
	max_difperp = max_difperp * cm^2 / second;
	min_difperp = min_difperp * cm^2 / second;
	dif_z_in = dif_z_in / ( (height*rhos0)^2/(radius*rhos0/cs0) );
	dif_z_out = dif_z_out / ( (height*rhos0)^2/(radius*rhos0/cs0) );
	rconduct = rconduct * cm;
	conduct_z_in = conduct_z_in / ( (height*rhos0)^2/(radius*rhos0/cs0) );
	conduct_z_out = conduct_z_out / ( (height*rhos0)^2/(radius*rhos0/cs0) );
	viscosity = viscosity * cm^2 / second;
	den_damp = den_damp / second;
	momentum_damp = momentum_damp / second;
	dt = dt * second;


	% generate grid
	dx = 2 * x_max / (nx-1);
	dz = 1 / (nz-1);

	% coordinates of the grids, including boundaries
	x = linspace(-x_max-dx, x_max+dx, nx+2);
	z = linspace(-dz, 1+dz, nz+2);
	% 'X' here means eXperiment
	xX = x * rhos0;
	zX = z * height * rhos0;

	% simulation protection
	max_difxdt = max_difperp * dt / dx^2;
	min_difxdt = min_difperp * dt / dx^2;

	save parameters.mat

elseif simulate_mode == 2
	tmp1 = ndiagnose;
	tmp2 = visual;
	tmp3 = enable_parallel;
	load parameters.mat
	ndiagnose = tmp1;
	visual = tmp2;
	enable_parallel = tmp3;
	simulate_mode = 2;
	last_file = get_last_file(data_path);
	last_diagnose = str2num(last_file(end-7: end-3));
end


