% This script is the input file for the CSDX device at UCSD
% Reference papers:
% Yan Z, Xu M, Diamond P H, et al. Intrinsic rotation from a residual stress at the boundary of a cylindrical laboratory plasma.[J]. Physical Review Letters, 2010, 104(6):065002.


clear; close all;

global dt nt_per_diagnose nx nz visual Tref denref init_uniform zbc_mode B0 mu ...
	local_nustar enable_parallel data_path diagnose_path fig_path dif_mode n_up_to

% Simulate mode selection
%  1: start from the very beginning of the simulation
%  2: continue from the last .mat file in ./  . Nothing should be 
%     changed except ndiagnose and visual
%  3: continue from the last .mat file in ../  . Used to simulate in 
%     different parameters
simulate_mode = 1;

% Automatic path management
if simulate_mode == 1
	if ~exist('code/', 'file')
		error('No code/ found. You might be in a further/ directory. Simulate_mode 1 unavailable.');
	end
	code_path = fullfile(pwd, 'code');
elseif simulate_mode == 2
	if ~exist('parameters.mat', 'file')
		error('No parameters.mat file found. Simulate_mode 2 unavailable.');
	end
	load('parameters.mat', 'code_path', 'data_path');
elseif simulate_mode == 3
	if ~exist('../parameters.mat', 'file')
		error('No parameters.mat file found in ../  ,Simulate_mode 3 unavailable');
	end
	load('../parameters.mat', 'code_path');
else
	error('Invalid simulate_mode');
end
addpath(code_path);

data_path = fullfile(pwd, 'data');
diagnose_path = fullfile(pwd, 'diagnose');
fig_path = fullfile(pwd, 'figs');
if ~exist(data_path, 'dir')
	mkdir('data')
end
if ~exist(diagnose_path, 'dir')
	mkdir('diagnose')
end
if ~exist(fig_path, 'dir')
	mkdir('figs')
end











% Missions
visual = 1;
enable_parallel = 1;  % solve the Poisson equaiton through parallel computing



% Normalizing units
Tref = 3;  % reference temperature in eV
denref = 10e12;  % reference density in cm^-3
B0 = 1e3;  % Bz in Gauss
mu = 40;  % ion mass over proton mass. mu_argon=40


% Time step
dt = 0.004*2.611e-5;  % time step width in second
nt_per_diagnose = 300;
%  total diagnose times you want in current directory, which may not be
%  finished since the program might be interrupted
ndiagnose = 1800;


% Device settings
radius = 10;  % the radius of the vacuum vessel in cm
height = 270;  % the height of the cylinder vessel in cm


% Grid  (grids at x, y, and z boundaries are not included)
%  we use a cubic mesh, with x\in[-x_max, x_max]cm, y\in[-x_max, x_max]cm, 
%  z\in[0,height]cm
x_max = 11;
%  number of grid points in both x and y directions(excludes boundaries)
%  has to be an even number
nx = 200; %146
nz = 22;


% Initial conditions
init_uniform = 0.004;
init_perturbation = 1e-14;


% Boundary conditions
%  x and y boundaries are given by the [calc] matrix
% z boundary conditions: 1 for free, 2 for periodic, 3 for sheath (not ready yet)
zbc_mode = 2;


% Source
%  rotational symmetric soure, lengths in cm
Te_tanhsrc_max = 5.5e-4/0.0028;
Te_tanhsrc_radius = 3;
Te_tanhsrc_incline = 0.55;
Te_gausssrc_magnitude = 0;%0.0005;
Te_gausssrc_sigma = 2.0;

den_tanhsrc_max = 5.5e-4/0.0028;  %5.45e-4
den_tanhsrc_radius = 3;
den_tanhsrc_incline = 0.55;
den_gausssrc_magnitude = 0;
den_gausssrc_sigma = 2.0;


% Diffusion
%  the Coulomb loarithm is always calculated using denref and Tref

%  dif_mode: model of perpendicular diffusion coefficient:
%  1: classical diffusion with den=0.5*denref and Te=0.5*Tref
%  2: updated using classical diffusion with local den and Te
%  3: two layers of constant diffusion coefficients:
%	      dif_perp_in for r<rdif, dif_perp_out for r>rdif
dif_mode = 1;
rdif = 3.2;  % in cm
dif_perp_in = 6e3;  % in cm^2/s
dif_perp_out = 6e3;
% limit dif_perp
max_difperp = 6e3;  % in cm^2/s
min_difperp = 3e2;

%  parallel diffusion coefficient for den, w and vi
dif_z_in = 0; %9.5e6;
dif_z_out = 0; %9.5e6;%2e4;
n_up_to = 3;  % turn off lowpass if <0


% Conduction
%  parallel thermal conductivity \chi
rconduct = 3.2;  % in cm
conduct_z_in = 188.8888e6;  % in cm^2/s
conduct_z_out = 188.8888e6;


% Viscisity (uniform)
% viscosity is the perpendicular diffusion coefficient of vorticity and
%  ion parallel momentum given by ion-ion collisionality
viscosity = 4.8e3;  % in cm^2/s


% Decay (uniform)
%  loss of density by recombination
%  loss of temperature is unclear and is set to be at the same rate of density damping
%  loss of momentum by ion-neutral collsion

den_damp = 5e3;
momentum_damp = 3.4e3;


% Resistivity (by electron-ion collision \nu)
%  local_nustar=0: inv_nustar is calculated using denref and Tref
%  local_nustar=1: inv_nustar is calculated using local den and Te
%   Coulomb logarithm is always calculated using denref and Tref
local_nustar = 1;

main;
