function generate_constants (height, radius, dx, dz, nx, nz, dt, ...
		rdif, dif_perp_in, dif_perp_out, dif_z_in, dif_z_out, ...
		rconduct, conduct_z_in, conduct_z_out, viscosity, ...
		den_damp, momentum_damp, Te_tanhsrc_max, Te_tanhsrc_radius, ...
		Te_tanhsrc_incline, Te_gausssrc_magnitude, Te_gausssrc_sigma, ...
		den_tanhsrc_max, den_tanhsrc_radius, den_tanhsrc_incline, ...
		den_gausssrc_magnitude, den_gausssrc_sigma)
% This function generates all constants used by the simulation
%  in this program ddz means the operator d/dz
%  the generated constants can be a scalar, or a scalar field 
%   excluding boundaries.
global ddzdt ddx convxdt difxdt difzdt src_Te src_den x z xX fig_path ...
	zX calc den_dampdt momentum_dampdt conductzdt viscosityxdt outside

x3d = repmat(reshape(x(2: end-1), [], 1, 1), 1, nx, nz);
y3d = repmat(reshape(x(2: end-1), 1, [], 1), nx, 1, nz);
r3d = sqrt(x3d.^2 + y3d.^2);
calc = r3d < radius;
x3d1 = repmat(reshape(x, [], 1, 1), 1, nx+2, nz+2);
y3d1 = repmat(reshape(x, 1, [], 1), nx+2, 1, nz+2);
r3d1 = sqrt(x3d1.^2 + y3d1.^2);
outside = r3d1 >= radius;


% Diffusion
difxdt = dif_perp_in * (r3d<rdif) + dif_perp_out * (r3d>=rdif);
difzdt = dif_z_in * (r3d<rdif) + dif_z_out * (r3d>=rdif);
difxdt = difxdt .* (dt/dx^2);
difzdt = difzdt * (dt/dz^2);

% Conduction
conductzdt = conduct_z_in * (r3d<rconduct) ...
	+ conduct_z_out * (r3d>=rconduct);
conductzdt = conductzdt .* (dt/dz^2);

% Viscosity
viscosityxdt = viscosity * (dt/dx^2);

% Decay
den_dampdt = den_damp * dt;
momentum_dampdt = momentum_damp * dt;

% Source
src_Te = calc .* Te_tanhsrc_max .* 0.5 .* (1 - ...
	tanh( (r3d-Te_tanhsrc_radius)./Te_tanhsrc_incline ) );
src_Te = src_Te + calc .* Te_gausssrc_magnitude .* ...
	normpdf(r3d, 0, Te_gausssrc_sigma);
src_den = calc .* den_tanhsrc_max .* 0.5 .* (1 - ...
	tanh( (r3d-den_tanhsrc_radius)./den_tanhsrc_incline ) );
src_den = src_den + calc .* den_gausssrc_magnitude .* ...
	normpdf(r3d, 0, den_gausssrc_sigma);
src_Te = src_Te * dt;
src_den = src_den * dt;
%Te_src_zdecay = repmat(reshape(exp(-(1-z(2:end-1))./Te_src_zdecay_len), ...
%	1, 1, []), nx, nx, 1);
%den_src_zdecay = repmat(reshape(exp(-(1-z(2:end-1))./den_src_zdecay_len), ...
%	1, 1, []), nx, nx, 1);
%src_Te = src_Te .* Te_src_zdecay;
%src_den = src_den .* den_src_zdecay;

% Derivatives
ddzdt = radius/height / (2*dz) * dt;
ddx = 1 / (2*dx);

% Convection
convxdt = radius / (2*dx) * dt;


close all;
figure;
subplot(2,2,1);  pcolor(xX(2:end-1), xX(2:end-1), src_Te(:,:,2)./dt);
colormap jet;  colorbar;  shading interp;
xlabel('y/cm'); ylabel('x/cm');
title('$$src\_pe/\left(T_{0}n_{0}/t_{0}\right)$$', 'interpreter', 'latex');

subplot(2,2,2);  pcolor(xX(2:end-1), xX(2:end-1), src_den(:,:,2)./dt);
colormap jet;  colorbar;  shading interp;
xlabel('y/cm');  ylabel('x/cm');
title('$$src\_den/\left(n_{0}/t_{0}\right)$$', 'interpreter', 'latex');

subplot(2,2,3);  pcolor(xX(2:end-1), xX(2:end-1), difxdt(:,:,2)./(dt/dx^2));
colormap jet;  colorbar;  shading interp;
xlabel('y/cm');  ylabel('x/cm');
title('$$dif\_perp/\left(\rho_{s}^{2}/t_{0}\right)$$', 'interpreter', 'latex');

fig_file = fullfile(fig_path, 'Sources and diffusion coefficients');
print(gcf, '-dpng', fig_file);
