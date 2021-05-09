% This script makes a movie to show the evolution of full density, electric potential, density 
%  profile and EXB velocity profile

clear;
global den phi Te

load('parameters.mat');
addpath(code_path);
last_file = get_last_file('data');
last_diag = str2num(last_file(end-7: end-4));

%---input---
start_diag = 707;
end_diag = last_diag;
zdiag = ceil(length(zX)/2);  % the x-y plane for contours
r_start = 0.4;  %cm
r_end = 8;
x_tick_contour = -10: 5: 10;
x_tick_plot = 0:2:8;
font_size = 20;
fig_position = [50,50,1500,600];
frame_rate = 8;
line_width = 2;
%-----------

build_grid;
rX = rhos0 * r;
nr_start = find(rX >= r_start, 1, 'first');
if nr_start < 2
    error(['r_start must be greater than ', num2str(dr*rhos0), ...
        ' for the current grid']);
end
nr_end = find(rX <= r_end, 1, 'last');
nrdiag = nr_end - nr_start + 1;
x_lim_plot = [x_tick_plot(1), x_tick_plot(end)];

[xX2d, yX2d] = ndgrid(xX(2:end-1), xX(2:end-1));  % excludes boundaries
rX2d = sqrt(xX2d.^2 + yX2d.^2);
out2d = rX2d > (radius*rhos0);  % get the indices of the outside grid points
r_axis = rX(nr_start: nr_end);
vx = zeros(nx+2, nx+2, nz+2);
vy = zeros(nx+2, nx+2, nz+2);

movie_name = fullfile('figs', ['tinymovie_t=', num2str(start_diag), '-', num2str(end_diag), '.avi']);
writerObj = VideoWriter(movie_name, 'Motion JPEG AVI');
writerObj.FrameRate = frame_rate;
open(writerObj);

for idiag = start_diag: end_diag
	disp(['making tiny movie, step ', num2str(idiag), ' of ', num2str(end_diag)]);
	data_file = fullfile('data', ['dat', sprintf('%4.4d', idiag)]);
	load(data_file, 'den', 'phi', 'Te');

	den_disp = den(2:end-1, 2:end-1, zdiag);
	den_disp(out2d) = nan;
    Te_disp = Te(2:end-1, 2:end-1, zdiag);
    Te_disp(out2d) = nan;
	phi_disp = phi(2:end-1, 2:end-1, zdiag);
	phi_disp(out2d) = nan;
    vx(2:end-1, 2:end-1, 2:end-1) = 1/(2*dx) * ...
        (phi(2:end-1, 1:end-2, 2:end-1) - phi(2:end-1, 3:end, 2:end-1));
    vy(2:end-1, 2:end-1, 2:end-1) = 1/(2*dx) * ...
        (phi(3:end, 2:end-1, 2:end-1) - phi(1:end-2, 2:end-1, 2:end-1));
	[vr, vtht] = vec2pol(vx, vy);
	denp = field2pol(den);
    Tep = field2pol(Te);
	den_profile = zonal_average(denp(nr_start: nr_end, :, :));
    Te_profile = zonal_average(Tep(nr_start: nr_end, :, :));
	zf_profile = zonal_average(vtht(nr_start: nr_end, :, :));
	tX_ms = idiag * dt * nt_per_diagnose * t0 * 1e3;

	fig = figure('Visible', 'off');  set(gcf, 'Position', fig_position);

	subplot(2,3,1);  pcolor(yX2d, xX2d, denref*den_disp);
	colormap jet;  colorbar;  shading interp;
	title('$$n/\rm cm^{-3}$$', 'interpreter', 'latex');
	xlabel('$$y$$/cm', 'interpreter', 'latex');
	ylabel('$$x$$/cm', 'interpreter', 'latex');
    set(gca, 'XTick', x_tick_contour);
    set(gca, 'YTick', x_tick_contour);
	set(gca, 'FontSize', font_size);
    
    subplot(2,3,2);  pcolor(yX2d, xX2d, Tref*Te_disp);
	colormap jet;  colorbar;  shading interp;
	title('$$T_{e}/\rm eV$$', 'interpreter', 'latex');
	xlabel('$$y$$/cm', 'interpreter', 'latex');
	ylabel('$$x$$/cm', 'interpreter', 'latex');
    set(gca, 'XTick', x_tick_contour);
    set(gca, 'YTick', x_tick_contour);
	set(gca, 'FontSize', font_size);

	subplot(2,3,3);  pcolor(yX2d, xX2d, Tref*phi_disp);
	colormap jet;  colorbar;  shading interp;
	title('$$\phi/\rm V$$', 'interpreter', 'latex');
	xlabel('$$y$$/cm', 'interpreter', 'latex');
	ylabel('$$x$$/cm', 'interpreter', 'latex');
    set(gca, 'XTick', x_tick_contour);
    set(gca, 'YTick', x_tick_contour);
	set(gca, 'FontSize', font_size);

	subplot(2,3,4);  plot(r_axis, denref*den_profile, 'b-', 'LineWidth', line_width);
	xlabel('$$r$$/cm', 'interpreter', 'latex');
	ylabel('$$\left<n\right>/\rm cm^{-3}$$', 'interpreter', 'latex');
    set(gca, 'XTick', x_tick_plot);
    set(gca, 'XLim',x_lim_plot);
	set(gca, 'FontSize', font_size);
    
    subplot(2,3,5);  plot(r_axis, Tref*Te_profile, 'b-', 'LineWidth', line_width);
	xlabel('$$r$$/cm', 'interpreter', 'latex');
	ylabel('$$\left<T_{e}\right>/\rm eV$$', 'interpreter', 'latex');
    set(gca, 'XTick', x_tick_plot);
    set(gca, 'XLim',x_lim_plot);
	set(gca, 'FontSize', font_size);
    

	subplot(2,3,6);  plot(r_axis, 1e-2*cs0*zf_profile, 'b-', 'LineWidth', line_width);
	xlabel(['$$r$$/cm', newline, 't = ', num2str(tX_ms), ' ms'], ...
		'interpreter', 'latex');
	ylabel('$$\left<v_{E,\theta}\right>/\rm \left(m/s\right)$$', 'interpreter', 'latex');
    set(gca, 'XTick', x_tick_plot);
    set(gca, 'XLim',x_lim_plot);
	set(gca, 'FontSize', font_size);

	pause(0.2);
	writeVideo(writerObj, getframe(fig));
	close;
end
close(writerObj);
