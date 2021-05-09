clear;  close all;

load('parameters.mat');
addpath(code_path);

%%
%---input---
font_size = 28;
fig_position = [50, 50, 600, 600];
start_diag = 1189;
end_diag = 1225;
r_diag = 3;  % cm
m_max = 10;
line_width = 2;
font_size = 15;
%-----------

build_grid;
rX = rhos0 * r;
ir_diag = find(rX > r_diag, 1, 'first');
n_diag = end_diag - start_diag + 1;
I_m_t = zeros(m_max, n_diag);

for idiag = start_diag: end_diag
	disp(['step ', num2str(idiag) ' of ' num2str(end_diag)]);
	data_file = fullfile('data', ['dat', sprintf('%4.4d', idiag)]);
	load(data_file, 'den', 'phi');
	fieldp = field2pol(phi);
	for iz = 2: nz+1
		field_m = squeeze(fft(fieldp(ir_diag, 1: end-1, iz), [], 2));
		I_m = 0.5 * (abs(field_m(2: m_max+1)) / (ntht/2)) .^ 2;
		I_m_t(:, idiag - start_diag + 1) = I_m_t(:, idiag-start_diag+1) + reshape(I_m, m_max, 1);
	end
	I_m_t(:, idiag-start_diag+1) = I_m_t(:, idiag-start_diag+1) ./ nz;
end

%%
close;
m = 1:m_max;
tX = t0 * dt * nt_per_diagnose * (start_diag: end_diag);
tX_ms = tX * 1e3;
fig = figure('name', 'm_components');
line_colors = jet(m_max);
for im = 1: m_max
	plot(tX_ms, log10(Tref^2 * I_m_t(im, :)), 'Color', line_colors(im, :), ...
		'lineWidth', line_width);  hold on;
end
hold off;
grid on;
colormap(line_colors);  cbar = colorbar;
set(cbar.Label, 'String', 'm');
set(gca, 'CLim', [0.5, m_max + 0.5])
xlabel('$$t \mathrm{(ms)}$$', 'interpreter', 'latex');
ylabel('$$\lg\left<\tilde{\phi}_m^2\right>$$', 'interpreter', 'latex');
set(gca, 'fontSize', font_size);

pic_name = fullfile('figs', ['m_components_phi_t=', sprintf('%4.4d', start_diag), '-', ...
	sprintf('%4.4d', end_diag)]);
print(gcf, '-dpng', [pic_name, '.png']);
