% This script gets the zonal_average of den, vi, w, Te and phi

clear;  close all;

load('parameters.mat');
addpath(code_path);

%%
%---input---
idiag = 1212;  %1190;
font_size = 15;
fig_position = [50, 50, 1200, 800];
line_width = 2;
%-----------

build_grid;
data_file = fullfile('data', ['dat', sprintf('%4.4d', idiag)]);
load(data_file)

denp = field2pol(den);
vip  = field2pol(vi);
wp   = field2pol(w);
Tep  = field2pol(Te);
phip = field2pol(phi);

den_profile = squeeze(zonal_average(denp));
vi_profile  = squeeze(zonal_average(vip));
w_profile   = squeeze(zonal_average(wp));
Te_profile  = squeeze(zonal_average(Tep));
phi_profile = squeeze(zonal_average(phip));

rX = rhos0 * r;
tX = idiag * dt * nt_per_diagnose * t0;
close;
fig = figure;
set(fig, 'Position', fig_position);
suptitle(['tind=', sprintf('%4.4d', idiag), ';  t=', num2str(tX), 'ms']);

subplot(2, 3, 1);  plot(rX, denref * den_profile, 'lineWidth', line_width);
xlabel('$$r$$ (cm)', 'interpreter', 'latex');
ylabel('$$n\ \rm\left(cm^{-3}\right)$$', 'interpreter', 'latex');
set(gca, 'fontSize', font_size)

subplot(2, 3, 2);  plot(rX, cs0 * vi_profile, 'lineWidth', line_width);
xlabel('$$r$$ (cm)', 'interpreter', 'latex');
ylabel('$$vi\ \rm\left(cm/s\right)$$', 'interpreter', 'latex');
set(gca, 'fontSize', font_size)

subplot(2, 3, 3);  plot(rX, cs0/rhos0 * w_profile, 'lineWidth', line_width);
xlabel('$$r$$ (cm)', 'interpreter', 'latex');
ylabel('$$w\ \rm\left(s^{-1}\right)$$', 'interpreter', 'latex');
set(gca, 'fontSize', font_size)

subplot(2, 3, 4);  plot(rX, Tref * Te_profile, 'lineWidth', line_width);
xlabel('$$r$$ (cm)', 'interpreter', 'latex');
ylabel('$$T_{e}\ \rm\left(eV\right)$$', 'interpreter', 'latex');
set(gca, 'fontSize', font_size)

subplot(2, 3, 5);  plot(rX, Tref * phi_profile, 'lineWidth', line_width);
xlabel('$$r$$ (cm)', 'interpreter', 'latex');
ylabel('$$\phi\ \rm\left(V\right)$$', 'interpreter', 'latex');
set(gca, 'fontSize', font_size)

profiles_path = fullfile(pwd, 'extracted_profiles');
if ~exist(profiles_path, 'dir')
	mkdir('extracted_profiles');
end

prof_name = fullfile(profiles_path, ['profiles_t=', sprintf('%4.4d', idiag)]);
print(fig, '-dpng', [prof_name, '.png']);
save(prof_name, 'rX', 'den_profile', 'vi_profile', 'w_profile', 'Te_profile', 'phi_profile');
