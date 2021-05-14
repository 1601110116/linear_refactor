clear all;  close all;

load('parameters.mat');
addpath(code_path);

%---input---
idiag = 71;
zdiag = 12;
l_max = 30;
%-----------

data_file = fullfile('data', ['dat', sprintf('%4.4d', idiag)]);
load(data_file);
den = den - filter_z(den, 0);

fig = figure();

build_grid;
denp = field2pol(den);
denp_pert = get_pert(denp(:, :, zdiag));
subplot(1,2,1);
pcolor(yp(:,:,zdiag) * rhos0, xp(:,:,zdiag) * rhos0, denp_pert);
colormap jet;  colorbar;  shading interp;
xlabel('y (cm)');
ylabel('x (cm)');

build_grid_2d;
den_m0 = filter_m(den(:, :, zdiag), 0, l_max);
den_pert = den(:, :, zdiag) - den_m0;
subplot(1,2,2);
pcolor(xX, xX, den_pert);
colormap jet;  colorbar;  shading interp;
xlabel('y (cm)');
ylabel('x (cm)');


