function h = zonal_fluctuation(g)
% This function returns the fluctuating part of field.
% h is g minus zonal averaged g
% f must be in polar grid
h = g - repmat(zonal_average(g), 1, size(g,2), size(g,3));
