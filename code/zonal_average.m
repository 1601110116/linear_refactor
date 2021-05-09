function h = zonal_average(g)
% This function returns the zonal average of field f. 
% f must be in polar grid, and includes z-boundaries
h = mean(mean(g(:, 1:end-1, 2:end-1), 3), 2);
