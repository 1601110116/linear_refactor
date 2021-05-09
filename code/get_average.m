function h = get_average(g)
% This function gets the azimuthal average of a given field in polar grid
h = mean(g(:, 1:end-1, :), 2);
