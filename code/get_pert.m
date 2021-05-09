function h = get_pert(g)
% This function gets the perturbed field of a given field in polar grid
h = g - repmat(get_average(g), 1, size(g, 2), 1);

