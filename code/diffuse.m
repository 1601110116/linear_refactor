function h = diffuse(g)
% This function calculates the sum of the diffusion terms of g
%  the resulting term h excludes boundaries
global difxdt difzdt den Te max_difxdt min_difxdt dif_mode dt cm second 
persistent consistent_dif

if dif_mode == 2
	consistent_dif = difxdt .* den(2:end-1, 2:end-1, 2:end-1) ...
		.* Te(2:end-1, 2:end-1, 2:end-1).^-0.5;
	consistent_dif(consistent_dif>max_difxdt) = max_difxdt;
	consistent_dif(consistent_dif<min_difxdt) = min_difxdt;
	h = consistent_dif .* (g(3:end, 2:end-1, 2:end-1) - 2*g(2:end-1, 2:end-1, 2:end-1) ...
		+ g(1:end-2, 2:end-1, 2:end-1)) + ...
		consistent_dif .* (g(2:end-1, 3:end, 2:end-1) - 2*g(2:end-1, 2:end-1, 2:end-1) ...
		+ g(2:end-1, 1:end-2, 2:end-1)) + ...
		difzdt .* (g(2:end-1, 2:end-1, 3:end) - 2*g(2:end-1, 2:end-1, 2:end-1) ...
		+ g(2:end-1, 2:end-1, 1:end-2));
else
	h = difxdt .* (g(3:end, 2:end-1, 2:end-1) - 2*g(2:end-1, 2:end-1, 2:end-1) ...
		+ g(1:end-2, 2:end-1, 2:end-1)) + ...
		difxdt .* (g(2:end-1, 3:end, 2:end-1) - 2*g(2:end-1, 2:end-1, 2:end-1) ...
		+ g(2:end-1, 1:end-2, 2:end-1)) + ...
		difzdt .* (g(2:end-1, 2:end-1, 3:end) - 2*g(2:end-1, 2:end-1, 2:end-1) ...
		+ g(2:end-1, 2:end-1, 1:end-2));
end