function sden (f, fi)
% This function advances den and applies its boundary conditions
global den den_aux calc src_den result init_uniform ve

result = f*den_aux(2:end-1, 2:end-1, 2:end-1) + fi*den(2:end-1, 2:end-1, 2:end-1) ...
	+ calc .* (-ddz(den.*ve) - convect(den) + diffuse(den) + recombine(den)) ...
	+ src_den;
[nan_ix, nan_iy, nan_iz] = ind2sub(size(result), find(isnan(result)));
[inf_ix, inf_iy, inf_iz] = ind2sub(size(result), find(isinf(result)));
if ~(isempty(nan_ix) && isempty(inf_ix))
	disp('sden: NaN or inf found, simulation paused for debugging');
	pause;
end
if f < 0.6
	den_aux = den;
end
den(2:end-1, 2:end-1, 2:end-1) = result;
den(den < init_uniform) = init_uniform;
den = zbcs(den);
den = lowpass_n(den);
