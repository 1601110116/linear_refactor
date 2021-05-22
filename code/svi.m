function svi (f, fi)
% This function advances vi and applies its boundary conditions
global vi vi_aux calc Te den result 

result = f*vi_aux(2:end-1, 2:end-1, 2:end-1) + fi*vi(2:end-1, 2:end-1, 2:end-1) ...
	+ calc .* ( -Te(2:end-1, 2:end-1, 2:end-1) ./ den(2:end-1, 2:end-1, 2:end-1) ...
	.* ddz(den) - convect(vi) + viscose(vi) + drag(vi));
[nan_ix, nan_iy, nan_iz] = ind2sub(size(result), find(isnan(result)));
[inf_ix, inf_iy, inf_iz] = ind2sub(size(result), find(isinf(result)));
if ~(isempty(nan_ix) && isempty(inf_ix))
    disp('svi: NaN or inf found, simulation paused for debugging');
    pause;
end
if f < 0.6
	vi_aux = vi;
end
vi(2:end-1, 2:end-1, 2:end-1) = result;
vi = zbcs(vi);
vi = lowpass_n(vi);
