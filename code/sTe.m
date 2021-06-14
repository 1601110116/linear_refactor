function sTe (f, fi)
% This function advances Te and applies its boundary conditions
global Te Te_aux calc src_Te result init_uniform ve

result = f*Te_aux(2:end-1, 2:end-1, 2:end-1) + fi*Te(2:end-1, 2:end-1, 2:end-1) ...
	+ calc .* (-0.6667*Te(2:end-1, 2:end-1, 2:end-1).*ddz(ve) - convect(Te) ...
	+ conduct(Te) + recombine(Te)) + src_Te;
[nan_ix, nan_iy, nan_iz] = ind2sub(size(result), find(isnan(result)));
[inf_ix, inf_iy, inf_iz] = ind2sub(size(result), find(isinf(result)));
if ~(isempty(nan_ix) && isempty(inf_ix))
    disp('sTe: NaN or inf found, simulation paused for debugging');
    pause;
end
if f < 0.6
	Te_aux = Te;
end
Te(2:end-1, 2:end-1, 2:end-1) = result;
Te = lowpass_n(Te);
Te(Te < init_uniform) = init_uniform;
Te = zbcs(Te);
