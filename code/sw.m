function sw (f, fi)
% This function advances w and applies its boundary conditions
global w w_aux calc jz result den vi

%result = f*w_aux(2:end-1, 2:end-1, 2:end-1) + fi*w(2:end-1, 2:end-1, 2:end-1) ...
%	+ calc .* ( ddz(jz) ./ den(2:end-1, 2:end-1, 2:end-1) ...
%   	- convect(w) + diffuse(w) + decay(w) );
result = f*w_aux(2:end-1, 2:end-1, 2:end-1) + fi*w(2:end-1, 2:end-1, 2:end-1) ...
	+ calc .* ( ddz(jz) ./ den(2:end-1, 2:end-1, 2:end-1) ...
   	- convect(w) + viscose(w) + drag(w) );
%result = f*w_aux(2:end-1, 2:end-1, 2:end-1) + fi*w(2:end-1, 2:end-1, 2:end-1) ...
%	+ calc .*( ddz(jz) + diffuse(w));
[nan_ix, nan_iy, nan_iz] = ind2sub(size(result), find(isnan(result)));
[inf_ix, inf_iy, inf_iz] = ind2sub(size(result), find(isinf(result)));
if ~(isempty(nan_ix) && isempty(inf_ix))
    disp('sw: NaN or inf found, simulation paused for debugging');
    pause;
end
if f < 0.6
	w_aux = w;
end
w(2:end-1, 2:end-1, 2:end-1) = result;
w = zbcs(w);

