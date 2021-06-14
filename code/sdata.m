function sdata ()
% This function solves other quantities to be saved other than the 
%  time-advancing quantities.
global Te vi jz ve phi vEx vEy inv_nustar ddx dt calc pe vdex vdey den mu ...
	ln_lambda t0 local_nustar

% ddz() multiplies a redundant dt, which is eliminated by the dt in the denominator
% jz at z boundaries are initially 0 and are never changed
% jz in this code denotes jz/Te
pe = den .* Te;
pe = lowpass_n(pe);
if local_nustar
	inv_nustar_dt = inv_nustar/dt * Te(2:end-1, 2:end-1, 2:end-1).^1.5 ...
		./ den(2:end-1, 2:end-1, 2:end-1);
	jz(2:end-1, 2:end-1, 2:end-1) = inv_nustar_dt .* calc ...
		.* (ddz(pe) - den(2:end-1, 2:end-1, 2:end-1) .* ddz(phi));
else
	jz(2:end-1, 2:end-1, 2:end-1) = inv_nustar/dt .* calc ...
		.* (ddz(pe) - den(2:end-1, 2:end-1, 2:end-1) .* ddz(phi));
end
jz = zbcs(jz);
jz = lowpass_n(jz);
%jz(2:end-1, 2:end-1, 2:end-1) = inv_nustar/dt * calc ...
%	.* (- ddz(phi));
ve = vi - jz./den;
ve = lowpass_n(ve);
vEx(2:end-1, 2:end-1, 2:end-1) = calc .* ddx .* (phi(2:end-1, 1:end-2, 2:end-1) - phi(2:end-1, 3:end, 2:end-1));
vEy(2:end-1, 2:end-1, 2:end-1) = calc .* ddx .* (phi(3:end, 2:end-1, 2:end-1) - phi(1:end-2, 2:end-1, 2:end-1));
vEx = zbcs(vEx);
vEy = zbcs(vEy);
vEx = lowpass_n(vEx);
vEy = lowpass_n(vEy);
%vdex(2:end-1, 2:end-1, 2:end-1) = calc .* ddx .* (pe(2:end-1, 3:end, 2:end-1) ...
%	- pe(2:end-1, 1:end-2, 2:end-1));
%vdey(2:end-1, 2:end-1, 2:end-1) = calc .* ddx .* (pe(1:end-2, 2:end-1, 2:end-1) ...
%	- pe(3:end, 2:end-1, 2:end-1));
