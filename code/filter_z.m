function f_filtered = filter_z(f, n)
% This function filters the last dimension of a 3D field using fft. 
%
% Input:
%  f: the field to filter. Its last dimension is filtered, usually z.
%     Note: f is assumed to include boundaries, and the parallel
%           boundary condition is periodic
%  n: the harmonic to keep. n=1 means the wavelength is the 
%     parallel size of the device.

f_n = fft(f(:, :, 2:end-1), [], 3);
res_n = zeros(size(f_n));
res_n(:, :, n+1) = f_n(:, :, n+1);
if n > 0
	res_n(:, :, end-(n-1)) = f_n(:, :, end-(n-1));
end
f_filtered = zeros(size(f));
f_filtered(:, :, 2:end-1) = ifft(res_n, [], 3);

f_filtered = zbcs(f_filtered);
%f_filtered(:, :, 1) = f_filtered(:, :, end-1);
%f_filtered(:, :, end) = f_filtered(:, :, 2);
